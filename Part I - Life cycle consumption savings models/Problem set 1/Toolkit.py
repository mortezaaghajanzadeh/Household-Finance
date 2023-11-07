import scipy.stats as st
import scipy as sp
import numpy as np
import numba
from scipy.spatial import KDTree
import matplotlib.pyplot as plt

def dsearchn(value, array):
    kdt = KDTree(value)
    return np.array([int(kdt.query(v)[1]) for v in array])


def discretize_assets_double_exp(amin, amax, n_a):
    # find maximum ubar of uniform grid corresponding to desired maximum amax of asset grid
    ubar = np.log(1 + np.log(1 + amax - amin))
    
    # make uniform grid
    u_grid = np.linspace(0, ubar, n_a)
    
    # double-exponentiate uniform grid and add amin to get grid from amin to amax
    return amin + np.exp(np.exp(u_grid) - 1) - 1


def discretize_assets_single_exp(amin, amax, n_a):
	# lower bound of exponential grid
	u_under = 1

	# upper bound of exponential grid
	u_bar = np.log(1 + amax - amin) + 1

	# linear grid for starter
	A_lin = np.linspace(u_under, u_bar, n_a) 

	# grid for assets
	A = amin + np.exp(A_lin - 1) - 1

	return A


# retirement value function

def retirement(N_a,a_grid, rr, β, γ, t_r, t_w, g_t, λ,vmin):
    Vr = np.zeros((N_a+1, t_r))
    Cr = np.zeros((N_a+1, t_r))
    Xr = np.zeros((N_a+1, t_r))

    Xr[1:,-1:] = a_grid 
    Cr[1:,-1:] = Xr[1:,-1:]
    Vr[1:,-1:] = Cr[1:,-1:]**(1-γ)/(1-γ)
    Vr[0,:] = vmin

    trend_pension = λ * g_t[t_w-1]

    for t in range(t_r-1,0, -1):
        t -= 1
        # I have to think more about the pension income
        pension = trend_pension

        Xp = a_grid * (1 + rr) + pension ## cash-on-hand tomorrow

        Cp = np.interp(Xp,Xr[:,t+1], Cr[:,t+1]) # interpolate consumption

        # Construct tomorrows value and tomorrows derivative wrt assets
        EV = β * np.interp(Xp,Xr[:,t+1], Vr[:,t+1])
        dV = β * Cp ** (-γ) * (1 + rr)

        Cr[1:,t:t+1] = dV ** (-1/γ) # Use FOC to find consumption policy
        Xr[1:,t:t+1] = Cr[1:,t:t+1] + a_grid # Implied cash on hand
        Vr[1:,t:t+1] = (Cp ** (1-γ) - 1 )/(1-γ) + EV
    return Vr, Cr, Xr

# working value function
def working(N_z, N_ω, N_a, a_grid, Z, ω, π, rw, rr, Xr, Cr, Vr, β, γ, t_w, g_t,λ,vmin):

    Vw = np.zeros((N_z * (N_a+1), t_w))
    Cw = np.zeros((N_z * (N_a+1), t_w))
    Xw = np.zeros((N_z * (N_a+1), t_w))

    for i in range(N_z):
        Vw[i * (N_a+1),:] = vmin
    pension = λ * g_t[t_w-1]

    for t in range(t_w, 0, -1):
        if t == t_w: # Last period working
            Xp = a_grid * (1 + rr) + pension
            Cp = np.interp(Xp,Xr[:,0], Cr[:,0])
            ## % Construct tomorrows value and tomorrows derivative wrt assets (no expectation)
            EV = β * np.interp(Xp,Xr[:,0], Vr[:,0])
            dV = β * Cp ** (-γ) * (1 + rr)
            for i in range(N_z):
                index = range(i * (N_a+1) + 1, (i + 1) * (N_a+1))
                Cw[index  ,t-1:t] = dV ** (-1/γ) # Use FOC to find consumption policy
                Xw[index ,t-1:t] =  Cw[index ,t-1:t] + a_grid
                Vw[index ,t-1:t] = (Cp ** (1-γ) - 1 )/(1-γ) + EV
        else:
            Cp = np.zeros((N_a, N_z * N_ω))
            Vp = np.zeros((N_a, N_z * N_ω))
            for i in range(N_z):
                for j in range(N_ω): # loop over transitory income states
                    Xp = a_grid * (1 + rw) + Z[i] * ω[j] * g_t[t] # Implied cash on hand tomorrow 
                    index = range(i * N_a + i, (i + 1) * (N_a + 1) )
                    index_2 = range(i * N_ω + j,i * N_ω + j + 1)
                    Cp[:,index_2] = np.interp(Xp,Xw[index,t], Cw[index,t])
                    Vp[:,index_2] = np.interp(Xp,Xw[index,t], Vw[index,t])
            dVp = Cp ** (-γ)
            # Construct tomorrows value and tomorrows derivative wrt assets
            EV = β * np.dot(Vp , π.T)
            dV = β * np.dot(dVp , π.T) * (1 + rw) # RHS of Euler equation
            for i in range(N_z):
                index = range(i * (N_a+1) + 1, (i + 1) * (N_a+1))
                Cw[index  ,t-1:t] = dV[:,i:i+1] ** (-1/γ) # Use FOC to find consumption
                Xw[index ,t-1:t] =  Cw[index ,t-1:t] + a_grid # Implied cash on hand
                Vw[index ,t-1:t] = (Cw[index ,t-1:t] ** (1-γ) - 1 )/(1-γ) + EV[:,i:i+1] 
    return Vw, Cw, Xw



# Simulation

def simulate_model(T, rw, rr, Xw, Cw, Y_lower, t_w, t_r, g_t, ρ_z, N, N_a, Xr, Cr, Z, λ, ε_z, ε_ω, A, Z0,start,plot=True):
    A_sim = np.zeros((N,T))
    Z_sim = np.zeros((N,T))
    ω_sim = np.zeros((N,T))
    income_sim = np.zeros((N,T))
    Zi_sim = np.zeros((N,T),dtype=int)
    X_sim = np.zeros((N,T))
    C_sim = np.zeros((N,T))

    ## Initial Values
    A_sim[:,0] = A
    Z_sim[:,0] = Z0
    ω_sim[:,0] = ε_ω[:,0]
    Zi_sim[:,0] = dsearchn(Z, Z_sim[:,0])
    for t in range(T):
        t += 1
        for i in range(N):
            if t <= t_w:
                income_sim[i,t-1] = max(Z_sim[i,t-1] * ω_sim[i,t-1] * g_t[t-1],Y_lower)
                X_sim[i,t-1] = A_sim[i,t-1] * (1+rw) + income_sim[i,t-1] #cash in hand today
                index = range(Zi_sim[i,t-1] * (N_a+1), (Zi_sim[i,t-1] + 1) * (N_a+1)) # index to pick out set of asset choices | on (nearest) income state
                C_sim[i,t-1] = np.interp(X_sim[i,t-1],Xw[index,t-1], Cw[index,t-1])
                # Next period possible incomes
                ω_sim[i,t] = ε_ω[i,t-1]
                Z_sim[i,t] = np.exp(ρ_z * np.log(Z_sim[i,t-1]) + ε_z[i,t-1])
                Zi_sim[i,t] = dsearchn(Z, Z_sim[i,t].reshape(1,1))
            else:
                income_sim[i,t-1] = λ * g_t[t_w-1] + Z_sim[i,t_w-1]
                X_sim[i,t-1] = A_sim[i,t-1] * (1+rr) + income_sim[i,t-1]
                C_sim[i,t-1] = np.interp(X_sim[i,t-1],Xr[:,t-t_w-1], Cr[:,t-t_w-1])
                ...
            if t < T:
                A_sim[i,t] = X_sim[i,t-1] - C_sim[i,t-1]
        if t == t_w:
            print("Simulated model for working age")
        if t == t_w + t_r:
            print("Simulated model for retirement age")
    if plot:
        C_mean = np.sum(C_sim,axis=0)/N
        plt.plot(range(start,start + T),C_mean)
        A_mean = np.sum(A_sim,axis=0)/N
        plt.plot(range(start,start + T),A_mean)
        I_mean = np.sum(income_sim,axis=0)/N
        plt.plot(range(start,start + T),I_mean)
        plt.plot(range(start,start + T),np.zeros(T),linestyle='--')
        plt.legend(['Consumption','Assets','Income'])
        
        plt.grid()
        
    return A_sim, Z_sim, ω_sim, income_sim, Zi_sim, X_sim, C_sim

# Plotting
def plot_policy_over_ages(X,C,T,label, N_a, start, t_base,line45=True, yaxis='Consumption'):
    index = range(0, T, 5)
    for i in index:
        plt.plot(X[:N_a+1,i],C[:N_a+1,i],label = 'age  ' + str(start + t_base +i))
    if line45:plt.plot(X[:int(N_a/3)+1,i],X[:int(N_a/3)+1,i],label = '45 degree line',linestyle='--')
    plt.xlabel('X')
    plt.ylabel(yaxis)
    plt.grid()
    plt.legend()
    plt.title('Policy function at {}'.format(label))
    plt.show()
    plt.close()

def plot_policy_over_states(X,C,Z,age,t_w,start,N_z,N_a,line45=True,yaxis='Consumption'):
    if age > t_w + start:
        print("Age is greater than working age, please use plot_policy_over_ages")
    else:
        label = 'Working'
        t = age - start
        index = range(0, N_z)
        for i in index:
            plt.plot(X[i * (N_a+1):(i+1) * (N_a+1),t],C[i * (N_a+1):(i+1) * (N_a+1),t],label = 'z = ' + str(round(Z[i][0],2)))
        if line45:     plt.plot(X[:int(N_a/3)+1,t],X[:int(N_a/3)+1,t],label = '45 degree line',linestyle='--')
        plt.legend()
        plt.grid()
        plt.xlabel('X')
        plt.ylabel(yaxis)
        
        plt.title('Policy function at {}'.format(label))
        plt.show()
        plt.close()

# Tau Chen's code

def tauchenhussey(N,mu,rho,sigma):
	import numpy as np
	""" 
	Function tauchenhussey

	Purpose:    Finds a Markov chain whose sample paths
				approximate those of the AR(1) process
					z(t+1) = (1-rho)*mu + rho * z(t) + eps(t+1)
				where eps are normal with stddev sigma

	Format:     {Z, Zprob} = TauchenHussey(N,mu,rho,sigma,m)

	Input:      N         scalar, number of nodes for Z
		    mu        scalar, unconditional mean of process
		    rho       scalar
		    sigma     scalar, std. dev. of epsilons
		    baseSigma scalar, std. dev. used to calculate Gaussian
			quadrature weights and nodes, i.e. to build the
			grid. I recommend that you use
                        baseSigma = w*sigma +(1-w)*sigmaZ where sigmaZ = sigma/sqrt(1-rho^2),
				and w = 0.5 + rho/4. Tauchen & Hussey recommend
				baseSigma = sigma, and also mention baseSigma = sigmaZ.

	Output:     Z       N*1 vector, nodes for Z
				Zprob   N*N matrix, transition probabilities

	Author:		Benjamin Tengelsen, Brigham Young University (python)
				Martin Floden, Stockholm School of Economics (original)
				January 2007 (updated August 2007)

	This procedure is an implementation of Tauchen and Hussey's
	algorithm, Econometrica (1991, Vol. 59(2), pp. 371-396)
	"""
	w = 0.5 + rho/4
	sigmaZ = sigma/np.sqrt(1-rho^2)
	baseSigma = w*sigma +(1-w)*sigmaZ

	Z     = sp.zeros((N,1))
	Zprob = sp.zeros((N,N))
	[Z,w] = gaussnorm(N,mu,baseSigma**2)
	for i in range(N):
		for j in range(N):
			EZprime    = (1-rho)*mu + rho*Z[i]
			Zprob[i,j] = w[j] * st.norm.pdf(Z[j],EZprime,sigma) / st.norm.pdf(Z[j],mu,baseSigma)
		
	for i in range(N):
		Zprob[i,:] = Zprob[i,:] / sum(Zprob[i,:])
		
	return Z.T,Zprob


def gaussnorm(n,mu,s2):
	""" 
	Find Gaussian nodes and weights for the normal distribution
	n  = # nodes
	mu = mean
	s2 = variance
	"""
	[x0,w0] = gausshermite(n)
	x = x0*sp.sqrt(2.*s2) + mu
	w = w0/sp.sqrt(sp.pi)
	return [x,w]

	
def gausshermite(n):
	"""
	Gauss Hermite nodes and weights following 'Numerical Recipes for C' 
	"""

	MAXIT = 10
	EPS   = 3e-14
	PIM4  = 0.7511255444649425

	x = sp.zeros((n,1))
	w = sp.zeros((n,1))

	m = int((n+1)/2)
	for i in range(m):
		if i == 0:
			z = sp.sqrt((2.*n+1)-1.85575*(2.*n+1)**(-0.16667))
		elif i == 1:
			z = z - 1.14*(n**0.426)/z
		elif i == 2:
			z = 1.86*z - 0.86*x[0]
		elif i == 3:
			z = 1.91*z - 0.91*x[1]
		else:
			z = 2*z - x[i-1]
		
		for iter in range(MAXIT):
			p1 = PIM4
			p2 = 0.
			for j in range(n):
				p3 = p2
				p2 = p1
				p1 = z*sp.sqrt(2./(j+1))*p2 - sp.sqrt(float(j)/(j+1))*p3
			pp = sp.sqrt(2.*n)*p2
			z1 = z
			z = z1 - p1/pp
			if sp.absolute(z-z1) <= EPS:
				break
		
		if iter>MAXIT:
			return 'too many iterations'
		x[i,0]     = z
		x[n-i-1,0] = -z
		w[i,0]     = 2./pp/pp
		w[n-i-1,0] = w[i]
	
	x = x[::-1]
	return [x,w]


def stationary_markov(Pi, tol=1E-14):
    ########### Finding the stationary distribution of a Markov chain might be helpful later
    # start with uniform distribution over all states
    n = Pi.shape[0]
    pi = np.full(n, 1/n)
    
    # update distribution using Pi until successive iterations differ by less than tol
    for _ in range(10_000):
        pi_new = Pi.T @ pi
        if np.max(np.abs(pi_new - pi)) < tol:
            return pi_new
        pi = pi_new