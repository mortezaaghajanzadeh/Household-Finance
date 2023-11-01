import scipy.stats as st
import scipy as sp
import numpy as np
import numba
from scipy.spatial import KDTree

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


def backward_iteration(Va, Pi, a_grid, y, r, beta, eis):
    # step 1: discounting and expectations
    Wa = (beta * Pi) @ Va
    
    # step 2: solving for asset policy using the first-order condition
    c_endog = Wa**(-eis)
    coh = y[:, np.newaxis] + (1+r)*a_grid
    
    g = np.empty_like(coh)
    for e in range(len(y)):
        g[e, :] = np.interp(coh[e, :], c_endog[e, :] + a_grid, a_grid)
        
    # step 3: enforcing the borrowing constraint and backing out consumption
    g = np.maximum(g, a_grid[0])
    c = np.maximum(coh - g, 1E-9) # consumption must be strictly positive
    # step 4: using the envelope condition to recover the derivative of the value function
    Va = (1+r) * c**(-1/eis)
    
    return Va, g, c

def policy_ss(Pi, a_grid, y, r, beta, eis, tol=1E-9):
    # initial guess for Va: assume consumption 5% of cash-on-hand, then get Va from envelope condition
    coh = y[:, np.newaxis] + (1+r)*a_grid
    c = np.maximum(0.05*coh, 1E-9) # consumption must be strictly positive
    Va = (1+r) * c**(-1/eis)
    
    # iterate until maximum distance between two iterations falls below tol, fail-safe max of 10,000 iterations
    for it in range(10_000):
        Va, g, c = backward_iteration(Va, Pi, a_grid, y, r, beta, eis)
        
        # after iteration 0, can compare new policy function to old one
        if it > 0 and np.max(np.abs(g - g_old)) < tol:
            return Va, g, c
        
        g_old = g
		


def get_lottery(g, a_grid):
    # step 1: find the i such that a' lies between gridpoints a_k and a_(k+1)
    index_k = np.searchsorted(a_grid, g) - 1
    
    # step 2: obtain lottery probabilities pi
    alpha = (a_grid[index_k+1] - g)/(a_grid[index_k+1] - a_grid[index_k])
    
    return index_k, alpha

@numba.njit
def forward_iteration(D, index_k, alpha, Pi):
    Dend = np.zeros_like(D)
    for e in range(index_k.shape[0]):
        for a in range(index_k.shape[1]):
            # send pi(e,a) of the mass to gridpoint i(e,a)
            Dend[e, index_k[e,a]] += alpha[e,a]*D[e,a]
            
            # send 1-pi(e,a) of the mass to gridpoint i(e,a)+1
            Dend[e, index_k[e,a]+1] += (1-alpha[e,a])*D[e,a]
            
    return Pi.T @ Dend


def distribution_ss(Pi, g, a_grid, tol=1E-10):
    index_k, alpha = get_lottery(g, a_grid)
    
    # as initial D, use stationary distribution for e, plus uniform over a
    pi = stationary_markov(Pi)
    D = pi[:, np.newaxis] * np.ones_like(a_grid) / len(a_grid)
    
    # now iterate until convergence to acceptable threshold
    for _ in range(20_000):
        D_new = forward_iteration(D, index_k, alpha, Pi)
        if np.max(np.abs(D_new - D)) < tol:
            return D_new
        D = D_new


def steady_state(a_grid, y_grid, Pi, beta, eis, A_bar, r_min, r_max):
	for _ in range(100):
		r = (r_min + r_max)/2
		va, g, c = policy_ss(Pi, a_grid, y_grid, r, beta, eis, tol = 1e-9)
		dist = distribution_ss(Pi, g, a_grid, tol = 1e-9)
		demand = sum(sum(dist * g))
		if demand > A_bar:
			r_max = r
		else:
			r_min = r
		
		if np.abs(demand - A_bar) < 1e-4:
			break

	return r, g, c, dist




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