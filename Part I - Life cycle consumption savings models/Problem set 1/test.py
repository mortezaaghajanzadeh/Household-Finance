#%%
import Toolkit as tk
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d as interp
import matplotlib.pyplot as plt
vmin = tk.vmin
epsilon = tk.epsilon
#%%
########## Model Parameters ##########
g_t = pd.read_excel('Income_profile.xlsx')['Y'].to_numpy().flatten() ## Income profile
β = 0.945 ## Discount factor
γ = 2.0 ## Risk Aversion
start = 25 ## starting age
t_w = 40 ## working age
t_r = 35 ## retirement age
T = t_w + t_r ## total periods
λ = 0.6  ## replacement rate
Y_lower = 4.8
rw = 0.02 ## Interest rate when working
rr = 0.02 ## Interest rate when retired

########## Income Process Parameters ##########
N_z = 3  ## number of states for the permanent income process
μ_z = 0.0 ## mean of the permanent income process
ρ_z = 1 ## persistence of the permanent income process
σ_η = 0.015   ## standard deviation of the permanent income process
N_ω = 3 ## number of states for the transitory income process
μ_ω = 0.0 ## mean of the transitory income process
ρ_ω = 0 ## iid process
σ_ω = 0.1  ## standard deviation of the transitory income process

########## Simulation Parameters ##########
N = 1000
μ_A = 1.916 
σ_A = 2.129
μ_z = 0.0
σ_z = 0.015
μ_ω = 0.0
########## Asset Grid ##########
a_max = 150 ## upper bound of the grid for assets
ϕ = 0  ## Borrowing Constraint
N_a = 15 ## number of grid points for assets


#%%
Z, ω, π_ω, π_z, π, A, Z0, ε_z, ε_ω,a_grid = tk.initialize(T,N_z, μ_z, ρ_z, σ_η,N_ω,μ_ω,ρ_ω,σ_ω,N_a,a_max,ϕ,N,μ_A,σ_A,σ_z)

## Shock in the interest rate
N_r = 5 ## number of states for the interest rate process
μ = 0.04 
r_f = 0.02
μ_r = np.log(1+r_f) + μ
σ_r = 0.18
r,π_r = tk.tauchenhussey(N_r, μ_r, 0, σ_r)
r = np.exp(r)[0].reshape(N_r,1) ## transitory income process
print(r)
r = r - 1


σ_r = 0.18
ε_r = np.exp(np.random.normal(-σ_r**2, σ_r, (N,T)))

N_α = 5
α_grid = np.linspace(0,1,N_α)
π = np.kron(np.kron(π_r, π_z), π_ω)


#%%
x = Xr[-1,5]
c = Cr[-1,5]
a = a_grid[-1]
α = 0.5
r_state = 0
z_state = 0
income = Z

def alpha_stocastic_income(α,x,c,a,r_state,z_state,π_r,π_z,income):
    A_t = x - c
    Y_t1 = income 
    X_t1 = Y_t1 +  A_t* (1 + r_f + α * (r - r_f)).T 
    c_t1 = X_t1 - a
    c_t1[c_t1<0] = epsilon
    M_t1 = β * c_t1 ** (-γ)  
    r_repeated = np.kron(r,np.ones((1,N_z))).T
    E =  ((M_t1 * (r_repeated - r_f)) @ π_r[r_state,:]).T @ π_z[z_state,:]
    return E

test_a = np.linspace(0,1,50)
res = []
for i in test_a:
    # res.append(alpha_stocastic_income(i,x,c,a,r_state,z_state,π_r,π_z,income))
    res.append(alpha(i,x,c,a,r_state,π_r,pension))
plt.plot(test_a,res)
fsolve(alpha,1,args = (x,c,a,r_state,π_r,pension))
# %%
α = 0.1
A_t = x - c
Y_t1 = pension 
X_t1 = Y_t1 +  A_t* (1 + r_f + α * (r - r_f)).T 
c_t1 = X_t1 - a
M_t1 = β * c_t1 ** (-γ)  
E =  (M_t1 * (r - r_f).T) @ π_r[r_state,:].T
E[0]
#%%
from scipy.optimize import fsolve
def alpha(α,x,c,a,r_index,π_r,income):
    A_t = x - c
    Y_t1 = income 
    X_t1 = Y_t1 +  A_t* (1 + r_f + α * (r - r_f)).T 
    c_t1 = X_t1 - a
    M_t1 = β * c_t1 ** (-γ)  
    E =  (M_t1 * (r - r_f).T) @ π_r[r_index,:].T
    return E[0]


def solve_alpha(x,c,a,r_index,π_r,income):
    alpha_0 , alpha_1 = alpha(0,x,c,a,r_index,π_r,income),alpha(1,x,c,a,r_index,π_r,income)
    if alpha_0 * alpha_1 > 0:
        if alpha_0 > 0: return 1
        else: return 0
    else:
        α = fsolve(alpha,0.5,args = (x,c,a,r_index,π_r,income))[0]
        if α >1: return 1
        else: return α

def solve_over_array(X,C,grid,r_index,π_r,income):
    αp = np.zeros((N_a,1))
    for n_a in range(0,N_a):
        αp[n_a] = solve_alpha(X[n_a],C[n_a],grid[n_a],r_index,π_r,income)
    return αp
    # αp = np.zeros((N_a,1))
    # for n_a in range(0,N_a):
    #     αp[n_a] = solve_alpha(Xp[n_a],Cp[n_a],a_grid[n_a],r_index = i,π_r = π_r,pension = pension)


def retirement_with_asset(r,π_r,γ,β,λ,g_t,t_w,t_r,N_a,a_grid,vmin):
    Vr = np.zeros((N_r * (N_a+1), t_r))
    Cr = np.zeros((N_r * (N_a+1), t_r))
    Xr = np.zeros((N_r * (N_a+1), t_r))
    αr = np.zeros((N_r * (N_a+1), t_r))

    # Set the last period
    for i in range(N_r):
        Vr[i*(N_a+1),:] = vmin
        index = range(i*(N_a+1)+1,(i+1)*(N_a+1))
        Xr[index,-1:] = a_grid 
        Cr[index,-1:] = Xr[index,-1:]
        Vr[index,-1:] = Cr[index,-1:]**(1-γ)/(1-γ)
        
    # backward iteration
    for t in range(t_r-1,0, -1):
        t -= 1
        # I have to think more about the pension income
        pension = λ * g_t[t_w-1]
        Cp = np.zeros((N_a,N_r))
        Vp = np.zeros((N_a,N_r))
        for i in range(N_r):
            Xp = a_grid * (1 + r[i]) + pension ## cash-on-hand tomorrow
            index = range(i*(N_a+1)+1,(i+1)*(N_a+1))
            Cp[:,i:i+1] = np.interp(Xp,Xr[index,t+1], Cr[index,t+1]) # interpolate consumption
            Vp[:,i:i+1] = np.interp(Xp,Xr[index,t+1], Vr[index,t+1]) # interpolate consumption
            for n_a in range(0,N_a):
                index = i*(n_a+1)+1
                αr[index,t:t+1] = solve_alpha(Xp[n_a],Cp[n_a,i:i+1],a_grid[n_a],r_index = i,π_r = π_r,income = pension)
        dVp = Cp ** (-γ) 
        EV = β * np.dot(dVp,π_r.T)
        for i in range(N_r):
            index = range(i*(N_a+1)+1,(i+1)*(N_a+1))
            dV = β * np.dot(dVp,π_r[i,:].T) * (1 + r[i])
            Cr[index,t:t+1] = dV.reshape(15,1) ** (-1/γ)
            Xr[index,t:t+1] = Cr[index,t:t+1] + a_grid
    return Cr,Xr,Vr,αr
Cr,Xr,Vr,αr = retirement_with_asset(r,π_r,γ,β,λ,g_t,t_w,t_r,N_a,a_grid,vmin)
    
αr[:,2]

#%%
# def working_with_asset(r,π_r,γ,β,λ,g_t,t_w,t_r,N_a,a_grid,vmin):
Vw = np.zeros(((N_r * N_z) * (N_a+1), t_w))
Cw = np.zeros(((N_r * N_z) * (N_a+1), t_w))
Xw = np.zeros(((N_r * N_z) * (N_a+1), t_w))
αw = np.zeros(((N_r * N_z) * (N_a+1), t_w))
for i in range(N_z):
    for j in range(N_r):
        Vw[(i+j) * (N_a+1),:] = vmin

pension = λ * g_t[t_w-1]
for t in range(t_w, 0, -1):
    if t == t_w: # Last period working
        Cp = np.zeros((N_a,1))
        Vp = np.zeros((N_a,1))
        counter = 0
        for i in range(N_r):
            Xp = a_grid * (1 + r[i]) + pension ## cash-on-hand tomorrow
            index = range(i*(N_a+1)+1,(i+1)*(N_a+1))
            Cp[:] = np.interp(Xp,Xr[index,0], Cr[index,0]) # interpolate consumption
            EV = β * np.interp(Xp,Xr[index,0], Vr[index,0])
            dV = β * Cp ** (-γ) * (1 + r[i])
            αp = solve_over_array(Xp,Cp,a_grid,r_index=i,π_r = π_r,income = pension)
            for j in range(N_z): 
                index = range(counter * (N_a+1) +1,(counter+1) * (N_a+1))
                counter += 1
                Cw[index,t-1:t] = dV ** (-1/γ)
                Xw[index,t-1:t] = Cw[index,t-1:t] + a_grid
                Vw[index,t-1:t] = (Cp ** (1-γ) - 1 )/(1-γ) + EV
                αw[index,t-1:t] = αp
                
    else:
        Cp = np.zeros((N_a, N_z * N_ω * N_r))
        Vp = np.zeros((N_a, N_z * N_ω * N_r))
        counter = 0
        for k in range(N_r):
            for i in range(N_z):
                for j in range(N_ω): # loop over transitory income states
                    Xp = a_grid * (1 + r[k]) + Z[i] * ω[j] * g_t[t] # Implied cash on hand tomorrow
                    index = range(counter * (N_a+1) +1,(counter+1) * (N_a+1))
                    index_2 = range(counter * N_ω +j,counter * N_ω +j+1)
                    Cp[:,index_2] = np.interp(Xp,Xw[index,t], Cw[index,t])
                    Vp[:,index_2] = np.interp(Xp,Xw[index,t], Vw[index,t])
                counter += 1
        dVp = Cp ** (-γ)
        # Construct tomorrows value and tomorrows derivative wrt assets
        EV = β * np.dot(Vp , π.T)
        counter = 0
        for k in range(N_r):
            dV = β * np.dot(dVp , π.T) * (1 + r[k])
            for i in range(N_z):
                index = range(counter * (N_a+1) +1,(counter+1) * (N_a+1))
                Cw[index  ,t-1:t] = dV[:,counter:counter+1] ** (-1/γ) # Use FOC to find consumption
                Xw[index ,t-1:t] =  Cw[index ,t-1:t] + a_grid # Implied cash on hand
                Vw[index ,t-1:t] = (Cw[index ,t-1:t] ** (1-γ) - 1 )/(1-γ) + EV[:,counter:counter+1]
                counter += 1
    ...

# αw[:,-1]
#%%
π.shape




#%%
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



#%%



#%%
A_t = x - c

Y_t1 = pension
X_t1 = Y_t1 + A_t * (1 + r_f + α * (r - r_f)).T 
X_t1
c_t1 = X_t1 - a
c_t1
c_t1 = np.interp(c_t1,Xr[:,t+1], Cr[:,t+1])
M_t1 = β * c_t1 ** (-γ)  
M_t1
E =  (M_t1 * (r - r_f).T) @ π_r[i,:].T
# E


# α_upper = 1
# α_lower = 0

# for _ in range(100):
#     res = alpha(α,x,c,a)
#     if abs(res)<1e-7:
#         print(α , "is the solution")
#         break
#     if res>0:
#         α_upper = α
#         α = (α_lower + α_upper)/2
#     else:
#         α_lower = α
#         α = (α_lower + α_upper)/2
#     print(α, res,α_upper,α_lower)

# %%

#%%
Xp =  x* (1 + r_f + α * (r - r_f)).T + pension
Xp

#%%
# Xr[1:,-1:] = a_grid + 0.01
# Cr[1:,-1:] = Xr[1:,-1:]
# Vr[1:,-1:] = Cr[1:,-1:]**(1-γ)/(1-γ)


# Vr[0,:] = vmin

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





#%%
# #%% Value Functions

# index = range(0, t_r, 10)
# for i in index:
#     plt.plot(Xr[:,i],Vr[:,i],label = 'age  ' + str(start + t_w +i))
# legend = []
# plt.legend()
# plt.title('Value Function at Retirement')
# plt.show()
# plt.close()
# index = range(0, t_w, 10)
# for i in index:
#     plt.scatter(Xw[:,i],Vw[:,i],label = 'age  ' + str(start + i))
# # plt.scatter(Xw[:,index],Vw[:,index])
# plt.legend()
# plt.title('Value Function at Working')
# plt.show()