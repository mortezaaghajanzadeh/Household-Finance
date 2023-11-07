#%%
import Toolkit as tk
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d as interp
import matplotlib.pyplot as plt


#%% Parameters
set_seed = 13990509
np.random.seed(set_seed)
########## Model Parameters ##########
β = 0.945 ## Discount factor
γ = 2.0 ## Risk Aversion
start = 25 ## starting age
t_w = 40 ## working age
t_r = 35 ## retirement age
T = t_w + t_r ## total periods
λ = 0.6 ## replacement rate

rw = 0.02 ## Interest rate when working
rr = 0.02 ## Interest rate when retired

########## Income Process Parameters ##########
## Trend
g_t = pd.read_excel('Income_profile.xlsx')['Y'].to_numpy().flatten() 
# g_t = np.exp(g_t)
trend_pension = λ * g_t[t_w-1]

## Permanent Income Process
N_z = 10 ## number of states for the permanent income process
μ_z = 0.0 ## mean of the permanent income process
ρ_z = 1 ## persistence of the permanent income process
σ_η = 0.015  ## standard deviation of the permanent income process
Z,π_z = tk.tauchenhussey(N_z, μ_z, ρ_z, σ_η)
Z = np.exp(Z)[0].reshape(N_z,1) ## permanent income process

## Transitory Income Process
N_ω = 10 ## number of states for the transitory income process
μ_ω = 0.0 ## mean of the transitory income process
ρ_ω = 0 ## iid process
σ_ω = 0.1 ## standard deviation of the transitory income process
ω,π_ω = tk.tauchenhussey(N_ω, μ_ω, ρ_ω, σ_ω)
ω = np.exp(ω)[0].reshape(N_ω,1) ## transitory income process

## Shock in the interest rate
N_r = 3 ## number of states for the interest rate process
μ = 0.04
r_f = 0.02
μ_r = np.log(r_f) + μ
σ_r = 0.18
r,π_r = tk.tauchenhussey(N_r, μ_r, 0, σ_r)
r = np.exp(r)[0].reshape(N_r,1) ## transitory income process
print(r)




π = np.kron(π_z, π_ω) ## transition matrix for the income process


########## Simulation Parameters ##########

N = 1000
μ_A = 1.916 
σ_A = 2.129
A = np.random.normal(μ_A - σ_A**2 / 2, σ_A, N)
A = np.exp(A)
A[A < 0] = 0

μ_z = 0.0
σ_z = 0.015
Z0 = np.exp(np.random.normal(μ_z - σ_z**2 / 2 , σ_z, N))
ε_z = np.random.normal(-σ_z**2/2, σ_z, (N,T))

μ_ω = 0.0
W0 = np.exp(np.random.normal(μ_ω - σ_ω**2 / 2 , σ_ω, N))
ε_ω = np.exp(np.random.normal(μ_ω-σ_ω**2, σ_ω, (N,T)))

Y_lower = 4.8



σ_r = 0.18
ε_r = np.exp(np.random.normal(-σ_r**2, σ_r, (N,T)))

########## Asset Grid ##########
a_max = 150 ## upper bound of the grid for assets
ϕ = 0  ## Borrowing Constraint
N_a = 15 ## number of grid points for assets
a_grid = np.linspace(ϕ, a_max, N_a).reshape(N_a,1) ## grid for assets
# a_grid = tk.discretize_assets_double_exp(ϕ, a_max, N_a).reshape(N_a,1)


N_α = 5
α_grid = np.linspace(0,1,N_α)

######### Thresholds ##########
vmin = -1.e10

# import time
# ### Retirement
# start_time = time.time()
# Vr, Cr, Xr = tk.retirement(N_a,a_grid, rr, β, γ, t_r, t_w, g_t, λ,vmin)
# print("--- %s seconds ---" % (time.time() - start_time))
# ### Working
# start_time = time.time()
# Vw, Cw, Xw = tk.working(N_z, N_ω, N_a, a_grid, Z, ω, π, rw, rr, Xr, Cr, Vr, β, γ, t_w, g_t,λ,vmin)
# print("--- %s seconds ---" % (time.time() - start_time))
# print("Found value functions and policy functions")
# ### Simulation
# start_time = time.time()
# A_sim, Z_sim, ω_sim, income_sim, Zi_sim, X_sim, C_sim =  tk.simulate_model(T, rw, rr, Xw, Cw, Y_lower, t_w, t_r, g_t, ρ_z, N, N_a, Xr, Cr, Z, λ, ε_z, ε_ω, A, Z0,start)
# print("--- %s seconds ---" % (time.time() - start_time))
# %%

Vr = np.zeros((N_r * (N_a+1), t_r))
Cr = np.zeros((N_r * (N_a+1), t_r))
Xr = np.zeros((N_r * (N_a+1), t_r))
Ar = np.zeros((N_r * (N_a+1), t_r))

# Set the last period
for i in range(N_r):
    Vr[i*(N_a+1),:] = vmin
    index = range(i*(N_a+1)+1,(i+1)*(N_a+1))
    Xr[index,-1:] = a_grid 
    Cr[index,-1:] = Xr[index,-1:]
    Ar[index,-1:] = Xr[index,-1:] - Cr[index,-1:]
    Vr[index,-1:] = Cr[index,-1:]**(1-γ)/(1-γ)

# backward iteration
for t in range(t_r-1,0, -1):
    t -= 1
    # I have to think more about the pension income
    pension = trend_pension
    
    Cp = np.zeros((N_a,N_r))
    Xp = np.zeros((N_a,N_r))
    for i in range(N_r):
        Xp = a_grid * (1 + r[i]) + pension ## cash-on-hand tomorrow
        index = range(i*(N_a+1)+1,(i+1)*(N_a+1))
        print(index)
        
        # Cp = np.interp(Xp,Xr[:,t+1], Cr[:,t+1]) # interpolate consumption

    break

    Xp = a_grid * (1 + rr) + pension ## cash-on-hand tomorrow

    Cp = np.interp(Xp,Xr[:,t+1], Cr[:,t+1]) # interpolate consumption

    # Construct tomorrows value and tomorrows derivative wrt assets
    EV = β * np.interp(Xp,Xr[:,t+1], Vr[:,t+1])
    dV = β * Cp ** (-γ) * (1 + rr)

    Cr[1:,t:t+1] = dV ** (-1/γ) # Use FOC to find consumption policy
    Xr[1:,t:t+1] = Cr[1:,t:t+1] + a_grid # Implied cash on hand
    Vr[1:,t:t+1] = (Cp ** (1-γ) - 1 )/(1-γ) + EV
#%%
α = 0.9
Ap = a_grid
np.multiply((a_grid*(1+r_f + α * (r-r_f).T ) + pension - Ap) ** (-γ),np.repeat((r-r_f).T,N_a,0) ) @ π_r.T

             



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