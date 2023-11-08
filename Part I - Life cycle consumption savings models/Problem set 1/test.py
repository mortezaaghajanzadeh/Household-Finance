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
N_z = 15  ## number of states for the permanent income process
μ_z = 0.0 ## mean of the permanent income process
ρ_z = 1 ## persistence of the permanent income process
σ_η = 0.015   ## standard deviation of the permanent income process
N_ω = 5 ## number of states for the transitory income process
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
Z, ω, π, A, Z0, ε_z, ε_ω,a_grid = tk.initialize(T,N_z, μ_z, ρ_z, σ_η,N_ω,μ_ω,ρ_ω,σ_ω,N_a,a_max,ϕ,N,μ_A,σ_A,σ_z)

## Shock in the interest rate
N_r = 100 ## number of states for the interest rate process
μ = 0.04
r_f = 0.02
μ_r = np.log(r_f) + μ
σ_r = 0.18
r,π_r = tk.tauchenhussey(N_r, μ_r, 0, σ_r)
r = np.exp(r)[0].reshape(N_r,1) ## transitory income process
# print(r)


σ_r = 0.18
ε_r = np.exp(np.random.normal(-σ_r**2, σ_r, (N,T)))

N_α = 5
α_grid = np.linspace(0,1,N_α)

# %%

Vr = np.zeros((N_r * (N_a+1), t_r))
Cr = np.zeros((N_r * (N_a+1), t_r))
Xr = np.zeros((N_r * (N_a+1), t_r))
Ar = np.zeros((N_r * (N_a+1), t_r))
αr = np.zeros((N_r * (N_a+1), t_r))

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
    pension = λ * g_t[t_w-1]
    break
    Cp = np.zeros((N_a,N_r))
    Vp = np.zeros((N_a,N_r))
    for i in range(N_r):
        Xp = a_grid * (1 + r[i]) + pension ## cash-on-hand tomorrow
        index = range(i*(N_a+1)+1,(i+1)*(N_a+1))
        Cp[:,i:i+1] = np.interp(Xp,Xr[:,t+1], Cr[:,t+1]) # interpolate consumption
        Vp[:,i:i+1] = np.interp(Xp,Xr[:,t+1], Vr[:,t+1]) # interpolate consumption

    dVp = Cp ** (-γ) 
    EV = β * np.dot(dVp,π_r.T)
    for i in range(N_r):
        index = range(i*(N_a+1)+1,(i+1)*(N_a+1))
        dV = β * np.dot(dVp,π_r[i,:].T) * (1 + r[i])
        Cr[index,t:t+1] = dV.reshape(15,1) ** (-1/γ)
        Xr[index,t:t+1] = Cr[index,t:t+1] + a_grid
    

#%%
α = 0.5
from scipy.optimize import fsolve
from scipy.optimize import least_squares
def alpha(α,x):
    Xp =  x* (1 + r_f + α * (r - r_f)).T + pension
    C_tempt = np.interp(Xp,Xr[:,t+1], Cr[:,t+1])
    return( (C_tempt ** (-γ)  * (r - r_f).T ) @ π_r[i,:].T)[0]
alpha(α,a_grid[10])

least_squares(alpha,0.5,args = (a_grid[5],),bounds=(0 - epsilon,1 + epsilon))
x = a_grid[5]
# α_upper = 1
# α_lower = 0

# for _ in range(100):
#     if abs(alpha(α,x))<1e-4:
#         print(α , "is the solution")
#         break
#     if alpha(α,x)>0:
#         α_upper = α
#         α = (α_lower + α_upper)/2
#     else:
#         α_lower = α
#         α = (α_lower + α_upper)/2
#     print(α, alpha(α,x),α_upper,α_lower)

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