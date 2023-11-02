#%%
import Toolkit as tk
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d as interp
import matplotlib.pyplot as plt

set_seed = 13990509
np.random.seed(set_seed)
#%% Parameters
########## Model Parameters ##########
β = 0.945 ## Discount factor
γ = 2.0 ## Risk Aversion
start = 25 ## starting age
t_w = 40 ## working age
t_r = 35 ## retirement age
T = t_w + t_r ## total periods
λ = 0.6 ## replacement rate

rw = 0.03 ## Interest rate when working
rr = 0.03 ## Interest rate when retired

########## Income Process Parameters ##########
## Trend
g_t = pd.read_excel('Income_profile.xlsx')['Y'].to_numpy().flatten() 
# g_t = np.exp(g_t)
trend_pension = λ * g_t[t_w-1]

## Permanent Income Process
N_z = 100 ## number of states for the permanent income process
μ_z = 0.0 ## mean of the permanent income process
ρ_z = 1 ## persistence of the permanent income process
σ_η = 0.015  ## standard deviation of the permanent income process
Z,π_z = tk.tauchenhussey(N_z, μ_z, ρ_z, σ_η)
Z = np.exp(Z)[0].reshape(N_z,1) ## permanent income process

## Transitory Income Process
N_ω = 5 ## number of states for the transitory income process
μ_ω = 0.0 ## mean of the transitory income process
ρ_ω = 0 ## iid process
σ_ω = 0.1 ## standard deviation of the transitory income process
ω,π_ω = tk.tauchenhussey(N_ω, μ_ω, ρ_ω, σ_ω)
ω = np.exp(ω)[0].reshape(N_ω,1) ## transitory income process

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

########## Asset Grid ##########
a_max = 150 ## upper bound of the grid for assets
ϕ = 0  ## Borrowing Constraint
N_a = 1000 ## number of grid points for assets
a_grid = np.linspace(ϕ, a_max, N_a).reshape(N_a,1) ## grid for assets
# a_grid = tk.discretize_assets_double_exp(ϕ, a_max, N_a).reshape(N_a,1)


######### Thresholds ##########
vmin = -1.e10

import time
### Retirement
start_time = time.time()
Vr, Cr, Xr = tk.retirement(N_a,a_grid, rr, β, γ, t_r, t_w, g_t, λ,vmin)
print("--- %s seconds ---" % (time.time() - start_time))
### Working
start_time = time.time()
Vw, Cw, Xw = tk.working(N_z, N_ω, N_a, a_grid, Z, ω, π, rw, rr, Xr, Cr, Vr, β, γ, t_w, g_t,λ,vmin)
print("--- %s seconds ---" % (time.time() - start_time))
print("Found value functions and policy functions")
### Simulation
start_time = time.time()
A_sim, Z_sim, ω_sim, income_sim, Zi_sim, X_sim, C_sim =  tk.simulate_model(T, rw, rr, Xw, Cw, Y_lower, t_w, t_r, g_t, ρ_z, N, N_a, Xr, Cr, Z, λ, ε_z, ε_ω, A, Z0)
print("--- %s seconds ---" % (time.time() - start_time))
# %%



#%% Policy Functions
tk.plot_policy_over_ages(Xr,Cr,t_r,'Retirement', N_a, start, t_w,line45=True, yaxis='Consumption')
tk.plot_policy_over_ages(Xw,Cw,t_w,'Retirement', N_a, start, 0,line45=True, yaxis='Consumption')
#%%
tk.plot_policy_over_states(Xw,Cw,Z,25,t_w,start,N_z,N_a,line45=True,yaxis='Consumption')


Sw = Xw - Cw
Sr = Xr - Cr

tk.plot_policy_over_ages(Xw,Sw,t_w,'Retirement', N_a, start, t_w,line45=False, yaxis='Consumption')
tk.plot_policy_over_ages(Xr,Sr,t_r,'Retirement', N_a, start, 0,line45=False,yaxis='Savings')

#%% Value Functions

index = range(0, t_r, 10)
for i in index:
    plt.plot(Xr[:,i],Vr[:,i],label = 'age  ' + str(start + t_w +i))
legend = []
plt.legend()
plt.title('Value Function at Retirement')
plt.show()
plt.close()
index = range(0, t_w, 10)
for i in index:
    plt.scatter(Xw[:,i],Vw[:,i],label = 'age  ' + str(start + i))
# plt.scatter(Xw[:,index],Vw[:,index])
plt.legend()
plt.title('Value Function at Working')
plt.show()