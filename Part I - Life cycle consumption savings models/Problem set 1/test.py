#%%
import Toolkit as tk
import numpy as np
import pandas as pd
#%%
########## Model Parameters ##########
β = 0.945 ## Discount factor
γ = 2.0 ## Risk Aversion
t_w = 40 ## working age
t_r = 35 ## retirement age
T = t_w + t_r ## total periods
g_t = pd.read_excel('Income_profile.xlsx')['Y'].to_numpy().flatten() ## growth rate of technology


rw = 0.02 ## Interest rate when working
rr = 0.02 ## Interest rate when retired

########## Income Process Parameters ##########
N_z = 3 ## number of states for the permanent income process
μ_z = 0.0 ## mean of the permanent income process
ρ_z = 1 ## persistence of the permanent income process
σ_z = 0.589 ## standard deviation of the permanent income process
Z,π_z = tk.tauchenhussey(N_z, μ_z, ρ_z, σ_z)


N = 3 ## number of states for the income process
rho = 1 ## persistence of the income process
mu = 0.0 ## mean of the income process
σ_η = 0.015 ## standard deviation of the income process


########## Simulation Parameters ##########

N = 1000
μ_A = 1.916 
σ_A = 2.129
A = np.exp(np.random.normal(μ_A, σ_A, N))


########## Asset Grid ##########
a_max = 150 ## upper bound of the grid for assets
ϕ = 0  ## Borrowing Constraint
N_a = 10 ## number of grid points for assets
a_grid = np.linspace(ϕ, a_max, N_a).reshape(N_a,1) ## grid for assets
######### Thresholds ##########
vmin = -1.e10

#%% Set up the problem

Vr = np.zeros((N_a+1, t_r))
Cr = np.zeros((N_a+1, t_r))
Xr = np.zeros((N_a+1, t_r))

Xr[1:,-1:] = a_grid + 0.01
Cr[1:,-1:] = Xr[1:,-1:]
Vr[1:,-1:] = Cr[1:,-1:]**(1-γ)/(1-γ)
Vr[0,:] = vmin
#%% Retirement

trend_pension = g_t[t_w-1]
for t in range(t_r-2, -1, -1):
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
