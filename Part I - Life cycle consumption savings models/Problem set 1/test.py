#%%
import Toolkit as tk
import numpy as np
import pandas as pd
#%% Parameters
########## Model Parameters ##########
β = 0.945 ## Discount factor
γ = 2.0 ## Risk Aversion
start = 25 ## starting age
t_w = 40 ## working age
t_r = 35 ## retirement age
T = t_w + t_r ## total periods
λ = 0.8

rw = 0.02 ## Interest rate when working
rr = 0.02 ## Interest rate when retired

########## Income Process Parameters ##########
## Trend
g_t = pd.read_excel('Income_profile.xlsx')['Y'].to_numpy().flatten()
# g_t = np.exp(g_t)

## Permanent Income Process
N_z = 3 ## number of states for the permanent income process
μ_z = 0.0 ## mean of the permanent income process
ρ_z = 1 ## persistence of the permanent income process
σ_η = 0.015 ## standard deviation of the permanent income process
Z,π_z = tk.tauchenhussey(N_z, μ_z, ρ_z, σ_η)
Z = np.exp(Z)[0].reshape(N_z,1) ## permanent income process

## Transitory Income Process
N_ω = 3 ## number of states for the transitory income process
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
ε_ω = np.random.normal(-σ_ω**2, σ_ω, (N,T))

Y_lower = 4.8

########## Asset Grid ##########
a_max = 150 ## upper bound of the grid for assets
ϕ = 0  ## Borrowing Constraint
N_a = 1000 ## number of grid points for assets
a_grid = np.linspace(ϕ, a_max, N_a).reshape(N_a,1) ## grid for assets
a_grid = tk.discretize_assets_double_exp(ϕ, a_max, N_a).reshape(N_a,1)


######### Thresholds ##########
vmin = -1.e10

#%% Retirement

Vr = np.zeros((N_a+1, t_r))
Cr = np.zeros((N_a+1, t_r))
Xr = np.zeros((N_a+1, t_r))

Xr[1:,-1:] = a_grid + 0.01
Cr[1:,-1:] = Xr[1:,-1:]
Vr[1:,-1:] = Cr[1:,-1:]**(1-γ)/(1-γ)
Vr[0,:] = vmin

trend_pension = λ * g_t[t_w-1]
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
#%% Working
Vw = np.zeros((N_z * (N_a+1), t_w))
Cw = np.zeros((N_z * (N_a+1), t_w))
Xw = np.zeros((N_z * (N_a+1), t_w))

for i in range(N_z):
    Vw[i * (N_a+1),:] = vmin

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
            for j in range(N_ω):
                index = range(i * N_a + i, (i + 1) * (N_a + 1) )
                index_2 = range(i * N_ω + j,i * N_ω + j + 1)
                Xp = a_grid * (1 + rw) + Z[i] * ω[j] * g_t[t] # Implied cash on hand tomorrow 
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
print("Found value function")
#%% Value Functions
import matplotlib.pyplot as plt

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
#%% Policy Functions
index = range(0, t_r, 10)
for i in index:
    plt.plot(Xr[:,i],Cr[:,i],label = 'age  ' + str(start + t_w +i))
plt.legend()
plt.title('Consumption Policy at Retirement')
plt.show()
plt.close()
index = range(0, t_w, 10)
for i in index:
    plt.plot(Xw[:int(N_a/2),i],Cw[:int(N_a/2),i],label = 'age  ' + str(start + i))
plt.legend()
plt.title('Consumption Policy at Working')
plt.show()
plt.close()



#%% Simulation
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
ω_sim[:,0] = W0
Zi_sim[:,0] = tk.dsearchn(Z, Z_sim[:,0])
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
            Zi_sim[i,t] = tk.dsearchn(Z, Z_sim[i,t].reshape(1,1))
        else:
            income_sim[i,t-1] = trend_pension + Z_sim[i,t_w-1]
            X_sim[i,t-1] = A_sim[i,t-1] * (1+rr) + income_sim[i,t-1]
            C_sim[i,t-1] = np.interp(X_sim[i,t-1],Xr[:,t-t_w-1], Cr[:,t-t_w-1])
            ...
        if t < T:
            A_sim[i,t] = X_sim[i,t-1] - C_sim[i,t-1]
    # if t == 1:
    #     break
    if t == t_w:
        print("Simulated model for working age")
    if t == t_w + t_r:
        print("Simulated model for retirement age")

# %%
C_mean = np.sum(C_sim,axis=0)/N
plt.plot(range(1,T+1),C_mean)
A_mean = np.sum(A_sim,axis=0)/N
plt.plot(range(1,T+1),A_mean)
#%%
A_sim[1,:]


