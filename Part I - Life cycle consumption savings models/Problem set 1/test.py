#%%
import Toolkit as tk
import numpy as np
#%%
########## Model Parameters ##########
β = 0.945 ## Discount factor
γ = 2.0 ## Risk Aversion
ϕ = 0 ## Borrowing Constraint
a_bar = 150 ## upper bound of the grid for assets
#######################################
########## Income Process Parameters ##########
N = 3 ## number of states for the income process
rho = 1 ## persistence of the income process
mu = 0.0 ## mean of the income process
σ_η = 0.015 ## standard deviation of the income process
#######################################




########## Discritization Parameters ##########
n_a = 100 ## number of grid points for assets
n_z = 3 ## number of grid points for income
#######################################

########## Convergence Parameters ##########
epsilon_v = 10**(-6) ## tolerance for the value function iteration
epsilon_d = 10**(-8) ## tolerance for the distribution iteration
epsilon_m = 10**(-3) ## tolerance for the market clearing
max_iter_v = 1000 ## maximum number of iterations
max_iter_d = 100 ## maximum number of iterations
max_iter_m = 100 ## maximum number of iterations
#######################################
#%%
a_grid_0 = tk.discretize_assets_single_exp(ϕ, a_bar, n_a)
z_grid,π = tk.tauchenhussey(n_z,mu,rho,σ_η)
#%%
z_grid