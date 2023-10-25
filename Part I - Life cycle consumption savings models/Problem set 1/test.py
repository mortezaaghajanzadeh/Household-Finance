#%%
import Toolkit as tk
#%%
########## Model Parameters ##########
β = 0.945 ## Discount factor
γ = 2.0 ## Risk Aversion
ϕ = 0 ## Borrowing Constraint
a_bar = 150 ## upper bound of the grid for assets
#######################################
########## Income Process Parameters ##########
N = 3 ## number of states for the income process
rho = 0.95 ## persistence of the income process
mu = 0.0 ## mean of the income process
sigma_eps = 0.015 ## standard deviation of the income process
#######################################

########## Discritization Parameters ##########
M = 1000 ## number of grid points for assets
N = 3 ## number of grid points for income
#######################################

########## Convergence Parameters ##########
epsilon_v = 10**(-6) ## tolerance for the value function iteration
epsilon_d = 10**(-8) ## tolerance for the distribution iteration
epsilon_m = 10**(-3) ## tolerance for the market clearing
max_iter_v = 1000 ## maximum number of iterations
max_iter_d = 100 ## maximum number of iterations
max_iter_m = 100 ## maximum number of iterations
#######################################
