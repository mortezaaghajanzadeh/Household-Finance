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
n




ϕ = 0 ## Borrowing Constraint
a_bar = 150 ## upper bound of the grid for assets
#######################################
########## Income Process Parameters ##########
N = 3 ## number of states for the income process
rho = 1 ## persistence of the income process
mu = 0.0 ## mean of the income process
σ_η = 0.015 ## standard deviation of the income process
#######################################
