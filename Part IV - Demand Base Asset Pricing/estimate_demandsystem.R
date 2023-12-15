# estimate_demandsystem.R
#
# (cc) Huebner
#
# Prepare estimation data for demand estimation and estimate
#
#
# Created on     March     21 2015
# Last modified  Dec       09 2023
#
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# ------ load packages
library(fst);
library(dplyr); 
library(lubridate); 
library(stringr)
library(skimr)
library(data.table);
library(statar);
library(Hmisc);
library(fixest);
library(ggplot2);
library(tidyr);
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
ownership_path <- './input/ownership_raw.fst' # 13F input file (lightly cleaned)
crsp_path <- './input/crsp_compu_chars.fst' # stock data (like in KY)
stocks_min <- 50 # minimum number of stocks for an institution
managers_min <- 10 # minimum number of owners per stock
write_file <- './output/estimation_sample.csv' # where to save sample
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# --- LOAD OWNERSHIP DATA

dt_in <- read_fst(ownership_path, as.data.table = T)
# num_inv_tot is number of investors a stock has at a point in time
# num_stocks is number of stocks any investor holds at a point in time
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# --- LOAD COMPUSTAT CHARACTERISTICS

# in KY: LNme (in million $) + LNbe + investment (Gat = growth of at) +
#   profitability + divA/be + beta (NOT INCLUDED HERE => equilibrium object!)
# KY winsorize investment and profitability and beta at 97.5 and 2.5,
#   and div at 97.5; this has already been done here
dt_crsp <- read_fst(crsp_path, as.data.table = T)
dt_crsp[, LN_be := log(be)]
dt_crsp[, c('be', 'me') := NULL] # me column already included in dt_in

# MERGE DATA
dt_in[, LN_me := log(me/1000^2)] # to get to million $
dt_in <- merge(dt_in, dt_crsp, all.x = T, by = c('permno', 'datem')) # left join
rm(dt_crsp); # free up memory (important when using more data than here)

# Data quality checks
dt_in <- dt_in[shares <= shrout*1000]
# drops 58 obs where investors has more shares than shrout
# Possibly due to short interest? => paper by Mainardi is more careful about this
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# --- OUTSIDE ASSET DEFINITION I
# part 1: NAs in characteristics or shrcd not 10 or 11
setorder(dt_in, datem, permno, mgrno)
dt_in[, inside := 0 + (!is.na(divA_be) & !is.na(profit) & !is.na(Gat) &
                       !is.na(LN_be) & !is.na(LN_me) & !is.na(ind_ff12) &
                       shrcd %in% c(10, 11))]

# part 2: not held by enough investors
dt_in[, num_inv_tot2 := .N, by = .(permno, datem)]
dt_in[num_inv_tot2 < managers_min, inside := 0]
dt_in[inside == 1]
# 11941039 obs with inside == 1
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# --- APPLY FILTERS I
# part 1: investor does not hold at least one outside asset
dt_in[, N_outside := .N - sum(inside), by = .(datem, mgrno)] # 39854
dt_in <- dt_in[N_outside > 0] # to 14.83mn obs left

# part 2: investor holds at least 50 inside assets
dt_in[, num_stocks := sum(inside == 1), by = .(mgrno, datem)] 
dt_in <- dt_in[num_stocks >= stocks_min] # from 14.83mn to 14.00mn
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# --- OUTSIDE ASSET DEFINITION II
# repeat the num_inv_tot < 10
dt_in[, num_inv_tot2 := .N, by = .(permno, datem)]
dt_in[num_inv_tot2 < managers_min, inside := 0]
dt_in[, c('num_inv_tot2', 'num_inv_tot', 'N_outside', 'num_stocks') := NULL]
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# --- DEFINE HH HOLDINGS
# also includes small institutions and institutions with few stocks filtered out
# basically residual of everything that is not in the data set otherwise
dt_hh <- unique(dt_in[, .(shares = (shrout*1000 - sum(shares)),
                          mgrno = 0), by = c('datem', 'permno') ])
dt_hh[, shares := ifelse(shares > 0, shares, 0)] # cannot have negative numbers 
tmp <- unique(dt_in[, -c('mgrno', 'shares')]) 
dt_hh <- merge(dt_hh, tmp, by = c('datem', 'permno'))
setcolorder(dt_hh, colnames(dt_in))
dt_in <- rbind(dt_in, dt_hh)
dt_in <- dt_in[shares != 0] # only mgrno 0 # to 14.07mn observations
rm(tmp, dt_hh)
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# --- DEFINE AUM
dt_in[, holdings := shares*abs(prc)]  

# outside assets
dt_outside <- dt_in[inside == 0, .(AUM_outside = sum(holdings)), by = c('mgrno', 'datem')]
# inside AUM
dt_inside <- dt_in[inside == 1, .(AUM_inside = sum(holdings)), by = c('mgrno', 'datem')]
# merge
dt_in <- merge(dt_in, dt_outside, by = c('mgrno', 'datem'), all.x = T)
dt_in <- merge(dt_in, dt_inside, by = c('mgrno', 'datem'), all.x = T)
dt_in[, AUM := AUM_inside + AUM_outside]
rm(dt_outside, dt_inside)
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# --- DEFINE LOG ODDS RATIO
dt_in[, log_odds := log(holdings/AUM) - log(AUM_outside/AUM)]
# 14068800 observations
# 11324502 observations after inside == 1
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# --- DEFINE K-Y INVESTMENT UNIVERSE
# expand the data to include permno in the opportunity set of an investor if
#   it is in a 3 year interval 
setorder(dt_in, mgrno, permno, dateq)
dt_in[, dateq := as.integer(dateq) ]
dt_in[, dateq_lead := shift(dateq, n=1L, type="lead", fill=NA), by = .(mgrno, permno) ]
dt_in[, dateq_n_rep := pmin(dateq_lead - dateq - 1, 11)]
# quarter difference between current hold date and next hold date (capped at 11qtr)
dt_in[is.na(dateq_n_rep), dateq_n_rep := 11]
# if missing, then current date is last stock holding and default to 11qtr
dt_in[, num:=1:.N]
# creates dateq_n_rep duplicates: 
dt_in <- dt_in[rep(num, (dateq_n_rep + 1))] 

dt_in[, ist_count := 0:(.N-1), by = c('mgrno', 'permno', 'datem')]
# add a counter on ist level (enumerate the duplicated)
dt_in[, datem := datem + ist_count*3]
# get the corresponding date for the duplicated holding datapoints
dt_in[, imputed := 0]
dt_in[ist_count != 0, imputed := 1]
# flag imputed columns for future dates: this means no actual holding but in universe
# cleaning up
dt_in <- dt_in[datem <= max(dt_in[imputed == 0][["datem"]]) ]
# filter imputed columns after sample end
dt_in[, c('dateq_n_rep', 'num', 'ist_count', 'dateq', 'dateq_lead') := NULL]
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# --- DEFINE K-Y INSTRUMENT FOR ME
# need to set AUM to zero for imputed
setnames(dt_in, "AUM", "AUM_old") # careful as the imputed AUM is from forward dates
dt_AUM <- unique(dt_in[imputed == 0, .(mgrno, datem, AUM=AUM_old)])  
dt_in <- merge(dt_in, dt_AUM, all.x = T, by = c('datem', 'mgrno'))
dt_in <- dt_in[!is.na(AUM)] # where the manager is not in 13F any more

# PRICE INSTRUMENT CALCULATION: see equation 19 nber WP KY
# Step 1: total AUM of manager divided by sum of all INSIDE assets (see above)
#   in investment universe of manager + outside asset
dt_in[, instrument_st := AUM/(1 + sum(inside)) , by = c('mgrno', 'datem')] 
dt_in[, num_stocks_IU := sum(inside), by = c('mgrno', 'datem')]

# now instead of summing over j \neq i, sum over all j and subtract i, then log
dt_in <- dt_in[inside == 1]
dt_in[, LN_me_IV  := log((sum(instrument_st, na.rm = T) - instrument_st)/1E6),
          by = c('datem', 'permno')]
dt_in <- dt_in[imputed == 0]
# 11324502 observations
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# --- final steps 
dt_in[, num_stocks := .N, by = .(mgrno, datem)]
dt_in[, num_inv_tot := .N, by = .(permno, datem)]

dt_in[, date := as.monthly(datem) ]
dt_in[, dateym := year(date)*100+month(date) ]
dt_in[, c("datem", "date", "instrument_st", "imputed", "inside") := NULL ]
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# --- filter to 2020 to speed up remainder of code

dt_in <- dt_in[dateym >= 202001]
# 2,915,606 observations
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# --- Estimate demand system by mgrno & date

# Define the IV regression formula
formula_main <- log_odds ~ LN_be + profit + Gat + divA_be | LN_me ~ LN_me_IV
formula_OLS  <- log_odds ~ LN_me + LN_be + profit + Gat + divA_be

# Run the IV regression by group using the fixest package
dt_estimates <- dt_in[, {
  fit <- fixest::feols(formula_main, data = .SD)
  fitOLS <- fixest::feols(formula_OLS, data = .SD)
  .(ela = 1 - coef(fit)['fit_LN_me'],
    r2 = fixest::r2(fit)['ar2'],
    elaOLS = 1 - coef(fitOLS)['LN_me'],
    r2OLS = fixest::r2(fitOLS)['ar2'],
    AUM = AUM[1])
}, by = .(mgrno, dateym)]

# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# --- Aggregate results

# Aggregate by date and plot elasticities
dt_plot <- 
dt_estimates[, .(Average = wtd.mean(pmax(ela, 0), weights = AUM),
                 Median = median(pmax(ela, 0)),
                 AverageOLS = wtd.mean(pmax(elaOLS, 0), weights = AUM),
                 MedianOLS = median(pmax(elaOLS, 0)),
                 AverageNoHH = wtd.mean(pmax(ela, 0), weights = AUM*(mgrno != 0)),
                 AverageNoHHOLS = wtd.mean(pmax(elaOLS, 0), weights = AUM*(mgrno != 0))), by = .(dateym)]

dt_plot <- 
  dt_plot %>%
  pivot_longer(cols = c('Average', 'Median', 'AverageOLS', 'MedianOLS',
                        'AverageNoHH', 'AverageNoHHOLS'),
               names_to='Measure',
               values_to='Elasticity')

ggplot(dt_plot, aes(x = dateym, y = Elasticity)) +
  geom_point(aes(col = Measure, shape = Measure)) +
  ylim(c(0, 0.75)) + 
  theme_bw()


# Aggregate by date and plot r2
dt_plot <- 
  dt_estimates[, .(Average = wtd.mean(r2, weights = AUM),
                   Median = median(r2),
                   AverageOLS = wtd.mean(r2OLS, weights = AUM),
                   MedianOLS = median(r2OLS),
                   AverageNoHH = wtd.mean(r2, weights = AUM*(mgrno != 0)),
                   AverageNoHHOLS = wtd.mean(r2OLS, weights = AUM*(mgrno != 0))), by = .(dateym)]

dt_plot <- 
  dt_plot %>%
  pivot_longer(cols = c('Average', 'Median', 'AverageOLS', 'MedianOLS',
                        'AverageNoHH', 'AverageNoHHOLS'),
               names_to='Measure',
               values_to='R2')

ggplot(dt_plot, aes(x = dateym, y = R2)) +
  geom_point(aes(col = Measure, shape = Measure)) +
  ylim(c(0, 0.75)) + 
  theme_bw()
