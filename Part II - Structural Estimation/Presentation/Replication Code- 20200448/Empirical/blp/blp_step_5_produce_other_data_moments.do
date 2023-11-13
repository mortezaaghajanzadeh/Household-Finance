global path ".."
global folder "$path\data" // Specify data directory.
global working "$path\temp" // Specify working directory.
global resultfolder "$path\output" // Specify result directory.
global matlabfolder "$path\blp\blp_matlab_data" // Specify working directory for matlab.
global intermediate_output "$path\blp\intermediate_output" // Specify result directory.


cd "$working"

 ************************************
// macro moments
 ************************************
* annualized 
capture: erase __000000.txt
 freduse GDPC1 FEDFUNDS USLSTL, clear
 g year=year(daten)
 collapse GDPC1 FEDFUNDS USLSTL, by(year)
 g ln_GDPC1=log(GDP)
keep if  year>=1994 & year<=2017
g log_FEDFUNDS=log(FEDFUNDS)
g log_delta= log(USLSTL) // Loan Loss Reserve (balance sheet item, not income statement)

 sort year
regress ln_GDPC1 year
predict ln_GDPC1_detrended , resid
tsset year

sum log_FEDFUNDS, d
g log_FEDFUNDS_mean=`r(mean)'
g log_FEDFUNDS_sd=`r(sd)'

reg log_FEDFUNDS L.log_FEDFUNDS
predict log_FEDFUNDS_shock , resid
g log_FEDFUNDS_rho= _b[L.log_FEDFUNDS]

sum USLSTL, d

collapse (mean) log_FEDFUNDS_mean log_FEDFUNDS_sd log_FEDFUNDS_rho
xpose, clear  varname

outsheet using  "$intermediate_output\moment_ffr.csv", comma nolabel replace

foreach sample of numlist 1/3 {
#delimit ;
if `sample'==1 local year1=1994; if `sample'==1 local year2=2017;if `sample'==1 local bank1=0; if `sample'==1 local bank2=1; //whole sample
if `sample'==2 local year1=1994; if `sample'==2 local year2=2005;if `sample'==2 local bank1=0; if `sample'==2 local bank2=1; //pre 2006
if `sample'==3 local year1=2006; if `sample'==3 local year2=2017;if `sample'==3 local bank1=0; if `sample'==3 local bank2=1; //post 2006
if `sample'==4 local year1=1994; if `sample'==4 local year2=2017;if `sample'==4 local bank1=1; if `sample'==4 local bank2=1; //big bank 
if `sample'==5 local year1=1994; if `sample'==5 local year2=2017;if `sample'==5 local bank1=0; if `sample'==5 local bank2=0; //small bank 
#delimit cr
********************************************************************************************
* Payout policy: equal weighted
********************************************************************************************
use "$folder\bank_repurchase_div_yield.dta", clear
keep if fyear>=`year1' & fyear<=`year2' 
g total_payout_ratio = share_rep2_yield + div_yield // repo yield: non-cash dividends (fama and french 2001) => Change in treasury stocks (is available since 1982 on compustat).
replace total_payout_ratio=0 if total_payout_ratio<0
winsor total_payout_ratio, p(.01) g(total_payout_ratio_w)
sum total_payout_ratio_w, d
winsor MB, g(MB_w) p(.01)
collapse (mean) total_payout_ratio_w MB_w 
xpose, clear  varname
outsheet using  "$intermediate_output\moment_payout_ew_`sample'.csv", comma nolabel replace

********************************************************************************************
* Payout policy: value weighted
********************************************************************************************
use "$folder\bank_repurchase_div_yield.dta", clear
keep if fyear>=`year1' & fyear<=`year2' 
g total_payout_ratio = share_rep2_yield + div_yield // repo yield: non-cash dividends (fama and french 2001) => Change in treasury stocks (is available since 1982 on compustat).
replace total_payout_ratio=0 if total_payout_ratio<0
winsor total_payout_ratio, p(.01) g(total_payout_ratio_w)
sum total_payout_ratio_w, d
winsor MB, g(MB_w) p(.01)
replace BE0=0 if BE0<0
bysort fyear: g sum_mcap_c=sum(mcap_c)
g mcap_c_share = mcap_c/sum_mcap_c
collapse (mean) total_payout_ratio_w MB_w [w=mcap_c]
xpose, clear  varname
outsheet using  "$intermediate_output\moment_payout_vw_`sample'.csv", comma nolabel replace
}



 ***********************************
// rates
 ************************************
freduse DGS5 DGS3 FEDFUNDS MORTGAGE30US DGS30, clear
g tq=qofd(daten)
collapse DGS5 DGS3 FEDFUNDS MORTGAGE30US DGS30, by(tq)
save rate_temp, replace

********************************************************************************************
// Call Report: create variables
********************************************************************************************

use "$folder\data_extracted.dta", clear
g rssdid = rssd
g dateq=date
format dateq $tq
keep rssdid dateq assets liabilities equity deposits foreigndep savdep intexp intexpdep domdepservicecharges /// 
loans ciloans intincloans loanleaselossprovision cash securities fedfundsrepoliab otherborrowedmoney ///
salaries exponpremises  netinc dividendoncommonstock loansleases* resloan* securities* tradingassets timedep* subordinateddebt fedfundsrepoasset transdep ///
nonintexp nonintinc intanddivincsecurities netloanchargeoffs numemployees
g tq=qofd(dateq)
format tq %tq
g quarter=quarter(dateq)
save callreports_small, replace


use callreports_small, clear
merge m:1 tq using rate_temp
keep if _m==3
drop _m
keep if  tq>=tq(1994q1) & tq<=tq(2017q4)
g dep_total=deposits+foreigndep //total deposits
g r_dep=intexpdep/dep_total*400 // deposit rates
g r_fee=domdepservicecharges/deposits*400  // service fee for deposits
g r_loan=intincloans/loans*400  // loan rates
g borrowing2=assets-dep_total-equity
g liquidity=securities+cash+tradingassets+fedfundsrepoasset
g netloanchargeoffs_loans = netloanchargeoffs/loans*400
g salaries_at=salaries/assets[_n-1]*100
g exponpremises_at=exponpremises/assets[_n-1]*100
g leverage=assets/equity

* size 
gsort tq  - assets
bysort tq: g asset_rank = _n/_N
g big = (asset_rank<=.01)

foreach var of varlist  r_dep r_loan r_fee netloanchargeoffs_loans  leverage {
		winsor `var' , gen(temp) p(.01)
		replace `var' = temp 
		drop temp

   } 

foreach var of varlist   salaries_at exponpremises_at {
		winsor `var' , gen(temp) p(.1)
		replace `var' = temp 
		drop temp

   } 

g s_dep=FEDFUNDS-r_dep+r_fee
g s_loan=r_loan-DGS5

sort rssdid tq
tsset  rssdid tq
 
g loan_to_dep = loans/dep_total
g loan_to_assets = loans/assets
g liquidity_to_assets = liquidity/assets
g borrowing_assets=borrowing2/assets
g borrowing_dep_total=borrowing2/dep_total
g deposits_assets=dep_total/assets
g nonintexp_assets=nonintexp/assets*4 // nonintexp_assets ratio
g nonintinc_assets=nonintinc/assets*4  // non-interest income

foreach var of varlist loan_to_dep loan_to_assets liquidity_to_assets borrowing_assets borrowing_dep_total deposits_assets nonintexp_assets nonintinc_assets {
		winsor `var' , gen(temp) p(.01)
		replace `var' = temp 
		drop temp
}

* reserve requirement
foreach var of varlist transdep savdep {
g `var'_share= `var'/deposits
replace  `var'_share=. if  `var'_share==0
}
g effective_reserve_ratio=transdep_share*10+savdep_share*1
sum effective_reserve_ratio, d

* maturity
g total = resloans_less_3m+resloans_3m_1y+resloans_1y_3y+resloans_3y_5y+resloans_5y_15y+resloans_over_15y+ ///
loansleases_less_3m+loansleases_3m_1y+loansleases_1y_3y+loansleases_3y_5y+loansleases_5y_15y+loansleases_over_15y



* asset duration
gen qij_sumj=securitiestreasury_less_3m+securitiestreasury_3m_1y+securitiestreasury_1y_3y+securitiestreasury_3y_5y+securitiestreasury_5y_15y+securitiestreasury_over_15y ///
+securitiesrmbs_less_3m+securitiesrmbs_3m_1y+securitiesrmbs_1y_3y+securitiesrmbs_3y_5y+securitiesrmbs_5y_15y+securitiesrmbs_over_15y ///
+securitiesothermbs_less_3y+securitiesothermbs_over_3y ///
+resloans_less_3m+resloans_3m_1y+resloans_1y_3y+resloans_3y_5y+resloans_5y_15y+resloans_over_15y ///
+loansleases_less_3m+loansleases_3m_1y+loansleases_1y_3y+loansleases_3y_5y+loansleases_5y_15y+loansleases_over_15y ///
+cash+fedfundsrepoasset

gen tqij_sumj=0.125*securitiestreasury_less_3m+0.625*securitiestreasury_3m_1y+2*securitiestreasury_1y_3y+4*securitiestreasury_3y_5y+10*securitiestreasury_5y_15y+20*securitiestreasury_over_15y ///
+0.125*securitiesrmbs_less_3m+0.625*securitiesrmbs_3m_1y+2*securitiesrmbs_1y_3y+4*securitiesrmbs_3y_5y+10*securitiesrmbs_5y_15y+20*securitiesrmbs_over_15y ///
+1.5*securitiesothermbs_less_3y+5*securitiesothermbs_over_3y ///
+0.125*resloans_less_3m+0.625*resloans_3m_1y+2*resloans_1y_3y+4*resloans_3y_5y+10*resloans_5y_15y+20*resloans_over_15y ///
+0.125*loansleases_less_3m+0.625*loansleases_3m_1y+2*loansleases_1y_3y+4*loansleases_3y_5y+10*loansleases_5y_15y+20*loansleases_over_15y ///
+0*cash+0*fedfundsrepoasset

gen assets_dur=tqij_sumj/qij_sumj
drop qij_sumj-tqij_sumj

* loan duration
gen qij_sumj= resloans_less_3m+resloans_3m_1y+resloans_1y_3y+resloans_3y_5y+resloans_5y_15y+resloans_over_15y ///
+loansleases_less_3m+loansleases_3m_1y+loansleases_1y_3y+loansleases_3y_5y+loansleases_5y_15y+loansleases_over_15y ///
+cash+fedfundsrepoasset

gen tqij_sumj=0.125*resloans_less_3m+0.625*resloans_3m_1y+2*resloans_1y_3y+4*resloans_3y_5y+10*resloans_5y_15y+20*resloans_over_15y ///
+0.125*loansleases_less_3m+0.625*loansleases_3m_1y+2*loansleases_1y_3y+4*loansleases_3y_5y+10*loansleases_5y_15y+20*loansleases_over_15y ///
+0*cash+0*fedfundsrepoasset

gen loan_dur=tqij_sumj/qij_sumj
drop qij_sumj-tqij_sumj


* prepayment adjusted security duration
gen qij_sumj=securitiestreasury_less_3m+securitiestreasury_3m_1y+securitiestreasury_1y_3y+securitiestreasury_3y_5y+securitiestreasury_5y_15y+securitiestreasury_over_15y ///
+securitiesrmbs_less_3m+securitiesrmbs_3m_1y+securitiesrmbs_1y_3y+securitiesrmbs_3y_5y+securitiesrmbs_5y_15y+securitiesrmbs_over_15y ///
+securitiesothermbs_less_3y+securitiesothermbs_over_3y ///
+cash+fedfundsrepoasset

gen tqij_sumj=0.125*securitiestreasury_less_3m+0.625*securitiestreasury_3m_1y+2*securitiestreasury_1y_3y+4*securitiestreasury_3y_5y+10*securitiestreasury_5y_15y+20*securitiestreasury_over_15y ///
+0.125*securitiesrmbs_less_3m+0.625*securitiesrmbs_3m_1y+2*securitiesrmbs_1y_3y+4*securitiesrmbs_3y_5y+10*securitiesrmbs_5y_15y+20*securitiesrmbs_over_15y ///
+1.5*securitiesothermbs_less_3y+5*securitiesothermbs_over_3y ///
+ 0*cash+0*fedfundsrepoasset

gen secu_dur=tqij_sumj/qij_sumj
drop qij_sumj-tqij_sumj


* prepayment adjusted asset duration
gen qij_sumj=securitiestreasury_less_3m+securitiestreasury_3m_1y+securitiestreasury_1y_3y+securitiestreasury_3y_5y+securitiestreasury_5y_15y+securitiestreasury_over_15y ///
+securitiesrmbs_less_3m+securitiesrmbs_3m_1y+securitiesrmbs_1y_3y+securitiesrmbs_3y_5y+securitiesrmbs_5y_15y+securitiesrmbs_over_15y ///
+securitiesothermbs_less_3y+securitiesothermbs_over_3y ///
+resloans_less_3m+resloans_3m_1y+resloans_1y_3y+resloans_3y_5y+resloans_5y_15y+resloans_over_15y ///
+loansleases_less_3m+loansleases_3m_1y+loansleases_1y_3y+loansleases_3y_5y+loansleases_5y_15y+loansleases_over_15y ///
+cash+fedfundsrepoasset

gen tqij_sumj=0.125*securitiestreasury_less_3m+0.625*securitiestreasury_3m_1y+2*securitiestreasury_1y_3y+4*securitiestreasury_3y_5y+10*securitiestreasury_5y_15y+20*securitiestreasury_over_15y ///
+4*securitiesrmbs_less_3m+4*securitiesrmbs_3m_1y+4*securitiesrmbs_1y_3y+4*securitiesrmbs_3y_5y+4*securitiesrmbs_5y_15y+4*securitiesrmbs_over_15y ///
+4*securitiesothermbs_less_3y+4*securitiesothermbs_over_3y ///
+4*resloans_less_3m+4*resloans_3m_1y+4*resloans_1y_3y+4*resloans_3y_5y+4*resloans_5y_15y+4*resloans_over_15y ///
+0.125*loansleases_less_3m+0.625*loansleases_3m_1y+2*loansleases_1y_3y+4*loansleases_3y_5y+10*loansleases_5y_15y+20*loansleases_over_15y ///
+0*cash+0*fedfundsrepoasset

gen assets_dur_prepayment_adjusted=tqij_sumj/qij_sumj
drop qij_sumj-tqij_sumj

* prepayment adjusted loan duration
gen qij_sumj= resloans_less_3m+resloans_3m_1y+resloans_1y_3y+resloans_3y_5y+resloans_5y_15y+resloans_over_15y ///
+loansleases_less_3m+loansleases_3m_1y+loansleases_1y_3y+loansleases_3y_5y+loansleases_5y_15y+loansleases_over_15y ///
+cash+fedfundsrepoasset

gen tqij_sumj=4*resloans_less_3m+4*resloans_3m_1y+4*resloans_1y_3y+4*resloans_3y_5y+4*resloans_5y_15y+4*resloans_over_15y ///
+0.125*loansleases_less_3m+0.625*loansleases_3m_1y+2*loansleases_1y_3y+4*loansleases_3y_5y+10*loansleases_5y_15y+20*loansleases_over_15y ///
+0*cash+0*fedfundsrepoasset

gen loan_dur_prepayment_adjusted=tqij_sumj/qij_sumj
drop qij_sumj-tqij_sumj


* prepayment adjusted security duration
gen qij_sumj=securitiestreasury_less_3m+securitiestreasury_3m_1y+securitiestreasury_1y_3y+securitiestreasury_3y_5y+securitiestreasury_5y_15y+securitiestreasury_over_15y ///
+securitiesrmbs_less_3m+securitiesrmbs_3m_1y+securitiesrmbs_1y_3y+securitiesrmbs_3y_5y+securitiesrmbs_5y_15y+securitiesrmbs_over_15y ///
+securitiesothermbs_less_3y+securitiesothermbs_over_3y ///
+cash+fedfundsrepoasset

gen tqij_sumj=0.125*securitiestreasury_less_3m+0.625*securitiestreasury_3m_1y+2*securitiestreasury_1y_3y+4*securitiestreasury_3y_5y+10*securitiestreasury_5y_15y+20*securitiestreasury_over_15y ///
+4*securitiesrmbs_less_3m+4*securitiesrmbs_3m_1y+4*securitiesrmbs_1y_3y+4*securitiesrmbs_3y_5y+4*securitiesrmbs_5y_15y+4*securitiesrmbs_over_15y ///
+4*securitiesothermbs_less_3y+4*securitiesothermbs_over_3y ///
+ 0*cash+0*fedfundsrepoasset

gen secu_dur_prepayment_adjusted=tqij_sumj/qij_sumj
drop qij_sumj-tqij_sumj


//Liabilities duration
gen qij_sumj=timedeple100k_less_3m+timedeple100k_3m_1y+timedeple100k_1y_3y+timedeple100k_over_3y ///
+timedepge100k_less_3m+timedepge100k_3m_1y+timedepge100k_1y_3y+timedepge100k_over_3y ///
+transdep+savdep+fedfundsrepoliab+subordinateddebt

gen tqij_sumj=0.125*timedeple100k_less_3m+0.625*timedeple100k_3m_1y+2*timedeple100k_1y_3y+5*timedeple100k_over_3y ///
+0.125*timedepge100k_less_3m+0.625*timedepge100k_3m_1y+2*timedepge100k_1y_3y+5*timedepge100k_over_3y ///
+0*transdep+0*savdep+0*fedfundsrepoliab+5*subordinateddebt

gen liabilities_dur=tqij_sumj/qij_sumj
drop qij_sumj-tqij_sumj


* asset share 
bysort tq: egen sum_assets=sum(assets)
g assets_share=assets/sum_assets
save call_report_sample, replace

collapse (mean) borrowing_assets borrowing_dep_total loan_to_assets liquidity_to_assets deposits_assets nonintexp_assets nonintinc_assets leverage  ///
assets_dur assets_dur_prepayment_adjusted loan_dur_prepayment_adjusted secu_dur_prepayment_adjusted liabilities_dur [w=assets_share]
xpose, clear  varname

********************************************************************************************
// Charge off transition matrix
*******************************************************************************************
use call_report_sample, clear
keep if  tq>=tq(1994q1) & tq<=tq(2017q4)
sort rssdid tq

reg netloanchargeoffs_loans FEDFUNDS [iw=assets_share]
predict epsilon, res
g beta_FEDFUNDS = _b[FEDFUNDS]


xtile group=epsilon [w=assets_share], nq(10)
tab group, sum(epsilon) nofreq 

sort rssdid tq
bysort rssdid: g F_group=group[_n+4]
tab  group F_group [aw=assets_share], nofreq row

tab  group F_group [aw=assets_share], nofreq row matcell(freq) matrow(names)
matrix list freq
matrix list names
putexcel set "$intermediate_output\trans", sheet("transition") modify


bysort rssdid: g F_epsilon=epsilon[_n+4]
outsheet rssdid tq FEDFUNDS epsilon F_epsilon assets_share  using  "$intermediate_output\charge_off_raw.csv", comma nolabel replace

putexcel A1=matrix(names) B1=matrix(freq/r(N)*10) 
collapse (mean) epsilon beta_FEDFUNDS [aw=assets_share], by(group) 
xpose, clear  varname
outsheet using  "$intermediate_output\moment_charge_off_vw.csv", comma nolabel replace


********************************************************************************************
// Charge off transition matrix
*******************************************************************************************
foreach sample of numlist 1/5 {
#delimit ;
if `sample'==1 local year1=1994; if `sample'==1 local year2=2017;if `sample'==1 local bank1=0; if `sample'==1 local bank2=1; //whole sample
if `sample'==2 local year1=1994; if `sample'==2 local year2=2005;if `sample'==2 local bank1=0; if `sample'==2 local bank2=1; //pre 2006
if `sample'==3 local year1=2006; if `sample'==3 local year2=2017;if `sample'==3 local bank1=0; if `sample'==3 local bank2=1; //post 2006
if `sample'==4 local year1=1994; if `sample'==4 local year2=2017;if `sample'==4 local bank1=1; if `sample'==4 local bank2=1; //big bank 
if `sample'==5 local year1=1994; if `sample'==5 local year2=2017;if `sample'==5 local bank1=0; if `sample'==5 local bank2=0; //small bank 
#delimit cr

use call_report_sample, clear
keep if year(dateq)>=`year1' & year(dateq)<=`year2' 
keep if big == `bank1' | big == `bank2'
keep if quarter==4
g year=year(dofq(tq))
tsset rssdid year
g netloanchargeoffs_loans_L1y = L1.netloanchargeoffs_loans
outsheet rssdid year FEDFUNDS netloanchargeoffs_loans netloanchargeoffs_loans_L1y nonintexp_assets nonintinc_assets assets_share ///
using  "$intermediate_output\net_charge_off_raw_`sample'.csv", comma nolabel replace
}


********************************************************************************************
// Balance sheet moments: value weighted
********************************************************************************************
foreach sample of numlist 1/5 {
#delimit ;
if `sample'==1 local year1=1994; if `sample'==1 local year2=2017;if `sample'==1 local bank1=0; if `sample'==1 local bank2=1; //whole sample
if `sample'==2 local year1=1994; if `sample'==2 local year2=2005;if `sample'==2 local bank1=0; if `sample'==2 local bank2=1; //pre 2006
if `sample'==3 local year1=2006; if `sample'==3 local year2=2017;if `sample'==3 local bank1=0; if `sample'==3 local bank2=1; //post 2006
if `sample'==4 local year1=1994; if `sample'==4 local year2=2017;if `sample'==4 local bank1=1; if `sample'==4 local bank2=1; //big bank 
if `sample'==5 local year1=1994; if `sample'==5 local year2=2017;if `sample'==5 local bank1=0; if `sample'==5 local bank2=0; //small bank 
#delimit cr

use call_report_sample, clear
keep if year(dateq)>=`year1' & year(dateq)<=`year2' 
keep if big == `bank1' | big == `bank2'
* ciloan fraction
g ciloans_total_loans = ciloans/loans
replace ciloans_total_loans=1 if ciloans_total_loans>1
replace ciloans_total_loans=0 if ciloans_total_loans==.


* spread sensitivity
sort rssdid tq
tsset rssdid tq
reg D.s_dep D.FEDFUNDS [pw=assets_share]
display _b[D.FEDFUNDS]
g s_dep_FFR=_b[D.FEDFUNDS]

reg D.s_loan D.FEDFUNDS [pw=assets_share]
display _b[D.FEDFUNDS]
g s_loan_FFR=_b[D.FEDFUNDS]

bysort rssdid: egen mean_borrowing_assets=mean(borrowing_assets)
g borrowing_assets_demean = borrowing_assets-mean_borrowing_assets
replace borrowing_assets_demean=. if quarter!=4 // to compute annualized standard deviation
g borrowing_assets_sd = borrowing_assets
replace borrowing_assets_sd=. if quarter!=4 // to compute annualized standard deviation

bysort rssdid: egen mean_borrowing_dep_total=mean(borrowing_dep_total)
g borrowing_dep_total_demean = borrowing_dep_total-mean_borrowing_dep_total
replace borrowing_dep_total_demean=. if quarter!=4 // to compute annualized standard deviation
g borrowing_dep_total_sd = borrowing_dep_total
replace borrowing_dep_total_sd=. if quarter!=4 // to compute annualized standard deviation

* collapse
local asset_vars borrowing_assets borrowing_dep_total loan_to_assets liquidity_to_assets deposits_assets nonintexp_assets nonintinc_assets effective_reserve_ratio ciloans_total_loans leverage
local duration_vars loan_dur secu_dur assets_dur assets_dur_prepayment_adjusted loan_dur_prepayment_adjusted secu_dur_prepayment_adjusted liabilities_dur 
local spread_vars s_dep s_loan loan_to_dep s_dep_FFR s_loan_FFR netloanchargeoffs_loans

collapse (mean) `asset_vars' `duration_vars' `spread_vars'  ///
(sd) borrowing_assets_sd borrowing_assets_demean_sd=borrowing_assets_demean borrowing_dep_total_sd borrowing_dep_total_demean_sd=borrowing_dep_total_demean ///
assets_dur_sd =assets_dur_prepayment_adjusted [w=assets_share]
xpose, clear  varname
outsheet using  "$intermediate_output\moment_bank_vw_`sample'.csv", comma nolabel replace

}


********************************************************************************************
// Spread momoment
********************************************************************************************
foreach sample of numlist 1/5 {
#delimit ;
if `sample'==1 local year1=1994; if `sample'==1 local year2=2017;if `sample'==1 local bank1=0; if `sample'==1 local bank2=1; //whole sample
if `sample'==2 local year1=1994; if `sample'==2 local year2=2005;if `sample'==2 local bank1=0; if `sample'==2 local bank2=1; //pre 2006
if `sample'==3 local year1=2006; if `sample'==3 local year2=2017;if `sample'==3 local bank1=0; if `sample'==3 local bank2=1; //post 2006
if `sample'==4 local year1=1994; if `sample'==4 local year2=2017;if `sample'==4 local bank1=1; if `sample'==4 local bank2=1; //big bank 
if `sample'==5 local year1=1994; if `sample'==5 local year2=2017;if `sample'==5 local bank1=0; if `sample'==5 local bank2=0; //small bank 
#delimit cr

*use callreports_small, clear
use call_report_sample, clear
keep if year(dateq)>=`year1' & year(dateq)<=`year2' 
keep if big == `bank1' | big == `bank2'

 
collapse (sum) assets liabilities equity dep_total foreigndep intexp intexpdep domdepservicecharges /// 
loans intincloans loanleaselossprovision cash securities fedfundsrepoliab otherborrowedmoney ///
salaries exponpremises  netinc dividendoncommonstock ///
loansleases* resloans_* transdep nonintexp nonintinc intanddivincsecurities tradingassets fedfundsrepoasset, by(date tq)
tsset tq
save callreports_agg, replace


use callreports_agg, clear
merge m:1 tq using rate_temp
keep if _m==3
drop _m
keep if  tq>=tq(1985q1) & tq<=tq(2017q4) // this allows FFRs to have more variation

g r_dep1=intexpdep/dep_total*4*100 // deposit rates
g r_loan1=(intincloans)/loans*4*100  // loan rates (not default adjusted)
g loss_prov=loanleaselossprovision/loans*4*100

replace r_dep1=. if r_dep1==0
g s_dep=FEDFUNDS-r_dep1
sum loss_prov  
g s_loan4=r_loan1-DGS5- `r(mean)' // not adjusting prepayment risk when calculate rate

sort tq
reg s_dep FEDFUNDS if tq>=tq(1985q1)
reg s_dep FEDFUNDS if tq>=tq(1994q1)
reg D.s_dep D.FEDFUNDS if tq>=tq(1985q1)
reg D.s_dep D.FEDFUNDS if tq>=tq(1994q1)

reg s_loan4 FEDFUNDS if tq>=tq(1985q1)
reg s_loan4 FEDFUNDS if tq>=tq(1994q1)
reg D.s_loan4 D.FEDFUNDS if tq>=tq(1985q1)
reg D.s_loan4 D.FEDFUNDS if tq>=tq(1994q1)

reg s_dep FEDFUNDS if tq>=tq(1994q1)
g s_dep_FFR=_b[FEDFUNDS]

reg s_loan4 FEDFUNDS if tq>=tq(1994q1)
g s_loan_FFR=_b[FEDFUNDS]

collapse (mean) s_dep_FFR s_loan_FFR
xpose, clear  varname
outsheet using  "$intermediate_output\moment_spread_FFR_sensitivity_`sample'.csv", comma nolabel replace

}


********************************************************************************************
* HHI 
********************************************************************************************
use "$folder\FDIC_SOD.dta", clear
g entity = rssdhcr
replace entity = rssdid if rssdhcr==0
g one=1
collapse (mean) one, by( entity stcntybr year)
collapse (sum) numbanks=one, by(  stcntybr year)
sum numbanks, d


use "$folder\FDIC_SOD.dta", clear
g entity = rssdhcr
replace entity = rssdid if rssdhcr==0
collapse (sum) depsumbr, by( entity stcntybr year)
g asset = depsumbr
egen double sasset=sum(asset), by(stcntybr year)
bysort stcntybr year: g double shsqrd= (asset/sasset)^2 
egen double HHI=sum(shsqrd), by(stcntybr year)
save hhi_bank_level, replace



foreach sample of numlist 1/5 {
#delimit ;
if `sample'==1 local year1=1994; if `sample'==1 local year2=2017;if `sample'==1 local bank1=0; if `sample'==1 local bank2=1; //whole sample
if `sample'==2 local year1=1994; if `sample'==2 local year2=2005;if `sample'==2 local bank1=0; if `sample'==2 local bank2=1; //pre 2006
if `sample'==3 local year1=2006; if `sample'==3 local year2=2017;if `sample'==3 local bank1=0; if `sample'==3 local bank2=1; //post 2006
if `sample'==4 local year1=1994; if `sample'==4 local year2=2017;if `sample'==4 local bank1=1; if `sample'==4 local bank2=1; //big bank 
if `sample'==5 local year1=1994; if `sample'==5 local year2=2017;if `sample'==5 local bank1=0; if `sample'==5 local bank2=0; //small bank 
#delimit cr

use hhi_bank_level, clear
keep if year>=`year1' & year<=`year2' 
collapse HHI  [iw=sasset]
g J_hat = 1/HHI
xpose, clear  varname
outsheet using  "$intermediate_output\HHI_`sample'.csv", comma nolabel replace
}

* regional variation
use hhi_bank_level, clear
collapse HHI  [iw=sasset], by(stcntybr year)
keep if year==2017
sum HHI, d
g J_hat = 1/HHI
sum J_hat, d


********************************************************************************************
// Standard error
********************************************************************************************
use "$folder\FDIC_SOD.dta", clear
keep if year>=1994 & year<=2017
g entity = rssdhcr 
replace entity = rssdid if rssdhcr==.|rssdhcr==0
g entity_name = namehcr 
replace entity_name = namefull if entity_name==""
keep year entity rssdid entity_name
duplicates drop year entity rssdid, force
save link_entity_rssdid, replace


use "$folder\bank_repurchase_div_yield.dta", clear
keep if fyear>=1994 & fyear<=2017
g total_payout_ratio = share_rep2_yield + div_yield // repo yield: non-cash dividends (fama and french 2001) => Change in treasury stocks (is available since 1982 on compustat).
replace total_payout_ratio=0 if total_payout_ratio<0
winsor total_payout_ratio, p(.01) g(total_payout_ratio_w)
sum total_payout_ratio_w, d
winsor MB, g(MB_w) p(.01)
keep conml entity fyear total_payout_ratio_w MB_w
save temp_bank_div_yield, replace


// merge with bank-level data
foreach sample of numlist 1/5 {

*local sample = 1
#delimit ;
if `sample'==1 local year1=1994; if `sample'==1 local year2=2017;if `sample'==1 local bank1=0; if `sample'==1 local bank2=1; //whole sample
if `sample'==2 local year1=1994; if `sample'==2 local year2=2005;if `sample'==2 local bank1=0; if `sample'==2 local bank2=1; //pre 2006
if `sample'==3 local year1=2006; if `sample'==3 local year2=2017;if `sample'==3 local bank1=0; if `sample'==3 local bank2=1; //post 2006
if `sample'==4 local year1=1994; if `sample'==4 local year2=2017;if `sample'==4 local bank1=1; if `sample'==4 local bank2=1; //big bank 
if `sample'==5 local year1=1994; if `sample'==5 local year2=2017;if `sample'==5 local bank1=0; if `sample'==5 local bank2=0; //small bank 
#delimit cr

use call_report_sample, clear
g year = year(date)

merge m:1 rssdid year using link_entity_rssdid
drop if _m==2 // it seems that some non-bank financial institutions such as charles schwarb are in the sample
drop _m
g fyear=year
merge m:1 entity fyear using temp_bank_div_yield
drop if _m==2 // it seems that some non-bank financial institutions such as charles schwarb are in the sample
drop _m




keep if year(dateq)>=`year1' & year(dateq)<=`year2' 
keep if big == `bank1' | big == `bank2'
* ciloan fraction
g ciloans_total_loans = ciloans/loans
replace ciloans_total_loans=1 if ciloans_total_loans>1
replace ciloans_total_loans=0 if ciloans_total_loans==.



bysort rssdid: egen mean_borrowing_assets=mean(borrowing_assets)
g borrowing_assets_demean = borrowing_assets-mean_borrowing_assets
replace borrowing_assets_demean=. if quarter!=4 // to compute annualized standard deviation
g borrowing_assets_sd = borrowing_assets
replace borrowing_assets_sd=. if quarter!=4 // to compute annualized standard deviation

bysort rssdid: egen mean_borrowing_dep_total=mean(borrowing_dep_total)
g borrowing_dep_total_demean = borrowing_dep_total-mean_borrowing_dep_total
replace borrowing_dep_total_demean=. if quarter!=4 // to compute annualized standard deviation
g borrowing_dep_total_sd = borrowing_dep_total
replace borrowing_dep_total_sd=. if quarter!=4 // to compute annualized standard deviation

g net_nonintexp_assets=nonintexp_assets-nonintinc_assets
keep if quarter==4

outsheet year rssdid total_payout_ratio borrowing_assets borrowing_dep_total ///
s_dep s_loan deposits_assets net_nonintexp_assets leverage MB_w  assets_share using  "$intermediate_output\standard_error_`sample'.csv", comma nolabel replace

collapse (mean) year rssdid total_payout_ratio_w borrowing_assets borrowing_dep_total ///
s_dep s_loan deposits_assets net_nonintexp_assets leverage MB_w  [iw=assets_share]
}


 !rmdir "$working"  /s /q // delete files in temp folder
