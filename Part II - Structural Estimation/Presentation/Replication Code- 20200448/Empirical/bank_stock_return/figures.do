// file path
global path ".."
global folder "$path\data" // Specify data directory.
global working "$path\temp" // Specify working directory.
global resultfolder "$path\output" // Specify result directory.


cd "$working"

use "$folder\callreport_ffiec_merged.dta", clear
keep rssdid dateq assets liabilities equity deposits foreigndep intexp intexpdep domdepservicecharges /// 
loans intincloans loanleaselossprovision cash securities fedfundsrepoliab otherborrowedmoney ///
salaries exponpremises  netinc dividendoncommonstock loansleases* resloans_* transdep ///
nonintexp nonintinc intanddivincsecurities tradingassets fedfundsrepoasset

g tq=qofd(dateq)
format tq %tq

save callreports_small, replace


// aggegate data
use callreports_small, clear
g dep_total=deposits+foreigndep //total deposits
g borrowing=fedfundsrepoliab+otherborrowedmoney //borrowing
g net_borrowing=fedfundsrepoliab+otherborrowedmoney-securities-cash // net borrowing 
 
collapse (sum) assets liabilities equity dep_total foreigndep intexp intexpdep domdepservicecharges /// 
loans intincloans loanleaselossprovision cash securities fedfundsrepoliab otherborrowedmoney ///
salaries exponpremises  netinc dividendoncommonstock ///
loansleases* resloans_* transdep nonintexp nonintinc intanddivincsecurities tradingassets fedfundsrepoasset, by(date tq)
tsset tq
save callreports_agg, replace


 ***********************************
// rates
 ************************************
freduse DGS5 DGS3 FEDFUNDS MORTGAGE30US DGS30, clear
g tq=qofd(daten)
collapse DGS5 DGS3 FEDFUNDS MORTGAGE30US DGS30, by(tq)
save rate_temp, replace

use callreports_agg, clear
merge 1:1 tq using rate_temp
keep if _m==3
drop _m

g r_dep1=intexpdep/dep_total*4*100 // deposit rates
g r_fee1=domdepservicecharges/dep_total*4*100   // service fee for deposits
g r_loan2=(intincloans-loanleaselossprovision)/loans*4*100  // loan rates (default adjusted)
g r_loan1=(intincloans)/loans*4*100  // loan rates (not default adjusted)
g r_security=intanddivincsecurities/securities*4*100 // security rates
g mc_dep1=nonintexp/dep_total*4*100 // marginal cost
g non_int=nonintinc/dep_total*4*100  // non-interest income
g s_loan_dep=r_loan1-r_dep1
g net_mc=mc_dep1-non_int
g loss_prov=loanleaselossprovision/loans*4*100
g nonintexp_assets = nonintexp/assets*4*100
g salaries_assets = salaries/assets*4*100
g exponpremises_assets = exponpremises/assets*4*100

replace r_dep1=. if r_dep1==0
g s_dep=FEDFUNDS-r_dep1+r_fee1
g s_loan2=r_loan2-DGS5
g s_loan1=r_loan1-DGS5
sum loss_prov  if  tq>=tq(1985q1)
g s_loan3=r_loan1-DGS5 - `r(mean)'
g s_loan4=r_loan1-DGS3 - `r(mean)'

g loan_dep_ratio = loans/dep_total
g s_mort = MORTGAGE30US - DGS30
keep if  tq>=tq(1985q1)

sum s_dep s_loan1 loss_prov net_mc
sum loan_dep_ratio 

sort tq
tsset tq
outsheet  tq FEDFUNDS s_dep s_loan3 DGS3 DGS5 r_dep1 r_loan1  loss_prov nonintexp_assets using  rates_loan_deposits.csv, comma nolabel replace



lpoly s_loan3 FEDFUNDS if tq>=tq(1985q1) & s_loan3>-1, ci ///
  ylab(, nogrid)   xtitle("Fed Funds Rate") ytitle("Loan Spreads") title("") ///
graphregion(color(white)) graphregion(margin(zero))   legend( off ) ylabel(-1(1)4) yscale(range(-1 4))
capture: graph export "$resultfolder\Figures\lpoly_loan_spread.pdf", replace 



 lpoly s_dep FEDFUNDS  if tq>=tq(1985q1), ci ///
  ylab(, nogrid)   xtitle("Fed Funds Rate") ytitle("Deposit Spreads") title("") ///
graphregion(color(white)) graphregion(margin(zero))  legend( off ) ylabel( -1(1)4)
capture: graph export "$resultfolder\Figures\lpoly_deposit_spread.pdf", replace 

