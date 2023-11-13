// file path
global path ".."
global folder "$path\data" // Specify data directory.
global working "$path\temp" // Specify working directory.
global resultfolder "$path\output" // Specify result directory.
global matlabfolder "$path\blp\blp_matlab_data" // Specify working directory for matlab.


cd "$working"


clear
cap rm "__000000.txt"
freduse FEDFUNDS UNRATE CPIAUCSL, clear 
gen tq=qofd(daten)
format %tq tq
g inf = (CPIAUCSL/CPIAUCSL[_n-1]-1)*100
collapse FEDFUNDS inf UNRATE, by(tq)
local serieslist FEDFUNDS
foreach var of local serieslist {
g `var'_chg=`var'-`var'[_n-1]
}
save "FEDFUNDS",replace


//******************************************************************************
//Kuttner shocks
//******************************************************************************
import excel using "$folder\daily FF surprises July 2008.xls", clear first 
rename A date
g tq = qofd(date)
format tq %tq
foreach var of varlist Change Surprise Expected {
replace `var'=`var'/100 
}

collapse (sum) Change Surprise Expected, by(tq)
tsset tq
tsfill // not every month has FOMC meetings

foreach var of varlist Change Surprise Expected {
replace `var'=0 if `var'==.
}
merge 1:1 tq using "FEDFUNDS"
keep if _m==3
drop _m


keep tq  Surprise Expected inf UNRATE

save "kuttner_shocks.dta",replace



insheet using "$folder\Gertler_Karadi_factor_data.csv", clear
g tm = ym(year,month)
g tq = qofd(dofm(tm))
collapse (sum) mp1_tc-ed4_tc, by(tq)
save "Gertler_Karadi_factor_data_tq", replace





 ***********************************
// rates
 ************************************
freduse DGS5 DGS3 FEDFUNDS MORTGAGE30US DGS30, clear
g tq=qofd(daten)
collapse DGS5 DGS3 FEDFUNDS MORTGAGE30US DGS30, by(tq)
save "rate_temp", replace

********************************************************************************************
// Call Report: create variables
********************************************************************************************
use "$folder\data_extracted.dta", clear
g rssdid = rssd
g dateq=date
format dateq $tq
keep rssdid dateq assets liabilities equity deposits foreigndep savdep intexp intexpdep domdepservicecharges /// 
loans ciloans intincloans intincnet loanleaselossprovision cash securities fedfundsrepoliab otherborrowedmoney ///
salaries exponpremises  netinc dividendoncommonstock loansleases* resloan* securities* tradingassets timedep* subordinateddebt fedfundsrepoasset transdep ///
nonintexp nonintinc intanddivincsecurities netloanchargeoffs numemployees
g tq=qofd(dateq)
format tq %tq
g quarter=quarter(dateq)
save callreports_small, replace


use callreports_small, clear
merge m:1 tq using "rate_temp"
keep if _m==3
drop _m
keep if  tq>=tq(1994q1) & tq<=tq(2017q4)
g dep_total=deposits+foreigndep //total deposits
g r_dep=intexpdep/dep_total*400 // deposit rates
g r_fee=domdepservicecharges/deposits*400  // service fee for deposits
g r_loan=intincloans/loans*400  // loan rates
g borrowing=fedfundsrepoliab+otherborrowedmoney //borrowing
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

foreach var of varlist  r_dep r_loan r_fee netloanchargeoffs_loans  leverage    {
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
g intincnet_assets=intincnet/assets*4

foreach var of varlist loan_to_dep loan_to_assets liquidity_to_assets borrowing_assets ///
borrowing_dep_total deposits_assets nonintexp_assets nonintinc_assets intincnet_assets {
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
save "call_report_sample", replace


use "call_report_sample", clear
keep if  tq>=tq(1994q1) & tq<=tq(2017q4)
merge m:1 tq using "kuttner_shocks.dta"
keep if _merge==3
drop _m

tsset rssdid tq


foreach var of varlist assets loans securities deposits borrowing2 equity  {
g ln_`var'=ln(`var')
}


foreach var of varlist  r_dep r_loan FEDFUNDS s_dep s_loan {
replace `var'=`var'/100
}
save call_report_local_projection, replace






use call_report_local_projection, clear
g net_nonintexp_assets = nonintexp_assets - nonintinc_assets 
g s_loan_ma = (L4.s_loan+s_loan+F4.s_loan)/3
label var FEDFUNDS "Fed funds rate"
tsset rssdid tq
g cons=1
set more off

foreach num of numlist 2 {
local contr ln_assets ln_loans ln_securities ln_deposits ln_borrowing2  ln_equity r_dep r_loan intincnet_assets nonintexp_assets nonintinc_assets 
local macro  
local horizon=`num'*4
local list ln_assets ln_loans ln_securities ln_deposits ln_borrowing2   
eststo clear
foreach var of local list {
eststo `var': reghdfe F`horizon'.`var'  `contr' `macro'  (D.FEDFUNDS=Surprise), ab(rssdid) cluster(rssdid tq)
estadd local FE "Yes"
}

esttab using "$resultfolder\Tables/LP_quantities_`num'y.tex", ///
replace b(3) se(3) star(* 0.10 ** 0.05 *** 0.01  ) ///
	mtitle("Assets"   "Loans" "Securities" "Deposits" "Borrowing" ) ///
ar2  br  label nocon  nonotes s( FE N r2_a, labels(`"Bank F.E."' `"Observations"' `"Adj. \(R^{2}\)"'  ) ///
fmt(0 %9.0fc 3))  width(\hsize) keep(D.FEDFUNDS)


local contr ln_assets ln_loans ln_securities ln_deposits ln_borrowing2  ln_equity r_dep r_loan intincnet_assets nonintexp_assets nonintinc_assets 
local macro  
local horizon=`num'*4
local list s_dep s_loan_ma net_nonintexp_assets // nonintexp_assets nonintinc_assets 
eststo clear
foreach var of local list {
eststo `var': reghdfe F`horizon'.`var'  `contr' `macro'  (D.FEDFUNDS=Surprise), ab(rssdid) cluster(rssdid tq)
estadd local FE "Yes"
}
 
esttab using "$resultfolder\Tables/LP_prices_`num'y.tex", ///
replace b(3) se(3) star(* 0.10 ** 0.05 *** 0.01  ) ///
	mtitle("Deposit spread" "Loan spread"     "Net non-int. exp." ) ///
ar2  br  label nocon  nonotes s( FE N r2_a, labels(`"Bank F.E."' `"Observations"' `"Adj. \(R^{2}\)"'  ) ///
fmt(0 %9.0fc 3))  width(\hsize) keep(D.FEDFUNDS)

}







* Reversal rate effects on lending

use "$folder\FDIC_SOD.dta", clear
collapse (sum) depsumbr, by( rssdid stcntybr year)
g asset = depsumbr
egen double sasset=sum(asset), by(stcntybr year)
bysort stcntybr year: g double shsqrd= (asset/sasset)^2 
egen double HHI=sum(shsqrd), by(stcntybr year)
collapse (mean) HHI [iw=depsumbr], by(rssdid year)
save hhi_rssdid, replace





use call_report_sample, clear
keep if  tq>=tq(1994q1) & tq<=tq(2017q4)
g year = year(dofq(tq))
merge m:1 rssdid year using "hhi_rssdid"
keep if _m==3
drop _m


foreach var of varlist assets loans securities deposits borrowing2 equity  {
g ln_`var'=ln(`var')
}


foreach var of varlist  r_dep r_loan FEDFUNDS s_dep s_loan {
replace `var'=`var'/100
}

tsset rssdid tq
g cons=1
set more off

g D_DFF = D.FEDFUNDS
g D_DFF_HHI = D.FEDFUNDS*HHI
g Low = FEDFUNDS<.02
g D_DFF_Low = D_DFF*HHI
g HHI_Low = Low*HHI


label var HHI_Low "HHI*Low"

local contr ln_assets ln_loans ln_securities ln_deposits ln_borrowing2  ln_equity r_dep r_loan intincnet_assets nonintexp_assets nonintinc_assets 
local macro  
local horizon 4
eststo clear
local var ln_equity
eststo : reghdfe F`horizon'.`var'   `macro' HHI  HHI_Low if tq>=tq(2000q1) & tq<=tq(2004q1) , ab(tq rssdid) 
estadd local control "No"
estadd local FE "Yes"
eststo : reghdfe F`horizon'.`var'  `contr' `macro' HHI  HHI_Low if tq>=tq(2000q1) & tq<=tq(2004q1) , ab(tq rssdid) 
estadd local control "Yes"
estadd local FE "Yes"
local var ln_loans
eststo : reghdfe F`horizon'.`var'  `macro' HHI  HHI_Low if tq>=tq(2000q1) & tq<=tq(2004q1) , ab(tq rssdid) 
estadd local control "No"
estadd local FE "Yes"
eststo : reghdfe F`horizon'.`var'  `contr' `macro' HHI  HHI_Low if tq>=tq(2000q1) & tq<=tq(2004q1) , ab(tq rssdid) 
estadd local control "Yes"
estadd local FE "Yes"
esttab using "$resultfolder\Tables/reversal_lending.tex", ///
replace b(3) se(3) star(* 0.10 ** 0.05 *** 0.01  ) ///
	mtitle( "Equity" "Equity" "Loan"   "Loan") ///
ar2    label nocon  nonotes s(control FE FE N r2_a, labels("Control" `"Bank F.E."' "Time F.E." `"Observations"' `"Adj. \(R^{2}\)"'  ) ///
fmt(0 0 0 %9.0fc 3))  width(\hsize) keep(HHI_Low)



 !rmdir "$working"  /s /q // delete files in temp folder

