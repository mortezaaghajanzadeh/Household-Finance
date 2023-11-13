global path ".."
global folder "$path\data" // Specify data directory.
global working "$path\temp" // Specify working directory.
global resultfolder "$path\output" // Specify result directory.
global matlabfolder "$path\blp\blp_matlab_data" // Specify working directory for matlab.


cd "$working"



*Macro
freduse CFBABSHNO /// corporate bond held by households (Billions)
HNOTSAQ027S /// Treasury securities held by households (Millions)
TFAABSHNO /// Total financial assets held by households (Billions)
CURRENCY /// currency in circulation (Billions)
MMMFFAQ027S /// MMF (Millions)
WRMFSL /// retail MMF (Billions)
TSDABSHNO /// Saving Deposits held by households (Billions)
NCBCDCA /// Checkable deposits and currency held by households (Billions )
DPSACBW027SBOG /// Deposits, All Commercial Banks (Billions )
NCBDBIQ027S  /// Nonfinancial corporate business; debt securities; liability, Level (Millions )
NCBLL /// Nonfinancial Corporate Business; Loans; Liability
TOTLL /// Total bank lending
CMDEBT /// HH debt (Billions)
BCNSDODNS /// corporate debt (Billions)
DFF DGS5 PI, clear

g year=year(daten)
g tq=qofd(daten)
keep if tq>=tq(1984q1) & tq<=tq(2017q4)

* change to billion
foreach var of varlist  HNOTSAQ027S MMMFFAQ027S NCBDBIQ027S NCBLL {
replace `var'=`var'/1000
}
collapse  CURRENCY TFAABSHNO WRMFSL MMMFFAQ027S HNOTSAQ027S CFBABSHNO DPSACBW027SBOG NCBCDCA TSDABSHNO NCBDBIQ027S NCBLL TOTLL DFF DGS5 CMDEBT BCNSDODNS PI, by(year)
g ShortBond=HNOTSAQ027S*.2+MMMFFAQ027S
g LongBond=CFBABSHNO+HNOTSAQ027S*.8
g TotalBond=HNOTSAQ027S + CFBABSHNO + MMMFFAQ027S


g cash_share=CURRENCY/(CURRENCY+DPSACBW027SBOG+HNOTSAQ027S*.2+MMMFFAQ027S)
g deposit_share = DPSACBW027SBOG/(CURRENCY+DPSACBW027SBOG+HNOTSAQ027S*.2+MMMFFAQ027S)
g bankloan_share=TOTLL/(CMDEBT+BCNSDODNS)

sum cash_share deposit_share bankloan_share DFF if year>=1994 & year<=2017
tsset year
save  macro_year, replace





* bank

use "$folder\callreport_ffiec_merged.dta", clear 
* asset duration
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

keep rssdid dateq assets liabilities equity intexpdep domdepservicecharges ///
deposits  transdep 	savdep foreigndep timedepge100k timedeple100k /// 
intexpdep intexptransdep intexpsavdep intexpfordep intexptimedepge100k intexptimedeple100k ///
loans ciloans reloans intincreloans intincciloans intincloans loanleaselossprovision cash securities fedfundsrepoliab otherborrowedmoney ///
salaries exponpremises  netinc dividendoncommonstock loansleases* resloans_* ///
nonintexp numemployees nonintinc assets_dur
g tq=qofd(dateq)
format tq %tq 
g year=year(dateq)
save callreports_small_1, replace

use callreports_small_1, clear
keep if year>=1994 & year<=2017
bysort rssdid year: g N=_N
drop if N<4

sort rssdid tq
foreach var of varlist intexpdep intexptransdep intexpsavdep intexpfordep domdepservicecharges intexptimedepge100k intexptimedeple100k ///
intincloans intincciloans intincreloans salaries exponpremises loanleaselossprovision nonintexp  {
bysort rssdid year: g cum_`var'=sum(`var')
replace `var'=cum_`var' if cum_`var'!=.
} //accumulative interest income in a year

keep if quarter(dofq(tq))==4
g dep_total=deposits+foreigndep //total deposits
g r_depo=intexpdep/dep_total[_n]*100  // deposit rates
g r_tran=intexptransdep/transdep[_n]*100  // transaction deposit rates
g r_savi=intexpsavdep/savdep[_n]*100  // saving deposit rates
g r_fore=intexpfordep/foreigndep[_n]*100  // foreigh deposit rates
g r_timg=intexptimedepge100k/timedepge100k[_n]*100  // time deposit rates gr than 10k
g r_timl=intexptimedeple100k/timedeple100k[_n]*100  // time deposit rates less than 10k
g r_fees=domdepservicecharges/dep_total[_n]*100 // service fee for deposits
g r_loan=intincloans/loans[_n]*100 // loan rates
g r_ciloan=intincciloans/ciloans[_n]*100 // ciloan rates
g r_reloans=intincreloans/reloans[_n]*100 // ciloan rates

g borrowing=fedfundsrepoliab+otherborrowedmoney //borrowing
g net_borrowing=fedfundsrepoliab+otherborrowedmoney-securities-cash // net borrowing 
g security_cash=securities+cash
g salaries_at=salaries/assets[_n-1]*100
g exponpremises_at=exponpremises/assets[_n-1]*100
g loss_prov=loanleaselossprovision/loans[_n-1]*100
g nonintexp_at=nonintexp/assets[_n-1]*100

drop if missing(r_dep)
drop if missing(r_loan)
drop if missing(salaries_at)
drop if missing(exponpremises_at)
drop if missing(loss_prov)
drop if missing(nonintexp_at)

foreach var of varlist  r_depo r_tran r_savi r_fore r_timg r_timl r_loan r_ciloan r_reloans r_fee salaries_at exponpremises_at loss_prov nonintexp_at {

		winsor `var' , gen(temp) p(.01)
		replace `var' = temp 
		drop temp

   }

   
sum r_depo r_tran r_savi r_fore r_timg r_timl r_loan r_ciloan r_reloan r_fee salaries_at exponpremises_at loss_prov nonintexp_at // the deposit rates are too high
save callreports_small_2, replace


* numbers of banks in each market in the raw data
use callreports_small_2, clear
unique rssdid
g one=1
collapse (sum) one, by(tq)
sum one

* branches
use "$folder\FDIC_SOD.dta", clear
g numbranch=1
collapse (sum) numbranch, by(rssdid year)
bysort year: egen total_numbranch=sum(numbranch)
g sharebranch=numbranch/total_numbranch
drop total_numbranch
save numbranch, replace

***************************
* merges
use callreports_small_2, clear
merge m:1 rssdid year using numbranch
keep if _m==3
drop _m
g numemp_numbran=numemployees/numbranch

foreach var of varlist  r_dep  sharebranch numemp_numbran exponpremises_at salaries_at nonintexp_at {
drop if missing(`var')
}


foreach var of varlist  numemp_numbran {
winsor `var', g(temp) p(.005)
drop `var'
rename temp `var'
}

//scale variables
replace numemp_numbran=numemp_numbran/100

* size 
gsort tq  - assets
bysort tq: g asset_rank = _n/_N
g big = (asset_rank<=.1)


save bank_demand, replace 






******************************
* Restrict sample
******************************

use bank_demand, clear

merge m:1 year using macro_year
keep if _m==3
drop _m

g dep_share = dep_total/(HNOTSAQ027S*.2+MMMFFAQ027S+DPSACBW027SBOG+CURRENCY*.3)/10^6
g loan_share=loans/(BCNSDODNS+CMDEBT)/10^6

g outsidebank=0

bysort rssdid: egen min_dep_share=min(dep_share)
sum min_dep_share, d

replace outsidebank=1 if min_dep_share<=10^-5 // keep banks larger than 1 bps
tab outsidebank
bysort rssdid: egen min_loan_share=min(loan_share)
replace outsidebank=1 if min_loan_share<=10^-6 // keep banks larger than 1 bps
tab outsidebank

bysort rssdid: egen mean_numbranch=mean(numbranch)
replace outsidebank=1 if mean_numbranch <= 10

drop N
bysort year: g N=_N
tab N
drop N
bysort rssdid: g T=_N 
tab T
replace outsidebank=1 if T<=4 // a bank at least has 4 year of observations

replace rssdid = -1 if outsidebank==1 
save temp, replace

// share of tiny banks
use temp, clear
collapse (sum) dep_share, by(tq outside)
collapse (mean) dep_share, by(outside)
display dep_share[1]/(dep_share[1]+dep_share[2])

// combines tiny banks to one bank
use temp, clear
keep if  outsidebank==1
collapse (mean) r_dep r_fee r_loan numemp_numbran exponpremises_at salaries_at ///
[iw=dep_share],  by(rssdid year outsidebank)
save outsidebank_1, replace

use temp, clear
keep if  outsidebank==1
collapse numbranch dep_share loan_share , by(rssdid year outsidebank)
merge 1:1 rssdid year using outsidebank_1
drop _m
save outsidebank, replace

// append to the other banks
use temp, clear
keep if  outsidebank==0
append using outsidebank
tab year
// seems non-linearity at the extreme
 foreach var of varlist   salaries_at exponpremises_at   {
		winsor `var' , gen(temp) p(.1)
		replace `var' = temp 
		drop temp
   }
* create variables
foreach var of varlist  numbranch numemp_numbran {
g log_`var'=log(`var')
}

save bank_demand_small_sample, replace


******************************
* BLP deposits
******************************
use bank_demand_small_sample, clear  
// add cash
g id_ind=2
append using macro_year
keep if year>=1994 & year<=2017


replace id_ind=1 if id_ind==.
replace dep_share=CURRENCY*.3/(HNOTSAQ027S*.2+MMMFFAQ027S+DPSACBW027SBOG+CURRENCY*.3) if id_ind==1


g DFF_IV=DFF if id_ind==1
replace DFF_IV=0 if id_ind==2

foreach var of varlist rssdid r_dep numbranch numemp_numbran exponpremises_at salaries_at {
replace `var'=0 if id_ind==1
}

egen id= group(id_ind rssdid)
egen t= group(year)

*replace sharebranch=sharebranch*100

replace r_fee=0 if r_fee==.
g price=DFF-r_dep+r_fee
replace log_numbranch=0 if log_numbranch==.
replace log_numemp_numbran=0 if log_numemp_numbran==.


sort t id_ind id
order   t id id_ind dep_share price DFF log_numbranch log_numemp_numbran exponpremises_at salaries_at
outsheet  t id id_ind dep_share price DFF log_numbranch log_numemp_numbran exponpremises_at salaries_at DFF_IV using  "$matlabfolder\dep_1.csv", comma nolabel replace
/*        1 2  3      4           5     6   7             8      	         9                10 		11  */
save dep_demand_all_banks, replace

/* pre 2005*/
use dep_demand_all_banks, clear
keep if year<=2005 
drop id t
egen id= group(id_ind rssdid)
egen t= group(year)
sort t id_ind id
order   t id id_ind dep_share price DFF log_numbranch log_numemp_numbran exponpremises_at salaries_at
outsheet  t id id_ind dep_share price DFF log_numbranch log_numemp_numbran exponpremises_at salaries_at DFF_IV using  "$matlabfolder\dep_2.csv", comma nolabel replace
/*        1 2  3      4           5     6   7             8      	         9                10 		11  */

/* post 2005*/
use dep_demand_all_banks, clear
keep if year>2005 
drop id t
egen id= group(id_ind rssdid)
egen t= group(year)
sort t id_ind id
order   t id id_ind dep_share price DFF log_numbranch log_numemp_numbran exponpremises_at salaries_at
outsheet  t id id_ind dep_share price DFF log_numbranch log_numemp_numbran exponpremises_at salaries_at DFF_IV using  "$matlabfolder\dep_3.csv", comma nolabel replace
/*        1 2  3      4           5     6   7             8      	         9                10 		11  */


/* big banks*/
use dep_demand_all_banks, clear
drop if outsidebank==1
keep if big==1 | id_ind==1
drop id
egen id= group(id_ind rssdid)
sort t id_ind id
order   t id id_ind dep_share price DFF log_numbranch log_numemp_numbran exponpremises_at salaries_at
outsheet  t id id_ind dep_share price DFF log_numbranch log_numemp_numbran exponpremises_at salaries_at DFF_IV using  "$matlabfolder\dep_4.csv", comma nolabel replace
/*        1 2  3      4           5     6   7             8      	         9                10 		11  */


/* small banks*/
use dep_demand_all_banks, clear
drop if outsidebank==1
keep if big==0 | id_ind==1
drop id
egen id= group(id_ind rssdid)
sort t id_ind id
order   t id id_ind dep_share price DFF log_numbranch log_numemp_numbran exponpremises_at salaries_at
outsheet  t id id_ind dep_share price DFF log_numbranch log_numemp_numbran exponpremises_at salaries_at DFF_IV using  "$matlabfolder\dep_5.csv", comma nolabel replace
/*        1 2  3      4           5     6   7             8      	         9                10 		11  */


******************
* Summary STATS
******************
use bank_demand_small_sample, clear
g id_ind=2

foreach var of varlist rssdid r_dep sharebranch numbranch numemp_numbran exponpremises_at salaries_at {
replace `var'=0 if id_ind==1
}

egen id= group(id_ind rssdid)
egen t= group(year)

replace sharebranch=sharebranch*100

replace log_numbranch=0 if log_numbranch==.
replace log_numemp_numbran=0 if log_numemp_numbran==.

g s_loan = r_loan-DGS5
g loan_to_dep = loans/dep_total
g loan_to_assets = loans/assets
g borrowing2=assets-dep_total-equity
g borrowing_dep_total=borrowing2/dep_total
g deposits_assets=dep_total/assets
g nonintexp_assets=nonintexp/assets // nonintexp_assets ratio: note that in this code, the income and expenses are already annualized
g nonintinc_assets=nonintinc/assets // non-interest income
g nonintexp_assets_net = nonintexp_assets-nonintinc_assets
g leverage=assets/equity

foreach var of varlist loan_to_dep-leverage {
winsor2 `var', cut(1 99) replace
}


replace numemp_numbran = numemp_numbran*100 //change the unit to 1 person 
replace dep_share=dep_share*100
replace loan_share=loan_share*100
replace nonintexp_assets_net=nonintexp_assets_net*100

lab var dep_share "Deposit shares"
lab var loan_share "Loan shares"
lab var r_dep "Deposit rates"
lab var r_loan "Loan rates"
lab var numbranch "No. of branches"
lab var numemp_numbran "No. of employees per branch"
lab var exponpremises_at "Expenses of fixed assets"
lab var salaries_at "Salaries"
lab var loan_to_dep "Loan-to-deposit ratio"
lab var loan_to_assets "Loan-to-asset ratio"
lab var borrowing_dep_total "Borrowing-to-deposit ratio"
lab var deposits_assets "Deposit-to-asset ratio"
lab var nonintexp_assets_net "Net noninterest expenses"
lab var leverage "Book leverage"
lab var assets_dur "Asset maturity"


bysort tq: egen total_assets = sum(assets)
g asset_share = assets/total_assets


estpost su  dep_share loan_share r_dep r_loan numbranch numemp_numbran exponpremises_at salaries_at ///
nonintexp_assets_net loan_to_dep borrowing_dep_total deposits_assets  leverage assets_dur, d
est store A

 esttab A using "$resultfolder\Tables\sum_stats.tex", replace ///
cells("mean(fmt(%9.3f)) sd(fmt(%9.3f))  p10(fmt(%9.3f)) p25(fmt(%9.3f)) p50(fmt(%9.3f)) p75(fmt(%9.3f)) p90(fmt(%9.3f)) ") ///
label booktabs nonum gaps noobs ///
collabels("mean" "sd" "p10" "p25" "p50" "p75" "p90") width(\hsize)


collapse (mean) dep_share loan_share r_dep r_loan numbranch numemp_numbran exponpremises_at salaries_at ///
nonintexp_assets_net loan_to_dep borrowing_dep_total deposits_assets  leverage assets_dur [iw=asset_share]



******************
* IV Robustness
******************
use "$folder\oesm_97_20_ma", clear
keep area area_title 
rename area msabr
rename area_title msanamb
destring msabr, replace

duplicates drop
save "msa_bls_raw", replace


use "msa_bls_raw", clear
sort msanamb
g msanamb_new=msanamb
order msabr msanamb_new
replace msanamb_new = subinstr(msanamb_new, " PMSA", "",.)
replace msanamb_new = subinstr(msanamb_new, " MSA", "",.)
replace msanamb_new = subinstr(msanamb_new, " Metropolitan Division", "",.)
replace msanamb_new = subinstr(msanamb_new, " NECTA Division", "",.)

*br msabr msanamb_new msanamb
duplicates drop msabr msanamb_new, force
bysort msanamb_new: g N=_N
keep if N>1
tab N
save "msa_bls_temp1", replace
 
 
use  "msa_bls_temp1", clear
keep if msabr<10^4 // msa code before 2005
rename msabr msabr_pre_2005
save   "msa_bls_pre_2005", replace

use  "msa_bls_temp1", clear
keep if msabr>=10^4 // msa code after 2005
rename msabr msabr_post_2005
duplicates drop msanamb_new, force
merge 1:1 msanamb_new using "msa_bls_pre_2005"
keep if _m==3
drop _m
keep msabr_post_2005 msabr_pre_2005
duplicates drop
save "link_pre_post_2005_msa", replace





use "$folder\FDIC_SOD.dta", clear
keep msabr msanamb 
duplicates drop
save "msa_sod", replace

use "msa_sod", clear
bysort msabr: keep if _n==1
merge 1:1 msabr using "$folder\msa_bls"
sort msanamb
keep if _m==3
keep msabr msanamb
save "link_msa_bls_sod", replace


use "$folder\oesm_97_20_ma", clear
tab occ_title if strpos( occ_title, "Loan")

replace occ_title="Loan Officers" if strpos( occ_title, "Loan")
keep if occ_title=="All Occupations" |occ_title=="Tellers" |occ_title=="Loan Officers" //|occ_title=="Loan Interviewers and Clerks"|occ_title=="Loan Officers and Counselors"| occ_title=="Loan and Credit Clerks"| occ_title=="Loan Interviewers"


keep area area_title a_mean h_mean year occ_title 
collapse a_mean h_mean, by(area area_title year occ_title )
sort year area occ_title


replace a_mean=. if a_mean==-99
replace h_mean=. if h_mean==-99

replace occ_title = subinstr(occ_title, " ", "", .)
reshape wide a_mean h_mean, i(area area_title year) j(occ_title) string
rename area msabr
rename area_title msanamb
destring msabr, replace

g msabr_pre_2005 = msabr
merge m:1 msabr_pre_2005 using "link_pre_post_2005_msa"
replace msabr = msabr_post_2005 if _m==3
drop _m

merge m:1 msabr using "link_msa_bls_sod"
keep if _m==3
drop _m
keep msabr year a_mean* h_mean*
save "msa_wage", replace



use "$folder\FDIC_SOD.dta", clear
merge m:1 msabr year using "msa_wage"
tab year _m
keep if _m==3
collapse (mean) a_mean* [iw=depsumbr], by(year rssdid)
save "bank_wage", replace




  * iv robust, deposit, all wage
use dep_demand_all_banks, clear
merge 1:1 year rssdid using "bank_wage"

foreach var of varlist salaries a_meanAllOccupations a_meanLoanOfficers a_meanTellers {
g ln_`var'=log(`var')
}


g log_dep_share= log(dep_share)
g yield = - price
egen time_industry=group(id_ind year)

label var  log_numbranch "Log number of branches"
label var  log_numemp_numbran "Log number of employees"
label var  yield "Yield sensitivity"

local iv a_meanTellers
keep if `iv'!=.
tab year

 
eststo clear
eststo:reghdfe log_dep_share   log_numbranch log_numemp_numbran (yield= `iv'), ab(time_industry rssdid) 
estadd local FE1 "Yes"
estadd local FE2 "Yes"
eststo:reghdfe log_dep_share   log_numbranch log_numemp_numbran (yield=salaries_at exponpremises_at ) , ab(time_industry rssdid) 
estadd local FE1 "Yes"
estadd local FE2 "Yes"
eststo:reghdfe log_dep_share   log_numbranch log_numemp_numbran (yield=salaries_at exponpremises_at  `iv') , ab(time_industry rssdid) 
estadd local FE1 "Yes"
estadd local FE2 "Yes"

esttab
esttab using  "$resultfolder\Tables\demand_estimation_wage_deposit_1997_2017.tex" ///
,  replace b(3) se(3) star(* 0.10 ** 0.05 *** 0.01  )   ///
ar2    label nocon  nonotes mtitle("Local wage" "Bank expenses"  "All") s(FE1 FE2  N r2_a, ///
label("Bank F.E." "Year-Sector F.E." "Observations" "Adj. $ R^2$") ///
 fmt(0 0  %9.0fc 3 )) width(\hsize)   
 
 

 
  * iv robust, loan, all wage
use bank_demand_small_sample, clear
g id_ind=2

foreach var of varlist rssdid r_dep sharebranch numbranch numemp_numbran exponpremises_at salaries_at {
replace `var'=0 if id_ind==1
}

egen id= group(id_ind rssdid)
egen t= group(year)

replace sharebranch=sharebranch*100

replace log_numbranch=0 if log_numbranch==.
replace log_numemp_numbran=0 if log_numemp_numbran==.

g s_loan = r_loan-DGS5
g log_loan_share= log(loan_share)
g yield = s_loan
egen time_industry=group(id_ind year)
merge 1:1 year rssdid using "bank_wage"

label var  log_numbranch "Log number of branches"
label var  log_numemp_numbran "Log number of employees"
label var  yield "Yield sensitivity"

local iv a_meanTellers
keep if `iv'!=.
tab year
winsor2 `iv' , cut(2.5 97.5) replace


set more off
eststo clear
eststo:reghdfe log_loan_share   log_numbranch log_numemp_numbran (yield= `iv'  ), ab(time_industry rssdid) 
estadd local FE1 "Yes"
estadd local FE2 "Yes"
eststo:reghdfe log_loan_share   log_numbranch log_numemp_numbran (yield= exponpremises_at salaries_at), ab(time_industry rssdid) 
estadd local FE1 "Yes"
estadd local FE2 "Yes"
eststo:reghdfe log_loan_share   log_numbranch log_numemp_numbran (yield=salaries_at exponpremises_at `iv') , ab(time_industry rssdid) 
estadd local FE1 "Yes"
estadd local FE2 "Yes"
esttab

esttab using  "$resultfolder\Tables\demand_estimation_wage_loan_1997_2017.tex" ///
,  replace b(3) se(3) star(* 0.10 ** 0.05 *** 0.01  )   ///
ar2    label nocon  nonotes mtitle("Local wage" "Bank expenses"  "All") s(FE1 FE2  N r2_a, ///
label("Bank F.E." "Year-Sector F.E." "Observations" "Adj. R-squared") ///
 fmt(0 0  %9.0fc 3 )) width(\hsize)   


 
!rmdir "$working"  /s /q

