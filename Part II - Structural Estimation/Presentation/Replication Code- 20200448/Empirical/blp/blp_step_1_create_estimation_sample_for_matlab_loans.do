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
keep rssdid dateq assets liabilities equity intexpdep domdepservicecharges ///
deposits  transdep 	savdep foreigndep timedepge100k timedeple100k /// 
intexpdep intexptransdep intexpsavdep intexpfordep intexptimedepge100k intexptimedeple100k ///
loans ciloans reloans intincreloans intincciloans intincloans loanleaselossprovision cash securities fedfundsrepoliab otherborrowedmoney ///
salaries exponpremises  netinc dividendoncommonstock loansleases* resloans_* ///
nonintexp numemployees 
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

g dep_share = dep_total/(HNOTSAQ027S*.2+MMMFFAQ027S+DPSACBW027SBOG+CURRENCY)/10^6
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
*drop if T<=4

replace rssdid = -1 if outsidebank==1 
save temp, replace

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
* BLP loan
******************************

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

* loan share
sort t id_ind id
order   t id id_ind loan_share s_loan DFF log_numbranch log_numemp_numbran exponpremises_at salaries_at 
outsheet  t id id_ind loan_share s_loan DFF log_numbranch log_numemp_numbran exponpremises_at salaries_at  using  "$matlabfolder\loan_1.csv", comma nolabel replace 
/*        1 2  3      4           5     6   7             8      	         9                10 		  */
save loan_demand_all_banks, replace



/* pre-2005*/
use loan_demand_all_banks, clear
keep if year<2005
drop id t
egen id= group(id_ind rssdid)
egen t= group(year)
sort t id_ind id
order   t id id_ind loan_share s_loan DFF log_numbranch log_numemp_numbran exponpremises_at salaries_at 
outsheet  t id id_ind loan_share s_loan DFF log_numbranch log_numemp_numbran exponpremises_at salaries_at using  "$matlabfolder\loan_2.csv", comma nolabel replace
/*        1 2  3      4           5     6   7             8      	         9                10 		11  */


/* post-2005*/
use loan_demand_all_banks, clear
keep if year>=2006
drop id t
egen id= group(id_ind rssdid)
egen t= group(year)
sort t id_ind id
order   t id id_ind loan_share s_loan DFF log_numbranch log_numemp_numbran exponpremises_at salaries_at 
outsheet  t id id_ind loan_share s_loan DFF log_numbranch log_numemp_numbran exponpremises_at salaries_at using  "$matlabfolder\loan_3.csv", comma nolabel replace
/*        1 2  3      4           5     6   7             8      	         9                10 		11  */


/* big banks*/
use loan_demand_all_banks, clear
keep if big==1 | id_ind==1
drop id t
egen id= group(id_ind rssdid)
egen t= group(year)
sort t id_ind id
order   t id id_ind loan_share s_loan DFF log_numbranch log_numemp_numbran exponpremises_at salaries_at 
outsheet  t id id_ind loan_share s_loan DFF log_numbranch log_numemp_numbran exponpremises_at salaries_at using  "$matlabfolder\loan_4.csv", comma nolabel replace
/*        1 2  3      4           5     6   7             8      	         9                10 		11  */


/* small banks*/
use loan_demand_all_banks, clear
keep if big==0 | id_ind==1
drop id t
egen id= group(id_ind rssdid)
egen t= group(year)
order   t id id_ind loan_share s_loan DFF log_numbranch log_numemp_numbran exponpremises_at salaries_at 
outsheet  t id id_ind loan_share s_loan DFF log_numbranch log_numemp_numbran exponpremises_at salaries_at using  "$matlabfolder\loan_5.csv", comma nolabel replace
/*        1 2  3      4           5     6   7             8      	         9                10 		11  */

