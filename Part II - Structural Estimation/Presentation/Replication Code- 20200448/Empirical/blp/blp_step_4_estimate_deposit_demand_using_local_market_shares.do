// file path
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
DFF DGS5, clear

g year=year(daten)
g tq=qofd(daten)
keep if tq>=tq(1984q1) & tq<=tq(2017q4)

* change to billion
foreach var of varlist  HNOTSAQ027S MMMFFAQ027S NCBDBIQ027S NCBLL {
replace `var'=`var'/1000
}
collapse  CURRENCY TFAABSHNO WRMFSL MMMFFAQ027S HNOTSAQ027S CFBABSHNO DPSACBW027SBOG NCBCDCA TSDABSHNO NCBDBIQ027S NCBLL TOTLL DFF DGS5 CMDEBT BCNSDODNS, by(year)
g ShortBond=HNOTSAQ027S*.2+MMMFFAQ027S
g LongBond=CFBABSHNO+HNOTSAQ027S*.8
g TotalBond=HNOTSAQ027S + CFBABSHNO + MMMFFAQ027S


g cash_share=CURRENCY/(CURRENCY+DPSACBW027SBOG+TotalBond)
g deposit_share = DPSACBW027SBOG/(CURRENCY+DPSACBW027SBOG+TotalBond)


g bond_share=NCBDBIQ027S/(NCBDBIQ027S+NCBLL)

g relative_market_size=(NCBDBIQ027S+NCBLL)/(CURRENCY+DPSACBW027SBOG+TotalBond)
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

* zip-county mapping
insheet using "$folder\ZIP_COUNTY_122017.csv", clear
keep zip county
destring zip, replace
save zip_county, replace

* total number of branches per bank
use "$folder\FDIC_SOD.dta", clear
g numbranch=1
collapse (sum) total_numbranch=numbranch, by(rssdid year )
save total_numbranch, replace


* CBs
use "$folder\FDIC_SOD.dta", clear
g numbranch=1
collapse (sum) depsumbr numbranch, by(rssdid year stcntybr)
merge m:1 year rssdid using total_numbranch // total number of branches nationally, for computing number of employee per branch
keep if _m==3
drop _m
merge m:1 rssdid year using callreports_small_2, ///
keepusing( r_depo r_fee r_loan salaries_at exponpremises_at numemployees)
keep if _m==3
drop _m
save cb_year, replace



* CB small: select banks with large enough market shares
use cb_year, replace
drop if depsumbr==0
bysort stcntybr year: egen sum_depsumbr=sum(depsumbr)
g share = depsumbr/sum_depsumbr
bysort stcntybr rssdid: egen min_share=min(share)
sum min_share, d

g outsidebank=(min_share<=1*10^-2) // keep banks larger than 1 bps
bysort stcntybr rssdid: g T=_N 
tab T
replace outsidebank=1 if T<=4 // a bank at least has 4 year of observations
replace rssdid = -1 if outsidebank==1 
save temp, replace

// combines tiny banks to one bank
use temp, clear
keep if  outsidebank==1
collapse (mean) r_depo r_fee r_loan salaries_at exponpremises_at numemployees total_numbranch ///
(sum) depsumbr numbranch  , by(rssdid stcntybr year outsidebank)
save outsidebank, replace

// append to the other banks
use temp, clear
keep if  outsidebank==0
append using outsidebank
bysort stcntybr year: g N=_N
tab N
keep rssdid stcntybr year outsidebank  r_depo r_fee r_loan salaries_at exponpremises_at numemployees total_numbranch ///
 depsumbr numbranch  
save cb_year_small, replace


use "$folder\county_info_BEA_2019", clear
split geoname, p(",")
drop if geoname2=="" // these are total values in a state
rename geoname1 county
rename geoname2 state
drop geoname3
keep if year>=1994 & year<=2017
bysort year: egen sum_personal_inc=sum(personal_inc)
g weight = personal_inc/sum_personal_inc
bysort geofips: egen mean_weight=mean(weight)
drop if missing(mean_weight)
set seed 1
*sample2 5, cluster(geofips)  
save county_sample, replace


//cash
use county_sample, clear
merge m:1 year using macro_year, nogenerate  keep(3)
g m = CURRENCY*10^6*mean_weight
g id_ind = 1
save county_cash, replace

//bond
use county_sample, clear
merge m:1 year using macro_year, nogenerate keep(3)
g m = (HNOTSAQ027S*.2+MMMFFAQ027S)*10^6 *mean_weight
g id_ind = 3
save county_bond, replace

//cb 
use county_sample, clear
merge m:1 year using macro_year, nogenerate assert(2 3) keep(3)
g stcntybr=geofips
destring stcntybr, replace
merge 1:m year stcntybr using cb_year
tab year _m
keep if _m==3

// combine small banks into 1
drop if depsumbr==0
bysort stcntybr year: egen sum_depsumbr=sum(depsumbr)
g share = depsumbr/sum_depsumbr
bysort stcntybr rssdid: egen min_share=min(share)
sum min_share, d
g outsidebank=(min_share<=.4*10^-1) // keep banks larger than 1 bps

bysort stcntybr rssdid: g T=_N 
tab T
replace outsidebank=1 if T<=4 // a bank at least has 4 year of observations
replace rssdid = -1 if outsidebank==1 
save temp, replace

// combines tiny banks to one bank
use temp, clear
keep if  outsidebank==1
collapse (mean) r_depo r_fee r_loan salaries_at exponpremises_at numemployees popu DFF total_numbranch ///
(sum) depsumbr numbranch  , by(rssdid geofips stcntybr year outsidebank)
save outsidebank, replace

// append to the other banks
use temp, clear
keep if  outsidebank==0
append using outsidebank
bysort stcntybr year: g N=_N
tab N
keep rssdid stcntybr geofips year outsidebank  r_depo r_fee r_loan salaries_at exponpremises_at ///
total_numbranch numemployees  depsumbr numbranch  popu DFF
g m = depsumbr
g id_ind = 2
save county_cb, replace






***************************
* merges
use county_cash, clear 
append using county_cb 
append using county_bond

bysort geofips year: egen sum_m = sum(m)
g dep_share = m/sum_m
drop if id_ind==3 // assign bond as outside good

g numemp_numbran=numemployees/total_numbranch // note that total_numbranch is the national total
g branch_density=numbranch/popu*1000000 // per 1 million people
replace r_fee=0 if r_fee==.

foreach var of varlist rssdid r_dep branch_density numemp_numbran exponpremises_at salaries_at {
replace `var'=0 if id_ind==1
}

g price = DFF -r_dep+r_fee 

tab year
bysort geofips year: g N=_N 
tab N
// seems non-linearity at the extreme
 foreach var of varlist   salaries_at exponpremises_at   {
		winsor `var' , gen(temp) p(.1)
		replace `var' = temp 
		drop temp
   }
sum numemp_numbran, d
winsor2 numemp_numbran, cut(1 99) replace
sum numemp_numbran, d

foreach var of varlist  numemp_numbran numbranch {
g log_`var'=log(`var')
replace log_`var'=0 if log_`var'==.
}

egen id= group(id_ind rssdid)
egen t= group(year geofips)
egen id_time= group(year)
egen id_county= group(geofips)

sort t id_ind id
order   t id id_ind dep_share price DFF branch_density log_numemp_numbran exponpremises_at salaries_at id_time id_county

save bank_demand_local_sample, replace


* reduced form results
use bank_demand_local_sample, clear
g log_dep_share= log(dep_share)
egen time_industry=group(id_ind year)
egen bank_county=group(id id_county)
g yield = - price

label var  log_numbranch "Log number of branches"
label var  log_numemp_numbran "Log number of employees"
label var  yield "Yield sensitivity"

set more off
eststo clear
local y log_dep_share
eststo:reghdfe `y'  log_numbranch log_numemp_numbran (yield=exponpremises_at salaries_at),  ab(id time_industry t )  vce(cluseter time)
estadd local FE1 "Yes"
estadd local FE2 "Yes"
estadd local FE3 "Yes"
esttab

esttab using  "$resultfolder\Tables\demand_local.tex", ///
 replace b(3) se(3) star(* 0.10 ** 0.05 *** 0.01  )   ///
ar2  br  label nocon  nonotes mtitle("Deposit") s(FE1 FE2 FE3 N r2_a, ///
label("Bank F.E." "Year-Sector F.E." "Year-County F.E." "Observations" "Adj. R-squared") ///
 fmt(0 0 0 %9.0fc 3 )) width(\hsize)   
 
 !rmdir "$working"  /s /q // delete files in temp folder
