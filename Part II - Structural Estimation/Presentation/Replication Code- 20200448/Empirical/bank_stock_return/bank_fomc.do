clear
set more off


// Kairong's file path
global path ".."
*global path "C:\Users\xiao\Dropbox\Monetary Transmission in the Banking System\Code\replication_packet"

global folder "$path\data" // Specify data directory.
global working "$path\temp" // Specify working directory.
global resultfolder "$path\output" // Specify result directory.


cd "$working"


* Treasury
capture: erase __000000.txt
freduse DFF FEDFUNDS DGS1 DGS2 DGS5 DGS10, clear
sort date
drop date
rename daten date
sort date
g Dy_10y=DGS10[_n]-DGS10[_n-1]
g Dy_5y=DGS5[_n]-DGS5[_n-1]
g Dy_2y=DGS2[_n]-DGS2[_n-1]
g Dy_1y=DGS1[_n]-DGS1[_n-1]

g Dterm_10_1=Dy_10y-Dy_1y
g Dterm_10_2=Dy_10y-Dy_2y
g Dterm_5_1=Dy_5y-Dy_1y
g Dterm_5_2=Dy_5y-Dy_2y

save trea_yield, replace

* FFR Future
import excel "$folder\fedfundsfutures.xlsx", first clear cellrange(A2)
g date=Code
sort date
g FFR_shock=-(CAFCS00-CAFCS00[_n-1])
keep date FFR_shock
save FFR_shock, replace

* HHI
use "$folder\FDIC_SOD", clear
g entity=rssdid
replace entity=rssdhcr if rssdhcr!=0
collapse (sum) depsumbr, by( cntynumb year entity)
egen double sum_depsumbr=sum(depsumbr), by(cntynumb year)
bysort cntynumb year: g double sum_depsumbr_sqrd= (depsumbr/sum_depsumbr)^2 
egen double HHI=sum(sum_depsumbr_sqrd), by(cntynumb year)
collapse (mean) HHI [iw=depsumbr], by(entity year)
save hhi, replace

* link bank to bank holding company
use rssdid rssdhcr year using "$folder\FDIC_SOD", clear
g entity=rssdid
replace entity=rssdhcr if rssdhcr!=0
duplicates drop 
save rssdid_entity, replace

* duration mismatch
use "$folder\data_extracted.dta", clear

//replace . with 0 so that calculation could be allowed later
foreach var of varlist _all{
	capture: replace `var'=0 if `var'==.
}
		


keep if qofd(date)>=tq(1997q2)
//Assets
/*
q -> dollar value
t -> duration
i -> ID
j -> category of assets or liabilities
*/
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




//Liabilities
gen qij_sumj=timedeple100k_less_3m+timedeple100k_3m_1y+timedeple100k_1y_3y+timedeple100k_over_3y ///
+timedepge100k_less_3m+timedepge100k_3m_1y+timedepge100k_1y_3y+timedepge100k_over_3y ///
+transdep+savdep+fedfundsrepoliab+subordinateddebt


gen tqij_sumj=0.125*timedeple100k_less_3m+0.625*timedeple100k_3m_1y+2*timedeple100k_1y_3y+5*timedeple100k_over_3y ///
+0.125*timedepge100k_less_3m+0.625*timedepge100k_3m_1y+2*timedepge100k_1y_3y+5*timedepge100k_over_3y ///
+0*transdep+0*savdep+0*fedfundsrepoliab+5*subordinateddebt


gen liabilities_dur=tqij_sumj/qij_sumj

drop qij_sumj-tqij_sumj
g tq=qofd(date)

capture: drop _m
g rssdid = rssd
merge m:1 rssdid year using rssdid_entity
keep if _m==3
drop _m
collapse assets_dur liabilities_dur [iw=assets], by(entity tq)
save bank_duration, replace





* link to stock price
insheet using "$folder\permco_rssdid_link.csv", clear
drop if dt_end<19940101

foreach var of varlist dt_end dt_start {
g year=int(`var'/10000)
g month = floor(`var'/100)-100*year
g quarter= floor(month/4)+1
g day = `var'- 100*(floor(`var'/100))
replace day=day-1 if day>31 
replace `var'= mdy(month,day,year)
format `var' %td
drop month quarter day year
}

save permco_rssdid_link, replace

* stock return
use "$folder\crsp_daily_stock_return", clear
joinby permco using permco_rssdid_link
replace dt_end=td(31dec2017) if dt_end==td(31dec2016)
keep if date>=dt_start & date<=dt_end
duplicates drop permno date, force
drop dt_*
save bank_level_ret, replace
 
* 
use bank_level_ret, clear
g year=year(date)
g tq=qofd(date)

merge m:1 date using "$folder\fomc_date.dta"
drop _merge
replace isfomc=0 if isfomc==.
keep if  isfomc==1

merge m:1 entity tq using "bank_duration"
keep if _m==3|_m==1
drop _merge

 
merge m:1 entity year using "hhi"
keep if _m==3
drop _merge


merge m:1 date using "trea_yield.dta"
keep if _m==3
drop _merge

merge m:1 date using "FFR_shock.dta"
keep if _m==3
drop _merge

merge m:1 date using  "$folder\FF3.dta"
keep if _m==3
drop _merge

replace ret = ret*100
g exret = ret-Mkt 

replace HHI=HHI*100

g duration_gap=assets_dur-liabilities_dur



bysort entity: egen mean_duration_gap=mean(duration_gap)
bysort entity: egen mean_HHI=mean(HHI)
sum mean_HHI, d
g High_HHI = mean_HHI>`r(p50)'
sum mean_duration_gap, d
g High_dura = mean_duration_gap>`r(p50)'

foreach var of varlist Dy_1y Dy_2y Dy_10y FFR_shock  Dterm_10_1 {
g duration_`var'=duration_gap*`var'
}
foreach var of varlist Dy_1y Dy_2y Dy_10y FFR_shock Dterm_10_1  {
g HHI_`var'=HHI*`var'
}
foreach var of varlist Dy_1y Dy_2y Dy_10y FFR_shock Dterm_10_1  {
g High_HHI_`var'=High_HHI*`var'
}

save bank_level_regression, replace


* regresssion

use bank_level_regression, clear
g crisis=(year>=2007 &year<=2009|year>=2000 &year<=2001)

keep if year(date)>=1994 & year(date)<=2017
keep if crisis==0

foreach var of varlist ret exret Mkt Dy_2y  Dterm_5_2 HHI_Dy_2y Dy_1y  Dterm_5_1 HHI_Dy_1y {
winsor `var', g(temp) p(.025)
drop `var'
rename temp `var'
}

g Low = (DFF<=2)
foreach var of varlist Mkt Dy_2y   Dterm_5_2 HHI_Dy_2y Dy_1y  Dterm_5_1 HHI_Dy_1y  {
g `var'_Low = `var'*Low
}


label var Dy_1y "Policy shock"
label var HHI "HHI"
label var HHI_Dy_1y "HHI*Policy shock"
label var Dterm_5_1 "$\Delta$ Term spread"
label var Dy_1y_Low "Low*Policy shock"
label var HHI_Dy_1y_Low "Low*HHI*Policy shock"

label var Mkt "Market return"
label var Dy_2y "Policy shock"
label var HHI_Dy_2y "HHI*Policy shock"
label var Dy_2y_Low "Low*Policy shock"
label var HHI_Dy_2y_Low "Low*HHI*Policy shock"

 * two year in one equation
local y ret
local x Dy_2y  Mkt Dterm_5_2
local x1 Dy_2y_Low  Mkt_Low Dterm_5_2_Low Low
local contr  HHI_Dy_2y  
local contr1 HHI_Dy_2y_Low Low
eststo clear
quietly: eststo: reg  `y' `x'  if DFF>2,  cluster(date)
estadd local FE "Yes"
quietly: eststo: reg  `y' `x'  if DFF<=2,  cluster(date)
estadd local FE "Yes"
quietly: eststo: reg  `y' `x' `x1' ,  cluster(date)
estadd local FE "Yes"
quietly: eststo: reg  `y' `x' `contr'  if DFF>2,  cluster(date)
estadd local FE "Yes"
quietly: eststo: reg  `y' `x' `contr'  if DFF<=2,  cluster(date)
estadd local FE "Yes"
quietly: eststo: reg  `y' `x' `x1' `contr' `contr1' ,  cluster(date)
estadd local FE "Yes"
esttab, order(Dy_2y Dy_2y_Low HHI_Dy_2y HHI_Dy_2y_Low)
esttab using "$resultfolder/Tables/bank_fomc_r1.tex",  ///
replace b(3) se(3) star(* 0.10 ** 0.05 *** 0.01  )   ///
mtitles("High" "Low" "All" "High"  "Low" "All")   ///
ar2    label nocon  nonotes  s( FE  N r2_a, ///
label("Control" "Observations" "Adj. $ R^2$") ///
 fmt(0 %9.0fc 3 )) width(\hsize)  keep(Dy_2y Dy_2y_Low HHI_Dy_2y HHI_Dy_2y_Low)  order(Dy_2y Dy_2y_Low HHI_Dy_2y HHI_Dy_2y_Low)




 * one year
 local y ret
local x Dy_1y  Mkt Dterm_5_1
local contr  HHI_Dy_1y  
eststo clear
quietly: eststo: reg  `y' `x'  if DFF>2,  cluster(date)
quietly: eststo: reg  `y' `x'  if DFF<=2,  cluster(date)
quietly: eststo: reg  `y' `x' `contr'  if DFF>2,  cluster(date)
quietly: eststo: reg  `y' `x' `contr'  if DFF<=2,  cluster(date)
esttab, order(Dy_1y HHI_Dy_1y Dterm_5_1 Mkt)
capture: esttab using "$resultfolder/Tables/bank_fomc_1y.tex", 		replace b(3) se(3) star(* 0.10 ** 0.05 *** 0.01  )   ///
mtitles("High FFR" "Low FFR"  "High FFR"  "Low FFR" )   ///
ar2  br  label nocon  nonotes  s(   N r2_a, ///
label("Observations" "Adj, R-squared") ///
 fmt( %9.0fc 3 )) width(\hsize)   order(Dy_1y HHI_Dy_1y Dterm_5_1 Mkt)

 


 
 
 

* effect of FFR shock on bank stocks
use bank_level_regression, clear
keep if year(date)>=1994 & year(date)<=2017
foreach var of varlist ret exret Mkt Dy_2y  Dterm_5_2 HHI_Dy_2y Dy_1y  Dterm_5_1 HHI_Dy_1y {
winsor `var', g(temp) p(.01)
drop `var'
rename temp `var'
}

reg ret FFR_shock

