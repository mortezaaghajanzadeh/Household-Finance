
// file path
global path ".."
global folder "$path\data" // Specify data directory.
global working "$path\temp" // Specify working directory.
global resultfolder "$path\output" // Specify result directory.


cd "$working"


clear
set more off



* Treasury
capture: erase __000000.txt
freduse DFF DGS1 DGS2 DGS5 DGS10, clear
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


* Kuttner shocks
import excel "$folder\daily FF surprises July 2008.xls", first clear 
g date=A
sort date
save kuttner_shock, replace

* Information shocks
insheet using "$folder\data_fig4_Jarocinski_Karadi.csv",  clear 
g tm= ym(year,period)
format tm %tm
keep ff4_hf sp500_hf mpshocksign cbishocksign tm
save karadi_shock, replace


use "$folder\FF3", clear

merge 1:1 date using "$folder\49industry.dta"
drop _merge

merge 1:1 date using "trea_yield.dta"
keep if _m==3
drop _merge

merge 1:1 date using "FFR_shock.dta"
keep if _m==3
drop _merge

merge 1:1 date using "kuttner_shock.dta"
drop _merge


merge 1:1 date using "$folder\fomc_date.dta"
drop _merge
replace isfomc=0 if isfomc==.

g tm = mofd(date)
merge m:1 tm using "karadi_shock.dta"
drop _merge

foreach var of varlist Agric-Other {
replace `var'=. if `var'==-99.99
}

g year=year(date)

sort date
keep if isfomc==1
keep if year(date)>=1994 & year(date)<=2017
reg Mkt Dy_1y Dterm_5_1, robust 
reg Mkt Dy_1y Dterm_5_1 if year>=2010, robust 
reg Mkt Dy_1y Dterm_5_1 if year<2010, robust 

reg Banks Dy_1y Dterm_5_1, robust 
reg Banks Dy_1y Dterm_5_1 if year>=2012, robust 
reg Banks Dy_1y Dterm_5_1 if year<2012, robust 

reg Mkt FFR_shock , robust 
reg Mkt FFR_shock  if year>=2010, robust 
reg Mkt FFR_shock  if year<2010, robust 

g Banks_exret=Banks - Mkt
g crisis=(year>=2007 &year<=2009|year>=2000 &year<=2001)
*g crisis=(year>=2007 &year<=2009)

 * karadi shock
 twoway  scatter Banks_exret mpshocksign if isfomc==1 & year>=1994 & crisis==0  & DFF<=2, mlabel(date)  ///
 || lfit Banks_exret mpshocksign if isfomc==1 & year>=1994 & crisis==0 & DFF<=2 ///
, ylab(, nogrid) graphregion(color(white))  ///
 ytitle("Banking industry excess return") xtitle("Policy shock on FOMC days") ///
      legend( off)
 capture: graph export "$resultfolder\Figures/bank_karadi_below.pdf", replace  

 * karadi shock	  
twoway  scatter Banks_exret mpshocksign if isfomc==1 & year>=1994 & crisis==0 & DFF>2, mlabel(date)  ///
 || lfit Banks_exret mpshocksign if isfomc==1 & year>=1994 & crisis==0 & DFF>2 ///
, ylab(, nogrid) graphregion(color(white))  ///
 ytitle("Banking industry excess return") xtitle("Policy shock on FOMC days") ///
      legend( off)
 capture: graph export "$resultfolder\Figures/bank_karadi_above.pdf", replace  

 


 * two year treasury yield
 twoway  scatter Banks_exret Dy_2y if isfomc==1 & year>=1994 & crisis==0 & DFF<=2, mlabel(date)  ///
 || lfit Banks_exret Dy_2y if isfomc==1 & year>=1994 & crisis==0 & DFF<=2 ///
, ylab(, nogrid) graphregion(color(white))  ///
 ytitle("Banking industry excess return") xtitle("Policy shock on FOMC days") ///
      legend( off)
 capture: graph export "$resultfolder\Figures/bank_treasury2y_below.pdf", replace  

 
  twoway  scatter Banks_exret Dy_2y if isfomc==1 & year>=1994 & crisis==0 & DFF>2, mlabel(date)  ///
 || lfit Banks_exret Dy_2y if isfomc==1 & year>=1994 & crisis==0 & DFF>2 ///
, ylab(, nogrid) graphregion(color(white))  ///
 ytitle("Banking industry excess return") xtitle("Policy shock on FOMC days") ///
      legend( off)
 capture: graph export "$resultfolder\Figures/bank_treasury2y_above.pdf", replace  

 
 
 
 * bar charts
 use "$folder\FF3", clear

merge 1:1 date using "$folder\49industry.dta"
drop _merge

merge 1:1 date using "trea_yield.dta"
keep if _m==3
drop _merge

merge 1:1 date using "FFR_shock.dta"
keep if _m==3
drop _merge

merge 1:1 date using "$folder\fomc_date.dta"
drop _merge
replace isfomc=0 if isfomc==.

foreach var of varlist Agric-Other {
replace `var'=. if `var'==-99.99
}

g year=year(date)

sort date
keep if year(date)>=1994 & year(date)<=2017

g crisis=(year>=2007 &year<=2009|year>=2000 &year<=2001)

keep if isfomc==1
foreach var of varlist Agric-Other {
winsor `var', g(temp) p(.01)
drop `var'
rename temp `var'
}
foreach var of varlist Agric-Other {
winsor `var', g(temp) p(.01)
drop `var'
rename temp `var'
}
foreach var of varlist Agric-Other {
reg `var' Mkt
g `var'_exret = `var' - _b[Mkt]*Mkt
}

save temp1, replace


use temp1, replace

foreach var of varlist Agric-Other Mkt {
reg `var' Dy_2y if isfomc==1 & year>=1994 & crisis==0 & DFF<=2 
g beta1_`var'=_b[Dy_2y]
g t1_`var'=_b[Dy_2y]/_se[Dy_2y] 
g se1_`var'=_se[Dy_2y]

}

foreach var of varlist Agric-Other Mkt {
reg `var' Dy_2y if isfomc==1 & year>=1994 & crisis==0 & DFF>2 
g beta2_`var'=_b[Dy_2y]
g t2_`var'=_b[Dy_2y]/_se[Dy_2y]
g se2_`var'=_se[Dy_2y]
}

save temp, replace

use temp, clear
keep if _n==1
keep *beta* t* se*
g id=1
reshape long beta1_ beta2_ t1_ t2_ se1_ se2_, i(id) j(Industry) string
egen Industry_n = group(Industry)
gsort beta1_
gen n_beta=_n 

graph hbar  beta1_, over(Industry, axis(off) sort(1)) blabel(group, color(black) ) ///
 ylab(-15(5)5, nogrid) graphregion(color(white)) ysize(10) xsize(6) ytitle("FFR<=2%")  
 capture: graph export "$resultfolder\Figures/industry_fomc_bar_lt2.pdf", replace  

graph hbar  beta2_, over(Industry, axis(off) sort(1)) blabel(group, color(black) ) ///
 ylab(-15(5)5, nogrid) graphregion(color(white)) ysize(10) xsize(6)  ytitle("FFR>2%")
 capture: graph export "$resultfolder\Figures/industry_fomc_bar_gt2.pdf", replace  




 
 
 	  
