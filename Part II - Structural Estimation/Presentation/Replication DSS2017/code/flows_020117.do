clear all
set matsize 10000
set maxvar  10000
set more off

cd "output"

****************************************
** TABLE 1, PANEL D -- No of Branches *
****************************************

use "..\data\sample_flows", clear
keep if year<=2013

collapse (count) branches=depsumbr, by(cert year)
collapse (mean) branches, by(cert)

sort cert
save "..\temp\cert_branches", replace

***************************
/**  TABLE 1, PANEL A  */
**************************

use "..\data\sample_flows", clear
keep if year<=2013

sort fips
merge m:1 fips using "..\data\avgherfdepcty"
drop if _merge==2
drop _merge
	
sort fips
merge m:1 fips using "..\data\fips_char"
drop if _merge==2
tab _merge
drop _merge	

keep if year==2000
drop if pop_2000_over65==. | pop_2000_ba==. | medinc_1999==. | avgherfdepcty==.

duplicates drop fips, force

gen highherf=.
sum avgherfdepcty, detail
replace highherf=1 if avgherfdepcty>=r(p50)
replace highherf=0 if avgherfdepcty<r(p50)

outreg2 using  "..\tables\table1_panel_a.txt", replace sum(log) keep(pop_2000_over65  pop_2000_ba medinc_1999 avgherfdepcty ) eqkeep(N mean sd)
outreg2 if  highherf==0  using  "..\tables\table1_panel_a.txt", append sum(log) keep(pop_2000_over65  pop_2000_ba medinc_1999 avgherfdepcty ) eqkeep(N mean sd)
outreg2 if highherf==1  using  "..\tables\table1_panel_a.txt", append sum(log) keep(pop_2000_over65  pop_2000_ba medinc_1999 avgherfdepcty ) eqkeep(N mean sd)


***************************
/**  PREPARE FLOW REGS  */
**************************
use "..\data\sample_flows", clear
keep if year<=2013

sort fips
merge m:1 fips using "..\data\avgherfdepcty"
drop if _merge==2
drop _merge
	
sort fips
merge m:1 fips using "..\data\fips_char"
drop if _merge==2
tab _merge
drop _merge	

* generate interactions with FF rate
gen d1_fftar_avgherfdepcty=d1_fftar*avgherfdepcty
gen d1_fftar_age=d1_fftar*(pop_2000_over65/100)
gen d1_fftar_income=d1_fftar*lnmedic1999
gen d1_fftar_edu_hs=d1_fftar*(pop_2000_highschool/100)
gen d1_fftar_edu_ba=d1_fftar*(pop_2000_ba/100)

gen zerolower=0
replace zerolower=1 if year>=2010
egen fipszero = group(fips zerolower)

drop if d1_lndep==. | d1_fftar_avgherfdepcty==.

forvalues i = 1(1)5 {
foreach fixedeffect in styear fips fipszero uninumbr year {
	sort `fixedeffect' uninumbr year
	by `fixedeffect': egen obs=count(uninumbr)
	drop if obs==1
	drop obs
}
}

***************************
/**  TABLE 1, PANEL B  */
**************************
gen highherf=.
sum avgherfdepcty, detail
replace highherf=1 if avgherfdepcty>=r(p50)
replace highherf=0 if avgherfdepcty<r(p50)

gen depsumbr_mil=depsumbr/1000

outreg2 using  "..\tables\table1_panel_b.txt", replace sum(log) keep(depsumbr_mil d1_lndep d1_fftar avgherfdepcty ) eqkeep(N mean sd)
outreg2 if  highherf==0  using  "..\tables\table1_panel_b.txt", append sum(log) keep(depsumbr_mil d1_lndep d1_fftar avgherfdepcty ) eqkeep(N mean sd)
outreg2 if highherf==1  using  "..\tables\table1_panel_b.txt", append sum(log) keep(depsumbr_mil d1_lndep d1_fftar avgherfdepcty ) eqkeep(N mean sd)

******************
/** TABLE 3 */
*****************

reghdfe d1_lndep d1_fftar_avgherfdepcty, absorb(styear uninumbr fips fipszero) vce(cluster fips) tolerance(0.0001) maxiterations(100000)
outreg2 using "..\tables\table3.txt", replace se bdec(3) bracket e(r2 r2_a df_r df_a )
reghdfe d1_lndep d1_fftar_avgherfdepcty, absorb(year uninumbr  fips fipszero) vce(cluster fips) tolerance(0.0001) maxiterations(100000)
outreg2 using "..\tables\table3.txt", append se bdec(3) bracket e(r2 r2_a df_r df_a )
reghdfe d1_lndep d1_fftar_avgherfdepcty, absorb(year fips fipszero) vce(cluster fips ) tolerance(0.0001) maxiterations(100000)
outreg2 using "..\tables\table3.txt", append se bdec(3) bracket e(r2 r2_a df_r df_a )

keep if tot_fips_year>1

forvalues i = 1(1)5 {
	foreach fixedeffect in bankyear styear fips fipszero uninumbr year {
		sort `fixedeffect' uninumbr year
		by `fixedeffect': egen obs=count(uninumbr)
		drop if obs==1
		drop obs
	}
}

reghdfe d1_lndep d1_fftar_avgherfdepcty, absorb(bankyear styear uninumbr fips fipszero) vce(cluster fips) tolerance(0.0001) maxiterations(100000)
outreg2 using "..\tables\table3.txt", append se bdec(3) bracket e(r2 r2_a df_r df_a )
reghdfe d1_lndep d1_fftar_avgherfdepcty, absorb(bankyear uninumbr fips  fipszero) vce(cluster fips) tolerance(0.0001) maxiterations(100000)
outreg2 using "..\tables\table3.txt", append se bdec(3) bracket e(r2 r2_a df_r df_a )
reghdfe d1_lndep d1_fftar_avgherfdepcty, absorb(year fips fipszero) vce(cluster fips) tolerance(0.0001) maxiterations(100000)
outreg2 using "..\tables\table3.txt", append se bdec(3) bracket e(r2 r2_a df_r df_a )


*************************
/** TABLE 5, PANEL B */
*************************

use "..\data\sample_flows", clear
keep if year<=2013

* average bank assets in 2010 dollars from call reports
sort fips
merge m:1 fips using "..\data\avgherfdepcty"
drop if _merge==2
drop _merge
	
sort fips
merge m:1 fips using "..\data\fips_char"
drop if _merge==2
tab _merge
drop _merge	

* generate interactions with FF rate
gen d1_fftar_avgherfdepcty=d1_fftar*avgherfdepcty
gen d1_fftar_age=d1_fftar*(pop_2000_over65/100)
gen d1_fftar_income=d1_fftar* lnmedic1999
gen d1_fftar_edu_ba=d1_fftar*(pop_2000_ba/100)

gen zerolower=0
replace zerolower=1 if year>=2010
egen fipszero = group(fips zerolower)

drop if d1_lndep==. | d1_fftar_avgherfdepcty==.
drop if d1_fftar_age==. | d1_fftar_income==. |  d1_fftar_edu_ba==.

foreach fixedeffect in styear fips fipszero uninumbr year {
	sort `fixedeffect' uninumbr year
	by `fixedeffect': egen obs=count(uninumbr)
	drop if obs==1
	drop obs
}



reghdfe d1_lndep d1_fftar_avgherfdepcty, absorb(styear uninumbr fips  fipszero) vce(cluster fips) tolerance(0.0001) maxiterations(100000)
outreg2 using "..\tables\table5_panelB.txt", replace se bdec(3) bracket e(r2 r2_a df_r df_a )
reghdfe d1_lndep d1_fftar_avgherfdepcty d1_fftar_age, absorb(styear uninumbr fips  fipszero) vce(cluster fips) tolerance(0.0001) maxiterations(100000)
outreg2 using "..\tables\table5_panelB.txt", append se bdec(3) bracket e(r2 r2_a df_r df_a )
reghdfe d1_lndep d1_fftar_avgherfdepcty d1_fftar_income, absorb(styear uninumbr fips  fipszero) vce(cluster fips) tolerance(0.0001) maxiterations(100000)
outreg2 using "..\tables\table5_panelB.txt", append se bdec(3) bracket e(r2 r2_a df_r df_a )
reghdfe d1_lndep d1_fftar_avgherfdepcty d1_fftar_edu_ba, absorb(styear uninumbr fips  fipszero) vce(cluster fips) tolerance(0.0001) maxiterations(100000)
outreg2 using "..\tables\table5_panelB.txt", append se bdec(3) bracket e(r2 r2_a df_r df_a )
reghdfe d1_lndep d1_fftar_avgherfdepcty d1_fftar_age d1_fftar_income d1_fftar_edu_ba, absorb(styear uninumbr fips  fipszero) vce(cluster fips) tolerance(0.0001) maxiterations(100000)
outreg2 using "..\tables\table5_panelB.txt", append se bdec(3) bracket e(r2 r2_a df_r df_a )

*******************************
******* FIGURE 2, Panel C  ***
******************************


use "..\data\sample_flows", clear
sort uninumbr year
drop if fips==.
by uninumbr: egen obs=count(year)
sort fips year uninumbr
by fips year: egen obsBranch=count(year)
collapse (min) obs obsBranch bkmo (mean) avgdeposit=depsumbr, by(uninumbr)
sort uninumbr
save "../temp/branch_char", replace


use "..\data\sample_flows", clear
sort uninumbr year
rename fips fips_original
by uninumbr: egen fips=mode(fips_original)
gen test=0 if fips~=. & fips_original~=.
replace test=1 if fips~=fips_original & fips~=. & fips_original~=.
tab test, missing
drop if test==1
collapse (min) fips, by(uninumbr)
sort uninumbr
save "../temp/bridge_uninumbr_fips", replace

use "..\data\sample_flows", clear
keep if year<=2013
graph drop _all

tsset uninumbr year
gen dlndep=d1.lndep*100

gen wdlndep = .
forvalues year = 1995(1)2013 {
	       winsor dlndep if year == `year', gen(temp) p(0.01) 
	       replace wdlndep = temp if year == `year'
	       drop temp
}

statsby _b _se, by(uninumbr) saving("../temp/branchlevel_coeff_flows_2013", replace):  reg wdlndep d1_fftar

graph drop _all

use "../temp/branchlevel_coeff_flows_2013", clear
capture rename _b_d1_fftar b_original
winsor b_original, p(0.05) gen(b_winsor)
keep uninumbr b_winsor b_original 

sort uninumbr
merge 1:1 uninumbr using "../temp/bridge_uninumbr_fips"
drop if _merge==2
drop _merge

sort uninumbr
merge 1:1 uninumbr using "../temp/branch_char"
drop if _merge==2
keep if obs>=20
drop _merge


collapse (median) b_original b_winsor, by(fips)
drop if b_original==0 | b_original==.

sort fips
merge m:1 fips using "..\data\avgherfdepcty"
keep if _merge==3
drop _merge

drop b_winsor
winsor b_original, gen(b_winsor) p(0.01)

xtile size = avgherfdepcty, nq(10)
replace size=size*10-5
label var size "Percentile County Avg. Herfindahl"

sort size
by size: egen avg_bwinsor=mean(b_winsor)
by size: egen avg_boriginal=mean(b_original)

graph twoway (scatter avg_bwinsor size, name(winsor) legend(off))  (lfit avg_bwinsor size)
graph export "../figures/figure4_panelc.pdf", replace

