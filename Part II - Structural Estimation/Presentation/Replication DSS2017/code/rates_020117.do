clear all
set matsize 10000
set maxvar  10000
set more off

cd "output"

local product "mm25k" 

use "..\data\all_products_q", clear

sort  uninumbr quar
tsset uninumbr quar

tab quar

gen d1_`product'=s1.apy_eoq_`product'
gen d1_12mcd10k=d1.apy_eoq_12mcd10k
replace d1_12mcd10k=d1_fftar-d1_12mcd10k

gen d1_fftar_avgherfdepcty=d1_fftar*avgherfdepcty

gen d1_fftar_rise=0
replace d1_fftar_rise=d1_fftar if d1_fftar>0 

gen zerolower=0
replace zerolower=1 if year>=2009
egen fipszero = group(fips zerolower)
egen certzero = group(cert zerolower)

replace d1_`product'=d1_fftar-d1_`product'
drop if d1_`product'==. | d1_fftar_avgherfdepcty==.


forvalues i = 1(1)5 {
	foreach fixedeffect in stquar fips fipszero uninumbr  {
		sort `fixedeffect' uninumbr quar
		by `fixedeffect': egen obs=count(uninumbr)
		drop if obs==1
		drop obs
	}
}

******************************
***** Table 1, Panel C ******
*****************************

gen highherf=.
sum avgherfdepcty, detail
replace highherf=1 if avgherfdepcty>=r(p50)
replace highherf=0 if avgherfdepcty<r(p50)
gen depsumbr_mil=depsumbr/1000


outreg2 if d1_mm25k~=. using  "..\tables\table1_panel_c.txt", replace sum(log) keep(depsumbr_mil d1_mm25k d1_12mcd10k  avgherfdepcty) eqkeep(N mean sd)
outreg2 if  highherf==0 & d1_mm25k~=. using  "..\tables\table1_panel_c.txt", append sum(log) keep(depsumbr_mil d1_mm25k d1_12mcd10k  avgherfdepcty) eqkeep(N mean sd)
outreg2 if highherf==1 & d1_mm25k~=. using  "..\tables\table1_panel_c.txt", append sum(log) keep(depsumbr_mil d1_mm25k d1_12mcd10k  avgherfdepcty) eqkeep(N mean sd)



**************************************
***** Table 2, Panel A, Col 4-6 ******
**************************************

reghdfe d1_`product' d1_fftar_avgherfdepcty, absorb(stquar uninumbr fips fipszero) cluster(fips) tolerance(0.0001) maxiterations(100000)
outreg2 using "../tables/table2_panelA.txt", replace se bdec(3) bracket e(r2 r2_a df_r df_a )
reghdfe d1_`product' d1_fftar_avgherfdepcty, absorb(uninumbr fips  fipszero quar) cluster(fips) tolerance(0.0001) maxiterations(100000)
outreg2 using "../tables/table2_panelA.txt", append se bdec(3) bracket e(r2 r2_a df_r df_a )
reghdfe d1_`product' d1_fftar_avgherfdepcty, absorb(fips  fipszero quar) cluster(fips) tolerance(0.0001) maxiterations(100000)
outreg2 using "../tables/table2_panelA.txt", append se bdec(3) bracket e(r2 r2_a df_r df_a )

keep if tot_fips_year>1

forvalues i = 1(1)5 {
	foreach fixedeffect in bankquar stquar fips fipszero uninumbr  {
		sort `fixedeffect' uninumbr quar
		by `fixedeffect': egen obs=count(uninumbr)
		drop if obs==1
		drop obs
	}
}


**************************************
***** Table 2, Panel A, Col 1-3 ******
**************************************

reghdfe d1_`product' d1_fftar_avgherfdepcty, absorb(bankquar stquar uninumbr fips fipszero) cluster(fips)  tolerance(0.0001) maxiterations(100000)
outreg2 using "../tables/table2_panelA.txt", append se bdec(3) bracket e(r2 r2_a df_r df_a )
reghdfe d1_`product' d1_fftar_avgherfdepcty , absorb(bankquar uninumbr fips fipszero) cluster(fips)  tolerance(0.0001) maxiterations(100000)
outreg2 using "../tables/table2_panelA.txt", append se bdec(3) bracket e(r2 r2_a df_r df_a )
reghdfe d1_`product' d1_fftar_avgherfdepcty, absorb(fips quar  fipszero) cluster(fips) tolerance(0.0001) maxiterations(100000)
outreg2 using "../tables/table2_panelA.txt", append se bdec(3) bracket e(r2 r2_a df_r df_a )



local product "12mcd10k"

use "..\data\all_products_q", clear

sort quar
merge m:1 quar using "..\data\rates12m"
drop if _merge==2
drop _merge

sort  uninumbr quar
tsset uninumbr quar

gen d1_swap12m_eoq=d1.swap12m_eoq

gen d1_`product'=s1.apy_eoq_`product'

gen d1_fftar_avgherfdepcty=avgherfdepcty*d1_tnote

gen zerolower=0
replace zerolower=1 if year>=2009
egen fipszero = group(fips zerolower)

replace d1_`product'=d1_fftar-d1_`product'

drop if d1_`product'==. | d1_fftar_avgherfdepcty==.

foreach fixedeffect in stquar fips fipszero uninumbr  {
	sort `fixedeffect' uninumbr quar
	by `fixedeffect': egen obs=count(uninumbr)
	drop if obs==1
	drop obs
}

**************************************
***** Table 2, Panel B, Col 4-6 ******
**************************************

reghdfe d1_`product' d1_fftar_avgherfdepcty, absorb(stquar uninumbr fips fipszero ) cluster(fips) tolerance(0.0001) maxiterations(100000)
outreg2 using "../tables/table2_panelB.txt", replace se bdec(3) bracket e(r2 r2_a df_r df_a )
reghdfe d1_`product' d1_fftar_avgherfdepcty, absorb(uninumbr fips quar fipszero ) cluster(fips) tolerance(0.0001) maxiterations(100000)
outreg2 using "../tables/table2_panelB.txt", append se bdec(3) bracket e(r2 r2_a df_r df_a )
reghdfe d1_`product' d1_fftar_avgherfdepcty, absorb(fips quar fipszero) cluster(fips) tolerance(0.0001) maxiterations(100000)
outreg2 using "../tables/table2_panelB.txt", append se bdec(3) bracket e(r2 r2_a df_r df_a )

keep if tot_fips_year>1

forvalues i = 1(1)5 {
	foreach fixedeffect in bankquar stquar fips fipszero uninumbr  {
		sort `fixedeffect' uninumbr quar
		by `fixedeffect': egen obs=count(uninumbr)
		drop if obs==1
		drop obs
	}
}

**************************************
***** Table 2, Panel B, Col 1-3 ******
**************************************

reghdfe d1_`product' d1_fftar_avgherfdepcty, absorb(bankquar stquar uninumbr fips fipszero) cluster(fips)  tolerance(0.0001) maxiterations(100000)
outreg2 using "../tables/table2_panelB.txt", append se bdec(3) bracket e(r2 r2_a df_r df_a )
reghdfe d1_`product' d1_fftar_avgherfdepcty , absorb(bankquar uninumbr fips  fipszero) cluster(fips)  tolerance(0.0001) maxiterations(100000)
outreg2 using "../tables/table2_panelB.txt", append se bdec(3) bracket e(r2 r2_a df_r df_a )
reghdfe d1_`product' d1_fftar_avgherfdepcty, absorb(fips quar  fipszero  ) cluster(fips) tolerance(0.0001) maxiterations(100000)
outreg2 using "../tables/table2_panelB.txt", append se bdec(3) bracket e(r2 r2_a df_r df_a )





****************
** TABLE 4  ****
****************

local product mm25k 

use "..\data\all_products_q", clear

sort  uninumbr quar
tsset uninumbr quar

gen d1_`product'=d1.apy_eoq_`product'

gen d1_ff_exp3_avgherfdepcty=d1_ff_exp3*avgherfdepcty
gen d1_ff_surp3_avgherfdepcty=d1_ff_surp3*avgherfdepcty

gen zerolower=0
replace zerolower=1 if year>=2009
egen fipszero = group(fips zerolower)

replace d1_`product'=d1_fftar-d1_`product'
drop if d1_`product'==. | d1_ff_exp3_avgherfdepcty==. | d1_ff_surp3_avgherfdepcty==.

forvalues i = 1(1)5 {
	foreach fixedeffect in stquar uninumbr fips fipszero  {
		sort `fixedeffect' uninumbr quar
		by `fixedeffect': egen obs=count(uninumbr)
		drop if obs==1
		drop obs
	}
}

reghdfe d1_`product' d1_ff_exp3_avgherfdepcty d1_ff_surp3_avgherfdepcty, absorb(stquar uninumbr fips fipszero) vce(cluster fips) tolerance(0.0001) maxiterations(100000)
outreg2 using "../tables/table4.txt", replace se bdec(3) bracket e(r2 r2_a df_r df_a )
reghdfe d1_`product' d1_ff_exp3_avgherfdepcty d1_ff_surp3_avgherfdepcty, absorb(uninumbr fips quar fipszero) vce(cluster fips) tolerance(0.0001) maxiterations(100000)
outreg2 using "../tables/table4.txt", append se bdec(3) bracket e(r2 r2_a df_r df_a )
reghdfe d1_`product' d1_ff_exp3_avgherfdepcty d1_ff_surp3_avgherfdepcty, absorb(fips quar fipszero) vce(cluster fips) tolerance(0.0001) maxiterations(100000)
outreg2 using "../tables/table4.txt", append bdec(3) bracket e(r2 r2_a df_r df_a )

keep if tot_fips_year>1


forvalues i = 1(1)5 {
	foreach fixedeffect in bankquar stquar uninumbr fips fipszero   {
		sort `fixedeffect' uninumbr quar
		by `fixedeffect': egen obs=count(uninumbr)
		drop if obs==1
		drop obs
	}
}

reghdfe d1_`product' d1_ff_exp3_avgherfdepcty d1_ff_surp3_avgherfdepcty , absorb(bankquar stquar uninumbr fips fipszero) vce(cluster fips) tolerance(0.0001) maxiterations(100000)
outreg2 using "../tables/table4.txt", append se bdec(3) bracket e(r2 r2_a df_r df_a )
reghdfe d1_`product' d1_ff_exp3_avgherfdepcty d1_ff_surp3_avgherfdepcty , absorb(bankquar uninumbr fips fipszero) vce(cluster fips) tolerance(0.0001) maxiterations(100000)
outreg2 using "../tables/table4.txt", append se bdec(3) bracket e(r2 r2_a df_r df_a )
reghdfe d1_`product' d1_ff_exp3_avgherfdepcty d1_ff_surp3_avgherfdepcty, absorb(fips quar fipszero) vce(cluster fips) tolerance(0.0001) maxiterations(100000)
outreg2 using "../tables/table4.txt", append se bdec(3) bracket e(r2 r2_a df_r df_a )



********************************
***** Table 5, Panel A      **
*******************************

local product "mm25k"

use "..\data\all_products_q", clear
sort fips
merge m:1 fips using "..\data\fips_char"
drop if _merge==2
tab _merge
drop _merge

sort  uninumbr quar
tsset uninumbr quar

gen d1_`product'=s1.apy_eoq_`product'
gen d1_fftar_avgherfdepcty=d1_fftar*avgherfdepcty
gen d1_fftar_age=d1_fftar*(pop_2000_over65/100)
gen d1_fftar_income=d1_fftar* lnmedic1999
gen d1_fftar_edu_hs=d1_fftar*(pop_2000_highschool/100)
gen d1_fftar_edu_ba=d1_fftar*(pop_2000_ba/100)
	 
gen zerolower=0
replace zerolower=1 if year>=2009
egen fipszero = group(fips zerolower)

replace d1_`product'=d1_fftar-d1_`product'
drop if d1_`product'==. | d1_fftar_avgherfdepcty==. | d1_fftar_age==. | d1_fftar_income==. | d1_fftar_edu_hs==. | d1_fftar_edu_ba==.


foreach fixedeffect in stquar fips fipszero uninumbr  {
	sort `fixedeffect' uninumbr quar
	by `fixedeffect': egen obs=count(uninumbr)
	drop if obs==1
	drop obs
}

reghdfe d1_`product' d1_fftar_avgherfdepcty d1_fftar_age, absorb(stquar uninumbr fips  fipszero) vce(cluster fips) tolerance(0.0001) maxiterations(100000)
outreg2 using "../tables/table5_panelA.txt", replace se bdec(3) bracket e(r2 r2_a df_r df_a )
reghdfe d1_`product' d1_fftar_avgherfdepcty d1_fftar_income, absorb(stquar uninumbr fips  fipszero) vce(cluster fips) tolerance(0.0001) maxiterations(100000)
outreg2 using "../tables/table5_panelA.txt", append se bdec(3) bracket e(r2 r2_a df_r df_a )
reghdfe d1_`product' d1_fftar_avgherfdepcty d1_fftar_edu_ba, absorb(stquar uninumbr fips  fipszero) vce(cluster fips) tolerance(0.0001) maxiterations(100000)
outreg2 using "../tables/table5_panelA.txt", append se bdec(3) bracket e(r2 r2_a df_r df_a )
reghdfe d1_`product' d1_fftar_avgherfdepcty d1_fftar_age d1_fftar_income d1_fftar_edu_ba, absorb(stquar uninumbr fips  fipszero) vce(cluster fips) tolerance(0.0001) maxiterations(100000)
outreg2 using "../tables/table5_panelA.txt", append se bdec(3) bracket e(r2 r2_a df_r df_a )




********************************
***** Figure 2, Panel A and B**
*******************************

use "..\data\all_products_q", clear
sort uninumbr quar
rename fips fips_original
by uninumbr: egen fips=mode(fips_original)
gen test=0 if fips~=. & fips_original~=.
replace test=1 if fips~=fips_original & fips~=. & fips_original~=.
tab test, missing
collapse (min) fips, by(uninumbr)
sort uninumbr
save "../temp/bridge_uninumbr_fips", replace


foreach product in "mm25k"  {
capture graph drop _all

use "..\data\all_products_q", clear 

tsset uninumbr quar
sort uninumbr quar

gen d1_`product'=s1.apy_eoq_`product'
replace d1_`product'=d1_fftar-d1_`product'

statsby _b _se, by(uninumbr) saving("../temp/branch_coeff_`product'_long", replace):  reg d1_`product' d1_fftar

use "../temp/branch_coeff_`product'_long.dta", clear
capture rename _b_d1_tnote b_original
capture rename _b_d1_fftar b_original
winsor b_original, p(0.01) gen(b_winsor)
keep uninumbr b_winsor b_original
save "../temp/deposit_beta", replace

sort uninumbr
merge 1:1 uninumbr using "../temp/bridge_uninumbr_fips"
drop if _merge==2
drop _merge

collapse (mean) b_original b_winsor, by(fips)
drop if b_original==0 | b_original==.
 
sort fips
merge m:1 fips using "..\data\avgherfdepcty"
keep if _merge==3
drop _merge

xtile size = avgherfdepcty, nq(10)
replace size=size*10-5
label var size "Percentile County Avg. Herfindahl"

sort size
by size: egen avg_bwinsor=mean(b_winsor)
by size: egen avg_boriginal=mean(b_original)


graph twoway (scatter avg_bwinsor size, ylabel(0.63 (0.02) 0.75) name(coeff_cd) legend(off))  (lfit avg_bwinsor size)
graph export "../figures/figure4_panel_a.pdf", replace
}



foreach product in "12mcd10k"  {
capture graph drop _all

use "..\data\all_products_q", clear 

tsset uninumbr quar
sort uninumbr quar

gen d1_`product'=s1.apy_eoq_`product'
replace d1_`product'=d1_fftar-d1_`product'

statsby _b _se, by(uninumbr) saving("../temp/branch_coeff_`product'_long", replace):  reg d1_`product' d1_tnote

use "../temp/branch_coeff_`product'_long.dta", clear
capture rename _b_d1_tnote b_original
capture rename _b_d1_fftar b_original
winsor b_original, p(0.01) gen(b_winsor)
keep uninumbr b_winsor b_original
save "../temp/deposit_beta", replace

sort uninumbr
merge 1:1 uninumbr using "../temp/bridge_uninumbr_fips"
drop if _merge==2
drop _merge

collapse (mean) b_original b_winsor, by(fips)
drop if b_original==0 | b_original==.
 
sort fips
merge m:1 fips using "..\data\avgherfdepcty"
keep if _merge==3
drop _merge

xtile size = avgherfdepcty, nq(10)
replace size=size*10-5
label var size "Percentile County Avg. Herfindahl"

sort size
by size: egen avg_bwinsor=mean(b_winsor)
by size: egen avg_boriginal=mean(b_original)


graph twoway (scatter avg_bwinsor size, name(coeff_cd) legend(off))  (lfit avg_bwinsor size)
graph export "../figures/figure4_panel_b.pdf", replace
}
