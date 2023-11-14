clear all
set matsize 10000
set maxvar  10000
set more off

cd "output"

use "..\data\all_products_w.dta", clear
keep if year<=2013
keep yyyyww
duplicates drop
sort yyyyww
gen time=_n
save "..\temp\bridge_time_yyyyww", replace


local length=5   
local product mm25k

// Analysis:
set more off
clear all
local fmax = `length'
local lmax = `length'
local N = 1+`fmax'+`lmax'
matrix week = J(`N',1,-`fmax')
forvalues i=2(1)`N' {
	matrix week[`i',1] = week[`i'-1,1]+1
}
matlist week

use "..\data\all_products_w.dta", clear
drop week
keep if year<=2013

sort yyyyww
merge m:1 yyyyww using "..\temp\bridge_time_yyyyww"
drop if _merge==2
drop _merge

sort  uninumbr time
tsset uninumbr time

gen tar_d1_ff_a=d1.ff_tar*avgherfdepcty

gen d1_`product'=d1.apy_eoq_`product'
gen d1_fftar=d1.ff_tar
replace d1_`product'=d1_fftar-d1_`product'

sort uninumbr time

foreach typ in tar  {

forvalues i=`fmax'(-1)1{
	by uninumbr: gen `typ'_f`i'_d1_ff_a = `typ'_d1_ff_a[_n+`i']
}

forvalues i=1/`lmax' {
	by uninumbr: gen `typ'_l`i'_d1_ff_a = `typ'_d1_ff_a[_n-`i']
}
}

reghdfe  d1_mm25k  tar_f1_d1_ff_a tar_f2_d1_ff_a tar_f3_d1_ff_a tar_f4_d1_ff_a tar_f5_d1_ff_a  tar_d1_ff_a tar_l1_d1_ff_a tar_l2_d1_ff_a tar_l3_d1_ff_a tar_l4_d1_ff_a tar_l5_d1_ff_a , absorb(time fips) tolerance(0.0001) maxiterations(100000)

matrix result = J(`N',2,0)
matlist result

local cum_beta = ""
local counter = 0

forvalues i=`fmax'(-1)1{

	local counter = `counter' + 1
	
	local cum_beta= "`cum_beta'" + "tar_f`i'_d1_ff_a"
	lincom "`cum_beta'" 
	display "`counter'"
	
	matrix result[`counter',1]=r(estimate)
	matrix result[`counter',2]=r(se)
		
	local cum_beta= "`cum_beta'" + "+ "
	display "`cum_beta'"


}



local counter = `counter' + 1
local cum_beta= "`cum_beta'" + "tar_d1_ff_a"
display "`cum_beta'"
lincom "`cum_beta'"
 
matrix result[`counter',1]=r(estimate)
matrix result[`counter',2]=r(se)



forvalues i=1(1)`lmax'{

	local counter = `counter' + 1
	
	local cum_beta= "`cum_beta'" + "+ tar_l`i'_d1_ff_a"
	lincom "`cum_beta'" 
	display "`counter'"
	
	matrix result[`counter',1]=r(estimate)
	matrix result[`counter',2]=r(se)
	display "`cum_beta'"
	


}

matlist result



clear
svmat week
svmat result


tsset week
label var week "Week"

gen sd1=result1+1.96*result2
gen sd2=result1-1.96*result2


twoway (tsline result1, tline(0) recast(connected) lpattern(solid) color(red) legend(off) name(exp) saving(exp, replace)) (tsline sd1 sd2,  lpattern(dash dash) color(blue blue))
graph export "../figures/figure5.pdf", replace

