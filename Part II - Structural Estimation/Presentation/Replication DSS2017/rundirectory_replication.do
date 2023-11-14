clear
set more off

**********************************
** NOTE: Set the root directory **
**********************************
cd "C:\Users\pschnabl\Dropbox\deposits\qje final submission\replication"


** deposit spreads  (public data)
do "code\analysis_cra_020117"
cd ..

** deposit spread  (public data)
do "code\depositspread_020117"
cd ..

** branch-level flow analysis  (public data)
do "code\flows_020117.do"
cd ..

** bank-level analysis (public data)
do "code\banks_020117.do"
cd ..

** branch-level rates (Ratewatch data)
* do "code\rates_020117.do"
* cd ..

** weekly rates results (Ratewatch data)
* do "code\rates_weekly_020117.do"
* cd ..
