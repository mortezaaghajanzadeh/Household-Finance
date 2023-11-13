INSTRUCTIONS ON HOW TO OBTAIN RESULTS AND FIGURES IN THE PAPER
``Bank Market Power and Monetary Policy Transmission: Evidence from a Structural Estimation" by 
Yifei Wang, Toni M. Whited, Yufeng Wu, and Kairong Xiao


*****************
    Model
*****************
The model-related caculation is organzied as follows:

0. the main folder contains the code to solve the baseline model together with the functions.
   run main.m --- it solves the model and  generate the model moments reported in Table 3. 
              --- it also performs the R_square check mentioned in Section 2.E.

1. the subfolder "Script" contains the code used to generate the tables and figures:
   run code F3_Concentration.m---it generates the relationship between FFR and deposit/loan spread and quanity under various bank concerntration in Figure 3.
   run code F4_CompareSpread.m---it generates the actual versus model-predicted deposit/loan rate comparison plots in Figure 4.
   run code F5_PlotSolution.m --- it generates the plots in Figure 5
   run code F6_ImpulseResponse.m---it generates the impulse response plots for positive/negative FFR shocks in Figure 6.
   run code T4_Determinants.m --- it generates the results reported in Table 4 pertaining to how different frictios contribute to monetary policy transmission.
   run code T5_Responsiveness.m---it generates the responsiveness of additional price/quantities to the FFR in Table 5.
   run code T9_BigSmallBanks.m---it generates the moments for big and smalls reported in Table 9 and examines how financial frictions influence the monetary transmission among big and small banks.
   run code T10_EarlyLate.m --- it generates the reuslts reported in Table 10 and examines how 1) operating costs, 2) number of banks & deposit/loan sensitivity, and 3) interest rate environment influence monetary tranmission in early and late sample periods.

   separetely, Estimate_SE.m calculates the standard errors for parameter estimates and moments conditions using the influence function approach.
               FFR_ChargeoffProcess.m calibrates the FFR and bank chargeoff process--- we do not estimate the parameters governing these processes in our estimation.

2. the subfolder "Results" contain all the saved results from running the above code---all the results are named following the scripts used to generate them.

3. the subfolder "Data" contains the data used in the code.

4. the subfolder "BLP" contains the code used speficially to calcuate the depositor and borrowers' demand functions.

5. the subfolder "Library" contains some basic functions such as demeaning and winsoring data---these functions are not coded by the authors.



*****************
    EMPIRICAL 
*****************

The empirical analysis uses Stata (.do files) and Matlab (.m files). To obtain results from the empirical analysis, run the following files with the corresponding programs:

0. The subfolder "blp" contains demand estimation and data moment, please run the following programs step by step:
	a. blp_step_1_create_estimation_sample_for_matlab_deposits.do:		create deposit demand estimation samples for matlab, create summary statistic table (Table 2), create robustness tables with alternative instruments (Table G1, G2)
	b. blp_step_1_create_estimation_sample_for_matlab_loans.do:		create loan demand estimation samples for matlab
	c. blp_step_2_estimate_deposit_demand.m:				estimate deposit demand  
	d. blp_step_3_estimate_loan_demand.m:					estimate loan demand
	e. blp_step_4_estimate_deposit_demand_using_local_market_shares.do: 	estimate demand using local market shares (Table D1)
	f. blp_step_5_produce_other_data_moments.do:				produce other data moments for SMD (Table 3 panel A, Table 4 actual momments, Table 9, Table 10 panel A)
	g. blp_step_6_create_latex_table_for_other_data_moment.m		collect the demand estimation results and other moments (Table F1, Table RP 9-13)



1. The subfolder "local_projection contains local projection
	a. local_projection.do: 						create additional momments using local projection (Table 5) and reversal effect in loans and book equity (Table 8)


2. The subfolder "var" contains VAR result
	a. Create_VAR_Moment.m:							estimate sensitivity of aggregate credit and loans to FFR using var (Table 4 last two rows)


3. The subfolder "bank_stock_return" contains bank stock return regression, deposit and loan spread graphs
	a. bank_fomc.do: 							estimate bank stock return response to monetary policy shocks (Table 7, Table H1)
	b. industry_fomc.do: 							estimate bank industry stock return response to monetary policy shocks (Figure 7, Figure H1, H2)
	c. figures.do:								plot relationship between deposit (loan) spreads and FFR (Figure 1)


4. The subfolder "temp" is the working folder of the Stata program. It will be cleared up after running each do file

5. The subfolder "output" contains the tables and figures produced by the codel. "Results.tex" is a collection of the figures and latex tables created by the above programs. 

6. The subfolder "data" contains the raw data use in the empirical analysis. The orignal data sources are decribed in Section I of the paper.
	a. "callreport_ffiec_merged" and "data_extracted" are Call Report Data from FFIEC
	b. "FDIC_SOD" is the Summary of Deposits from FDIC
	c. "crsp_daily_stock_return" is the daily stock returns from CRSP
	d. "fedfundsfutures" is the fed fund future prices retrieved from Datastream
	e. "bank_repurchase_div_yield" is the bank dividend and share repurchase data retrieved from Compustat
	f. "49industry" and "FF3" are the returns 49 industry portfolios and three factor returns from Ken French's website 
	g. "oesm_97_20_ma" is the local wage data from BLS, "county_info_BEA_2019" is the county level income data from BEA
	h. "msa_bls" is the MSA code file from BLS
	i. "ZIP_COUNTY_122017" is the county code file from HUD
	j. "permco_rssdid_link" is the link file that links the bank ID in call report to CRSP
	k. "fomc_date" is the FOMC date from the Federal Reserve website
	l. "daily FF surprises July 2008" is the fed fund surprises from Kenneth N. Kuttner's website
	m. "data_fig4_Jarocinski_Karadi" is the monetary policy shocks from Jarocinski and Karadi (2020)


Additional note: 
The Stata program uses the following packages "freduse" "reghdfe" "winsor2" "esttab". If one's Stata program does not have the above ado files. Please type in "ssc install freduse" to install the "freduse" package. The same applies to other packages.