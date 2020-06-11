*********************************************************************************
*
*	Network Chart Production, Other
*	Collin Jones 
*	Purpose: This do-file is an archive of most other plots that we have produced
*		while writing the paper. This code is not meant to be run start-to-finish, but as
*		a general resource for recreating any of these older paper plots.
*
*********************************************************************************

use ../output/Contagion_Data, clear


foreach var of varlist BHFN6631 BHFN6636 BHCK2948 BHDM6631 BHDM6636{
	replace `var' = `var'/1000000
	}
label var BHFN6631 "Noninterest-bearing Foreign Depos"
label var BHFN6636 "Interest-bearing Foreign Depos"
label var BHDM6631 "Noninterest-bearing Domestic Depos"
label var BHDM6636 "Interest-bearing Domestic Depos"

/*
preserve
keep if mkmv_id == "064057"

line BHFN6631 BHFN6636 BHDM6631 BHDM6636 BHCK2948 qt_dt, ///
	xtitle(" ", size(medsmall)) xlabel($charts_start(4)$charts_end, angle(45) labsize(small)) ///
	ytitle(" ",  height(6)) subtitle("", position(11) size(small)) ///
	title(" ", pos(11) size(medlarge)) graphregion(color(white) style(none)) ///
	legend(r(6) region(style(none)) size(small))
graph export ../output/BoNY_depos.pdf, replace

restore

preserve
keep if mkmv_id == "857473"

line BHFN6631 BHFN6636 BHDM6631 BHDM6636 BHCK2948 qt_dt, ///
	xtitle(" ", size(medsmall)) xlabel($charts_start(4)$charts_end, angle(45) labsize(small)) ///
	ytitle(" ",  height(6)) subtitle("", position(11) size(small)) ///
	title(" ", pos(11) size(medlarge)) graphregion(color(white) style(none)) ///
	legend(r(6) region(style(none)) size(small))
graph export ../output/State_Street_depos.pdf, replace

restore
*/

foreach var of varlist liab_in_y15deposits ///
	liab_in_actual liab_in_possible liab_in RISKM370 y_15_in_compat BHCK2948 tot_uninsur_depos_y9c tot_insur_depos ///
	y9c_all_deposits_minusFDIC y9c_all_deposits_minusFDICfor y9c_all_deposits{

	replace `var' = `var'/1000000
}


local perc_plot BHDMB993 BHCKB995 BHCK3548 BHCK4062 BHCK3049 ///
	BHCK2309 BHCK2332 BHCK2333 BHCKC699 BHCKB557 BHCKB984 BHCK3547

local perc_plot_y15 y_15_intra_deposits RISKY833 RISKM365 RISKM366 RISKM367 RISKM368

foreach var of varlist `perc_plot' `perc_plot_y15' {

	replace `var' = `var'/1000000
}


label var y_15_in_compat "Y-15 Intra-Financial, Compatible"
label var y_15_intra_deposits "Y-15 Intra-Financial Deposits"
label var y9c_all_deposits_minusFDIC "Y-9C All Deposit Categories, Less FDIC"
label var y9c_all_deposits_minusFDICfor "Y-9C All Deposit Categories, Less FDIC, Less Foreign"
label var y9c_all_deposits "Y-9C All Deposits"
label var perc_deposits_in_y15 "Percent All Deposits In System"
label var tot_insur_depos "FDIC Insured"
label var liab_in_y15deposits "Inside Liabilities, with Deposits from Y-15"

label var BHDMB993 "Fed Funds in Domestic Offices"
label var BHCK3547 "Derivatives with Negative Fair Value"
label var BHCKB995 "Security Repos"
label var BHCK3548 "Trading Liabilities"
label var BHCK4062 "Subordinated Notes and Debts"
label var BHCK3049 "Net Deferred Taxes"
label var BHCK2309 "Commercial Paper"
label var BHCK2332 "Other borrowed money, <1year"
label var BHCK2333 "Other borrowed money, >1year"
label var BHCKC699 "Subordinated Notes and Debts, trusts"
label var BHCKB557 "Allowance for Losses on off-balance sheet exposure"
label var BHCKB984 "Other"
label var RISKY833 "Borrowings from Financial Institutions"
label var RISKM365 "Unused portion of credit lines from financial institutions"
label var RISKM366 "Net negative exposure of SFTS with financial institutions"
label var RISKM367 "Net negative Derivatives"
label var RISKM368 "Net negative Derivatives: Potential Future"


local perc_plot BHDMB993 BHCKB995 BHCKB557 BHCK3547

local perc_plot_y15 RISKY833 RISKM365 RISKM366 RISKM367 RISKM368

levelsof tkr if !missing(RISKM370), local(to_gen) clean


/*
foreach tic in `to_gen'{
	twoway (line liab_in liab_in_actual liab_in_possible BHCK2948 qt_dt if tkr == "`tic'" & qt_dt >= $charts_start & qt_dt <= $charts_end) || (scatter liab_in_y15deposits RISKM370 qt_dt if tkr == "`tic'" & qt_dt >= $charts_start & qt_dt <= $charts_end), ///
	xtitle(" ", size(medsmall)) xlabel($charts_start(4)$charts_end, angle(45) labsize(small)) ylabel(, valuelabel labsize(small)) ///
	ytitle(" ",  height(6)) subtitle("", position(11) size(small)) ///
	title(" ", pos(11) size(medlarge)) graphregion(color(white) style(none)) ///
	legend(r(13) region(style(none)) size(small))
	graph export ../output/y15_compare_`tic'.pdf, replace

	scatter `perc_plot' qt_dt if tkr == "`tic'" & !missing(RISKM370) & qt_dt >= $charts_start & qt_dt <= $charts_end, ///
	xtitle(" ", size(medsmall)) xlabel(211(4)$charts_end, angle(45) labsize(small)) ylabel(, valuelabel labsize(small)) ///
	ytitle(" ",  height(6)) subtitle("", position(11) size(small)) ///
	title("Y-9C", pos(11) size(medlarge)) graphregion(color(white) style(none)) ///
	legend(r(6) holes(5 6) region(style(none)) size(small)) name(me_breakdown_`tic', replace)
	
	scatter `perc_plot_y15' qt_dt if tkr == "`tic'" & !missing(RISKM370) & qt_dt >= $charts_start & qt_dt <= $charts_end, ///
	xtitle(" ", size(medsmall)) xlabel(211(4)$charts_end, angle(45) labsize(small)) ylabel(, valuelabel labsize(small)) ///
	ytitle(" ",  height(6)) subtitle("", position(11) size(small)) ///
	title("Y-15", pos(11) size(medlarge)) graphregion(color(white) style(none)) ///
	legend(r(6) region(style(none)) size(small)) name(y15_breakdown_`tic', replace)
	
	graph combine me_breakdown_`tic' y15_breakdown_`tic', r(1) xcommon ycommon
	graph export ../output/liab_in_breakdown_`tic'.pdf, replace
	


twoway (line  y9c_all_deposits_minusFDIC tot_insur_depos tot_uninsur_depos_y9c y9c_all_deposits_minusFDICfor qt_dt if tkr == "`tic'" & qt_dt >= $charts_start & qt_dt <= $charts_end) || ///
	(scatter  perc_deposits_in_y15 qt_dt if tkr == "`tic'", yaxis(2)) || (scatter y_15_intra_deposits qt_dt if tkr == "`tic'" & qt_dt >= $charts_start & qt_dt <= $charts_end, yaxis(1)), ///
	xtitle(" ", size(medsmall)) xlabel($charts_start(4)$charts_end, angle(45) labsize(small)) ylabel(, valuelabel labsize(small)) ///
	ytitle(" ",  height(6)) subtitle("", position(11) size(small)) ylabel(0 .25 .5 .75 1, axis(2)) ///
	title(" ", pos(11) size(medlarge)) graphregion(color(white) style(none)) ///
	legend(r(6) region(style(none)) size(small))
	graph export ../output/y15_compare_depos_`tic'.pdf, replace
	}
	

graph close
twoway (scatter beta qt_dt if tkr == "STT") || (scatter beta qt_dt if tkr == "BK"), ///
	xtitle(" ", size(medsmall)) xlabel($charts_start(4)$charts_end, angle(45) labsize(small)) ylabel(, valuelabel labsize(small)) ///
	ytitle(" ",  height(6)) subtitle("", position(11) size(small)) ///
	title(" ", pos(11) size(medlarge)) graphregion(color(white) style(none)) ///
	legend(order(1 "State Street" 2 "Bank of NY, Mellon") r(1) region(style(none)) size(small))
graph export ../output/STT_BK_betas.pdf, replace
*/
twoway (scatter beta_max_1 qt_dt if  max_bank == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end, mcolor(blue) mlabgap(*5) mlabsize(vsmall) mlabel(tkr) mlabpos(11) mlabangle(rvertical)) || ///
	(scatter beta_max_2 qt_dt if  max_bank == 2 & qt_dt >= $charts_start & qt_dt <= $charts_end, mcolor(red) mlabgap(*5) mlabsize(vsmall) mlabel(tkr) mlabpos(11) mlabangle(rvertical)) || ///
	(scatter beta_max_3 qt_dt if  max_bank == 3 & qt_dt >= $charts_start & qt_dt <= $charts_end, mcolor(green) mlabgap(*5) mlabsize(vsmall) mlabel(tkr) mlabpos(6) mlabangle(rvertical)) || ///
	(scatter beta_max_4 qt_dt if  max_bank == 4 & qt_dt >= $charts_start & qt_dt <= $charts_end, mcolor(maroon) mlabgap(*5) mlabsize(vsmall) mlabel(tkr) mlabpos(6) mlabangle(rvertical)), ///
	ytitle("",  height(6)) xtitle("", size(medsmall)) ylabel(.4 .5 .6 .7 .8 .9) ///
	title(" ", pos(11) size(medlarge)) graphregion(color(white) style(none)) /// 
	xlabel($charts_start(4)$charts_end, angle(45) labsize(medsmall)) ///
	legend(order(1 "Highest {&beta}" 2 "Second Highest {&beta}" 3 "Third Highest {&beta}" 4 "Fourth Highest {&beta}") r(1) region(style(none)) size(small)) ///
	subtitle("Inside Liabilities/Total Liabilities", position(11) size(small)) name(beta_first, replace)
graph export ../output/betas.pdf, as(pdf) replace

twoway (scatter beta_max_fullsample_1 qt_dt if  max_bank_fullsample == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end, mcolor(blue) mlabgap(*5) mlabsize(vsmall) mlabel(tkr) mlabpos(11) mlabangle(rvertical)) || ///
	(scatter beta_max_fullsample_2 qt_dt if  max_bank_fullsample == 2 & qt_dt >= $charts_start & qt_dt <= $charts_end, mcolor(red) mlabgap(*5) mlabsize(vsmall) mlabel(tkr) mlabpos(11) mlabangle(rvertical)) || ///
	(scatter beta_max_fullsample_3 qt_dt if  max_bank_fullsample == 3 & qt_dt >= $charts_start & qt_dt <= $charts_end, mcolor(green) mlabgap(*5) mlabsize(vsmall) mlabel(tkr) mlabpos(6) mlabangle(rvertical)) || ///
	(scatter beta_max_fullsample_4 qt_dt if  max_bank_fullsample == 4 & qt_dt >= $charts_start & qt_dt <= $charts_end, mcolor(maroon) mlabgap(*5) mlabsize(vsmall) mlabel(tkr) mlabpos(6) mlabangle(rvertical)), ///
	ytitle("",  height(6)) xtitle("", size(medsmall)) ylabel(.4 .5 .6 .7 .8 .9) ///
	title(" ", pos(11) size(medlarge)) graphregion(color(white) style(none)) /// 
	xlabel($charts_start(4)$charts_end, angle(45) labsize(medsmall)) ///
	legend(order(1 "Highest {&beta}" 2 "Second Highest {&beta}" 3 "Third Highest {&beta}" 4 "Fourth Highest {&beta}") r(1) region(style(none)) size(small)) ///
	subtitle("Inside Liabilities/Total Liabilities", position(11) size(small)) name(beta_first, replace)
graph export ../output/betas_full.pdf, as(pdf) replace


//bys qt_dt: egen beta_max_all_`num'	= max(beta) if max_bank_all == .
//replace max_bank_all = `num' if beta_max_all_`num'== beta

*Our "Main" plot - loss ratio implied by gamma = 20%, highest beta used
twoway (tsline l_R_20_1 if max_bank == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end)  , ///
	xtitle("", size(medsmall)) xlabel($charts_start(4)$charts_end, angle(45) labsize(medsmall)) ///
	ytitle("",  height(6)) ylabel(10 20 30 40) ///
	title("Network Vulnerability Index", pos(11) size(medlarge)) 
graph export ../output/Vulnerability_Index.pdf, as(pdf) replace
//graph export ../output/Vulnerability_Index.emf, as(emf) replace

*Loss ratio without bankruptcy costs, highest beta used
twoway (tsline l_R_nogamma  if max_bank == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end)  , ///
	xtitle("", size(medsmall)) xlabel($charts_start(4)$charts_end, angle(45) labsize(medsmall)) ///
	ytitle("",  height(6)) ///
	title("Network Vulnerability Index", pos(11) size(medlarge)) 
graph export ../output/Vulnerability_Index_nogamma.pdf, as(pdf) replace
//graph export ../output/Vulnerability_Index_nogamma.emf, as(emf) replace

*Loss ratio with gamma = 1%, highest beta used
twoway (tsline l_R_1_1 if max_bank == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end)  , ///
	xtitle("", size(medsmall)) xlabel($charts_start(4)$charts_end, angle(45) labsize(medsmall)) ///
	ytitle("",  height(6)) ///
	title("Network Vulnerability Index", pos(11) size(medlarge)) 
graph export ../output/Vulnerability_Index_gamma1perc.pdf, as(pdf) replace
//graph export ../output/Vulnerability_Index_gamma1perc.emf, as(emf) replace

*Loss ratio with gamma 50% of maximum possible, highest beta used
twoway (tsline l_R_maxratio_500 if max_bank == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end)  , ///
	xtitle("", size(medsmall)) xlabel($charts_start(4)$charts_end, angle(45) labsize(medsmall)) ///
	ytitle("",  height(6)) ///
	title("Network Vulnerability Index", pos(11) size(medlarge)) 
graph export ../output/Vulnerability_Index_gamma50ofmax.pdf, as(pdf) replace
//graph export ../output/Vulnerability_Index_gamma50ofmax.emf, as(emf) replace


**************************************************************************************************
*
*					Plots with Flow of Funds
*
***************************************************************************************************

local to_check GS MS BAC2 BSC LEHMQ
local code_plot
foreach tic in `to_check'{
	local code_plot `code_plot' test_perc_weights_`tic'
}


*Index with insurance companies and mutual funds included - % estimated from FFunds categories
twoway (tsline logl_R_20_1 index_w_insurance nvi_benchmark if max_bank == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end)  , ///
	legend(order(1 "Y-9C BHC" 2 "Insurance Companies Added" 3 "Insurance Companies, Dealers Added") r(2) region(style(none)) size(small)) ///
	ytitle("",  height(6)) xtitle("", size(medsmall)) ///
	title(" ", pos(11) size(medlarge)) graphregion(color(white) style(none)) /// 
	xlabel($charts_start(4)$charts_end, angle(45) labsize(medsmall))
graph export ../output/index_w_additions.pdf, as(pdf) replace

twoway (tsline sum_delta_c_all  if max_bank == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end)   ///
	(tsline  sum_c_all if max_bank == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end, yaxis(2))  , ///
	ytitle("",  height(6)) xtitle("", size(medsmall)) ///
	subtitle("%", position(11) size(small)) ///
	title(" ", pos(11) size(medlarge)) graphregion(color(white) style(none)) ///
	xlabel($charts_start(4)$charts_end, angle(45) labsize(medsmall)) ///
	legend(order(1 "{&Sigma}{subscript:i}{&delta}{subscript:i}c{subscript:i} " 2 "{&Sigma}{subscript:i}c{subscript:i}") r(1) region(style(none)) size(small))
graph export ../output/NVI_components_aggloss.pdf, as(pdf) replace

*Plots of index with range of different insurance, mutual assets %out
//twoway (tsline l_R_20_1 index_w_insurance index_w_insurance40 index_w_insurance50 index_w_insurance60 index_w_insurance70 index_w_insurance80 index_w_insurance90 index_w_insurance100 if max_bank == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end)  , ///
//	xtitle("", size(medsmall)) xlabel($charts_start(4)$charts_end, angle(45) labsize(medsmall)) ///
//	ytitle("",  height(6)) ///
//	title("Network Vulnerability Index", pos(11) size(medlarge)) 
//graph export ../output/index_insurance_grades.pdf, as(pdf) replace

//twoway (tsline l_R_20_1 index_w_mutual40 index_w_mutual50 index_w_mutual60 index_w_mutual70 index_w_mutual80 index_w_mutual90 index_w_mutual100 if max_bank == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end)  , ///
//	xtitle("", size(medsmall)) xlabel($charts_start(4)$charts_end, angle(45) labsize(medsmall)) ///
//	ytitle("",  height(6)) ///
//	title("Network Vulnerability Index", pos(11) size(medlarge)) 
//graph export ../output/index_mutual_grades.pdf, as(pdf) replace

//twoway (tsline l_R_20_1 index_w_dealers40 index_w_dealers50 index_w_dealers60 index_w_dealers70 index_w_dealers80 index_w_dealers90 index_w_dealers100 if max_bank == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end)  , ///
//	xtitle("", size(medsmall)) xlabel($charts_start(4)$charts_end, angle(45) labsize(medsmall)) ///
//	ytitle("",  height(6)) ///
//	title("Network Vulnerability Index", pos(11) size(medlarge)) 
//graph export ../output/index_dealers_grades.pdf, as(pdf) replace

*Estimated Outside Assets, from Flow of Funds
//tsline sum_c insurance_out mutual_out depository_out if max_bank ==1 & qt_dt >= $charts_start & qt_dt <= $charts_end, title("Outside Assets", size(med)) ytitle(" ") xtitle(" ") graphregion(color(white) style(none)) legend(order(1 "BHC Sample" 2 "Insurance" 3 "Mutual Funds" 4 "Banks, FFunds") size(small) r(3))
//graph export ../output/assets_out.pdf, as(pdf) replace

*Total Assets, from Flow of Funds
tsline sum_assets total_insurance total_mutual total_depository bhc_ffunds total_dealer if max_bank ==1 & qt_dt >= $charts_start & qt_dt <= $charts_end, title("All Assets", size(med)) ytitle(" ") xtitle(" ") graphregion(color(white) style(none)) legend(order(1 "BHC Sample" 2 "Insurance" 3 "Mutual Funds" 4 "Private Depository, FFunds" 5 "W/ Holding Companies" 6 "Broker Dealers") size(small) r(3))
graph export ../output/assets_all.pdf, as(pdf) replace

/*

preserve
	keep if keep_obs == 1
	keep qt_dt delta mkmv_id name tkr
	rename delta delta_v9
	save ../temp/temp.dta, replace

	use ../output/Contagion_Data_v8.dta, clear
	rename prob_default_phys delta_v8
	keep if !missing(delta_v8)
	keep qt_dt delta_v8 mkmv_id
	merge 1:1 qt_dt mkmv_id using ../temp/temp.dta

	//levelsof nameshort, local(big_banks)
	sort qt_dt
	format qt_dt  %tqCCYY
	foreach bank in C JPM BAC WFC {
		line delta_v8 delta_v9 qt_dt if tkr == "`bank'" & qt_dt >= $charts_start & qt_dt <= $charts_end, title("`bank'", size(small)) /// 
		ytitle(" ") legend(order(1 "EDF8" 2 "EDF9") size(small) r(1)) ///
		xtitle(" ", size(medsmall)) xlabel($charts_start(4)$charts_end, angle(45) labsize(small)) ylabel(,labsize(small)) ///
		graphregion(color(white) style(none)) name(`bank', replace)
		local file_label = subinstr("`bank'", " ", "_", 30)
		graph export ../output/KMV_`file_label'.pdf, as(pdf) replace
}
grc1leg JPM C BAC WFC, r(2) legendfrom(JPM) graphregion(color(white) style(none))
graph export ../output/KMV8v9.pdf, as(pdf) replace
restore
*/


***************************************************************************************************
*
*					Some Robustness Exercises
*
***************************************************************************************************

twoway (scatter beta_max_1 qt_dt if max_bank == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end, mlabgap(*3) mlabsize(vsmall) mlabel(tkr) mlabpos(11) mlabangle(rvertical)) || ///
	(scatter beta_max_fullsample_1 qt_dt if max_bank_fullsample == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end, mlabgap(*5) mlabsize(vsmall) mlabel(tkr) mlabpos(6) mlabangle(rvertical)), ///
	ytitle("",  height(6)) xtitle("", size(medsmall)) ylabel(.6 .7 .8 .9 1) ///
	title(" ", pos(11) size(medlarge)) graphregion(color(white) style(none)) /// 
	xlabel($charts_start(4)$charts_end, angle(45) labsize(medsmall)) ///
	legend(order(1 "{&beta}{superscript:+}: Top 20 BHCs" 2 "{&beta}{superscript:+}: Full Sample") r(1) region(style(none)) size(small)) ///
	subtitle("Inside Liabilities/Total Liabilities", position(11) size(small))
//graph play ytitle_offset
graph export ../output/scatter_beta_full.pdf, as(pdf) replace

//	bys qt_dt: egen beta_max_1_perc`perc' = max(liab_in_frac_temp) if keep_obs == 1
//	bys qt_dt: egen beta_max_1_perc`perc'_actual
//	gen liab_in_frac_temp`perc' =  (liab_in + (`perc'/100)*liab_in_unc)/BHCK2948
twoway (scatter beta_max_1_perc50_actual qt_dt if liab_in_frac_temp50 == beta_max_1_perc50_actual & qt_dt >= $charts_start & qt_dt <= $charts_end, mlabgap(*3) mlabsize(tiny) mlabel(tkr) mlabpos(11) mlabangle(rvertical)) || ///
	(scatter beta_max_1_perc30_actual qt_dt if liab_in_frac_temp30 == beta_max_1_perc30_actual & qt_dt >= $charts_start & qt_dt <= $charts_end, mlabgap(*3) mlabsize(tiny) mlabel(tkr) mlabpos(11) mlabangle(rvertical)) || ///
	(scatter beta_max_1_perc70_actual qt_dt if liab_in_frac_temp70 == beta_max_1_perc70_actual & qt_dt >= $charts_start & qt_dt <= $charts_end, mlabgap(*3) mlabsize(tiny) mlabel(tkr) mlabpos(11) mlabangle(rvertical)), ///
	ytitle("",  height(6)) xtitle("", size(medsmall)) ylabel(.6 .7 .8 .9 1) ///
	title(" ", pos(11) size(medlarge)) graphregion(color(white) style(none)) /// 
	xlabel($charts_start(4)$charts_end, angle(45) labsize(medsmall)) ///
	legend(order(1 "{&beta}{superscript:+}: 50% In" 2 "{&beta}{superscript:+}: 30% In" 3 "{&beta}{superscript:+}: 70% In") r(1) region(style(none)) size(small)) ///
	subtitle("Inside Liabilities/Total Liabilities", position(11) size(small))
//graph play ytitle_offset
graph export ../output/scatter_beta_perc_unc.pdf, as(pdf) replace

twoway (scatter beta_max_1 qt_dt if max_bank == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end, mlabgap(*3) mlabsize(vsmall) mlabel(tkr) mlabpos(11) mlabangle(rvertical)) || ///
	(scatter beta_max_fullsample_1 qt_dt if max_bank_fullsample == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end, mlabgap(*5) mlabsize(vsmall) mlabel(tkr) mlabpos(6) mlabangle(rvertical)), ///
	ytitle("",  height(6)) xtitle("", size(medsmall)) ylabel(.6 .7 .8 .9 1) ///
	title(" ", pos(11) size(medlarge)) graphregion(color(white) style(none)) /// 
	xlabel($charts_start(4)$charts_end, angle(45) labsize(medsmall)) ///
	legend(order(1 "{&beta}{superscript:+}: Top 20 BHCs" 2 "{&beta}{superscript:+}: Full Sample") r(1) region(style(none)) size(small)) ///
	subtitle("Inside Liabilities/Total Liabilities", position(11) size(small))
//graph play ytitle_offset
graph export ../output/scatter_beta_full.pdf, as(pdf) replace

*For fixed default prob
twoway (tsline nvi_benchmark if max_bank == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end) || ///
	(tsline nvi_benchmark_deltafix  if max_bank == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end), ///
	xtitle(" ", size(medsmall)) xlabel($charts_start(4)$charts_end, angle(45) labsize(small)) ylabel(,labsize(small)) ///
	ytitle(" ",  height(6)) subtitle("%", position(11) size(small)) ylabel(0 5 10 15 20 25 30,labsize(small)) ///
	legend(order(1 "Benchmark NVI" 2 "with Default Prob Fixed at 6%") r(2) region(style(none)) size(small)) ///
	title(" ", pos(11) size(medlarge)) graphregion(color(white) style(none)) name(nvi_deltafix, replace)
//graph play ytitle_offset
graph export ../output/nvi_deltafixed.pdf, as(pdf) replace


*With Fixed Default Prob

*For Sample inclusions
twoway (tsline nvi_benchmark if max_bank == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end) || ///
	(tsline index_y9c_only  if max_bank == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end) || ///
	(tsline index_w_insurance if max_bank == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end) || ///
	(tsline index_w_dealer if max_bank == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end) || /// 
	(tsline index_w_other if max_bank == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end) || ///
	(tsline index_w_reit if max_bank == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end),  ///
	xtitle(" ", size(medsmall)) xlabel($charts_start(4)$charts_end, angle(45) labsize(small)) ylabel(,labsize(small)) ///
	ytitle(" ",  height(6)) subtitle("%", position(11) size(small)) ///
	legend(order(1 "Benchmark NVI" 2 "BHC Only" 3 "BHC, Insurance Companies" 4 "BHC, Broker Dealers" 5 "BHC, Other" 6 "BHC, REITs") r(4) region(style(none)) size(small)) ///
	title(" ", pos(11) size(medlarge)) graphregion(color(white) style(none)) name(nvi_betanum, replace)
//graph play ytitle_offset
graph export ../output/nvi_sample_inclusions.pdf, as(pdf) replace


//	gen index_w_insurance_w_dealers`perc_label' = 100*(sum_delta_c + `perc'*(total_insurance)*delta_insurance + `perc'*(total_dealers_compustat)*delta_dealers)/((1-(1+20/100)*beta_max_1)*(sum_c+`perc'*total_dealers_compustat+`perc'*total_insurance))

//global gammas 1 5 10 15 20 

//foreach var of varlist index_w_ff_gam*{
//	replace `var' = 50 if `var' >= 50
//}
/*
twoway (tsline nvi_benchmark index_w_insurance_w_dealer10 index_w_insurance_w_dealer30 index_w_insurance_w_dealer50 index_w_insurance_w_dealer70 index_w_insurance_w_dealer90  if max_bank == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end),  ///
	xtitle(" ", size(medsmall)) xlabel($charts_start(4)$charts_end, angle(45) labsize(small)) ylabel(,labsize(small)) ///
	ytitle(" ",  height(6)) subtitle("%", position(11) size(small)) ///
	legend(order(1 "Default" 2 "10% Out" 3 "30% Out" 4 "50% Out" 5 "70% Out" 6 "90% Out") r(2) region(style(none)) size(small)) ///
	title(" ", pos(11) size(medlarge)) graphregion(color(white) style(none))
//graph play ytitle_offset
graph export ../output/robustness_percffunds.pdf, as(pdf) replace
*/
replace nvi_benchmark8 = . if nvi_benchmark8 == 0
twoway (tsline nvi_benchmark nvi_benchmark8 if max_bank == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end),  ///
	xtitle(" ", size(medsmall)) xlabel($charts_start(4)$charts_end, angle(45) labsize(small)) ylabel(,labsize(small)) ///
	ytitle(" ",  height(6)) subtitle("%", position(11) size(small)) ///
	legend(order(1 "Benchmark NVI (Moody's Expected Default Frequency, Version 9)" 2 "NVI using Moody's Expected Default Frequency, Version 8") r(2) region(style(none)) size(small)) ///
	title(" ", pos(11) size(medlarge)) graphregion(color(white) style(none))
//graph play ytitle_offset
graph export ../output/robustness_edf8vs9.pdf, as(pdf) replace

foreach var of varlist avg_delta8_assets delta_dealer8 delta_insurance8 delta_reit8 delta_other8{
	replace `var' = . if `var' == 0
}

foreach var of varlist delta_*8 avg_delta8_assets avg_delta_assets_full delta*full{
	replace `var' = . if `var' == 0
}
line avg_delta_assets avg_delta8_assets qt_dt if max_bank == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end, ///
	xtitle(" ", size(medsmall)) xlabel($charts_start(8)$charts_end, angle(45) labsize(small)) ylabel(,labsize(small)) ///
	ytitle(" ",  height(6)) subtitle("% per year", position(11) size(small)) ///
	legend(order(1 "Moody's  Expected Default Frequency, Version 9 (Benchmark)" 2 "Moody's  Expected Default Frequency, Version 8") r(2) region(style(none)) size(vsmall)) ///
	title("BHCs", size(medlarge)) graphregion(color(white) style(none)) name(y9c, replace)

line delta_dealer delta_dealer8 qt_dt if max_bank == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end, ///
	xtitle(" ", size(medsmall)) xlabel($charts_start(8)$charts_end, angle(45) labsize(small)) ylabel(,labsize(small)) ///
	ytitle(" ",  height(6)) subtitle("% per year", position(11) size(small)) ///
	legend(order(1 "Dealers, EDF9" 2 "Dealers, EDF8") r(2) region(style(none)) size(vsmall)) ///
	title("Broker Dealers", size(medlarge)) graphregion(color(white) style(none)) name(dealers, replace)
	
line delta_insurance delta_insurance8 qt_dt if max_bank == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end, ///
	xtitle(" ", size(medsmall)) xlabel($charts_start(8)$charts_end, angle(45) labsize(small)) ylabel(,labsize(small)) ///
	ytitle(" ",  height(6)) subtitle("% per year", position(11) size(small)) ///
	legend(order(1 "Insurance, EDF9" 2 "Insurance, EDF8") r(2) region(style(none)) size(vsmall)) ///
	title("Insurance Companies", size(medlarge)) graphregion(color(white) style(none)) name(insurance, replace)
	
line delta_reit delta_reit8 qt_dt if max_bank == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end, ///
	xtitle(" ", size(medsmall)) xlabel($charts_start(8)$charts_end, angle(45) labsize(small)) ylabel(,labsize(small)) ///
	ytitle(" ",  height(6)) subtitle("% per year", position(11) size(small)) ///
	legend(order(1 "REITs, EDF9" 2 "REITs, EDF8") r(2) region(style(none)) size(vsmall)) ///
	title("REITs", size(medlarge)) graphregion(color(white) style(none)) name(reit, replace)
	
line delta_other delta_other8 qt_dt if max_bank == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end, ///
	xtitle(" ", size(medsmall)) xlabel($charts_start(8)$charts_end, angle(45) labsize(small)) ylabel(,labsize(small)) ///
	ytitle(" ",  height(6)) subtitle("% per year", position(11) size(small)) ///
	legend(order(1 "EDF9" 2 "EDF8") r(1) region(style(none)) size(vsmall)) ///
	title("Other", size(medlarge)) graphregion(color(white) style(none)) name(other, replace)

grc1leg y9c dealers insurance reit other, title(" ", pos(11) size(medlarge)) graphregion(color(white) style(none)) legendfrom(y9c)
graph export ../output/edf8vs9_probs.pdf, as(pdf) replace

line avg_delta_assets avg_delta_assets_full qt_dt if max_bank_fullsample == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end, ///
	xtitle(" ", size(medsmall)) xlabel($charts_start(8)$charts_end, angle(45) labsize(small)) ylabel(,labsize(small)) ///
	ytitle(" ",  height(6)) subtitle("% per year", position(11) size(small)) ///
	legend(order(1 "Benchmark" 2 "Balanced Panel") r(1) region(style(none)) size(small)) ///
	title("BHCs", size(medlarge)) graphregion(color(white) style(none)) name(y9c, replace)

line delta_dealer delta_dealer_full qt_dt if max_bank_fullsample == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end, ///
	xtitle(" ", size(medsmall)) xlabel($charts_start(8)$charts_end, angle(45) labsize(small)) ylabel(,labsize(small)) ///
	ytitle(" ",  height(6)) subtitle("% per year", position(11) size(small)) ///
	legend(order(1 "Dealers, Benchmark" 2 "Dealers, Balanced Panel") r(1) region(style(none)) size(small)) ///
	title("Broker Dealers", size(medlarge)) graphregion(color(white) style(none)) name(dealers, replace)
	
line delta_insurance delta_insurance_full qt_dt if max_bank_fullsample == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end, ///
	xtitle(" ", size(medsmall)) xlabel($charts_start(8)$charts_end, angle(45) labsize(small)) ylabel(,labsize(small)) ///
	ytitle(" ",  height(6)) subtitle("% per year", position(11) size(small)) ///
	legend(order(1 "Insurance, Benchmark" 2 "Insurance, Balanced Panel") r(1) region(style(none)) size(small)) ///
	title("Insurance Companies", size(medlarge)) graphregion(color(white) style(none)) name(insurance, replace)
	
line delta_reit delta_reit_full qt_dt if max_bank_fullsample == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end, ///
	xtitle(" ", size(medsmall)) xlabel($charts_start(8)$charts_end, angle(45) labsize(small)) ylabel(,labsize(small)) ///
	ytitle(" ",  height(6)) subtitle("% per year", position(11) size(small)) ///
	legend(order(1 "Insurance, Benchmark" 2 "Insurance, Balanced Panel") r(1) region(style(none)) size(small)) ///
	title("REITs", size(medlarge)) graphregion(color(white) style(none)) name(reit, replace)
	
line delta_other delta_other_full qt_dt if max_bank_fullsample == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end, ///
	xtitle(" ", size(medsmall)) xlabel($charts_start(8)$charts_end, angle(45) labsize(small)) ylabel(,labsize(small)) ///
	ytitle(" ",  height(6)) subtitle("% per year", position(11) size(small)) ///
	legend(order(1 "Insurance, Benchmark" 2 "Insurance, Balanced Panel") r(1) region(style(none)) size(small)) ///
	title("Other", size(medlarge)) graphregion(color(white) style(none)) name(other, replace)

grc1leg y9c dealers insurance reit other, title(" ", pos(11) size(medlarge)) graphregion(color(white) style(none))
graph export ../output/full_vs_benchmark_probs.pdf, as(pdf) replace

twoway (tsline nvi_benchmark if max_bank == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end) || ///
	(tsline nvi_benchmark_full  if max_bank_fullsample == 1  & qt_dt >= $charts_start & qt_dt <= $charts_end), ///
	xtitle(" ", size(medsmall)) xlabel($charts_start(4)$charts_end, angle(45) labsize(small)) ylabel(,labsize(small)) ///
	ytitle(" ",  height(6)) subtitle("%", position(11) size(small)) ///
	legend(order(1 "Benchmark NVI" 2 "Balanced Panel Only") r(3) region(style(none)) size(small)) ///
	title(" ", pos(11) size(medlarge) color(black)) graphregion(color(white) style(none)) name(nvi_betas_selection, replace)
graph export ../output/robustness_fullsample.pdf, as(pdf) replace

*Max Connectivities using FR-Y15
twoway (scatter beta_y15 qt_dt if tkr == "STT" & qt_dt >= yq(2013, 1) & qt_dt <= $charts_end, mcolor(red) mlabgap(*5) mlabsize(tiny) mlabel(tkr) mlabpos(11) mlabangle(rvertical) mlabcolor(black)) || ///
	(scatter beta_y15 qt_dt if tkr == "BK" & qt_dt >= yq(2013, 1) & qt_dt <= $charts_end, mcolor(red) mlabgap(*5) mlabsize(tiny) mlabel(tkr) mlabpos(11) mlabangle(rvertical) mlabcolor(black)) || ///
	(scatter beta_y15 qt_dt if tkr == "JPM" & qt_dt >= yq(2013, 1) & qt_dt <= $charts_end, mcolor(red) mlabgap(*5) mlabsize(tiny) mlabel(tkr) mlabpos(11) mlabangle(rvertical) mlabcolor(black)) || ///
	(scatter beta_y15 qt_dt if tkr == "GS" & qt_dt >= yq(2013, 1) & qt_dt <= $charts_end, mcolor(red) mlabgap(*5) mlabsize(tiny) mlabel(tkr) mlabpos(11) mlabangle(rvertical) mlabcolor(black)) || ///
	(scatter beta qt_dt if tkr == "JPM" & qt_dt >= yq(2013, 1) & qt_dt <= $charts_end, mcolor(green) mlabgap(*5) mlabsize(tiny) mlabel(tkr) mlabpos(11) mlabangle(rvertical) mlabcolor(black)) || ///
	(scatter beta qt_dt if tkr == "GS" & qt_dt >= yq(2013, 1) & qt_dt <= $charts_end, mcolor(green) mlabgap(*5) mlabsize(tiny) mlabel(tkr) mlabpos(11) mlabangle(rvertical) mlabcolor(black)) || ///
	(scatter beta qt_dt if tkr == "STT" & qt_dt >= yq(2013, 1) & qt_dt <= $charts_end, mcolor(green) mlabgap(*5) mlabsize(tiny) mlabel(tkr) mlabpos(11) mlabangle(rvertical) mlabcolor(black)) || ///
	(scatter beta qt_dt if tkr == "BK" & qt_dt >= yq(2013, 1) & qt_dt <= $charts_end, mcolor(green) mlabgap(*5) mlabsize(tiny) mlabel(tkr) mlabpos(11) mlabangle(rvertical) mlabcolor(black)), ///	
	xtitle(" ", size(medsmall)) xlabel(, angle(45) labsize(small)) ylabel(0 .25 .5 .75 1,labsize(small)) ///
	ytitle(" ",  height(6)) subtitle("Inside Liabilities/Total Liabilities", position(11) size(small)) ///
	legend(order(1 "Y-15" 5 "Y-9C") r(4) region(style(none)) size(small)) ///
	title(" ", pos(11) size(medlarge) color(black)) graphregion(color(white) style(none)) name(scatter_betas_y15, replace)
graph export ../output/scatter_connect_y15.pdf, as(pdf) replace

twoway (scatter beta_y15_compat qt_dt if tkr == "STT" & qt_dt >= yq(2013, 1) & qt_dt <= $charts_end, mcolor(red) mlabgap(*5) mlabsize(tiny) mlabel(tkr) mlabpos(11) mlabangle(rvertical) mlabcolor(black)) || ///
	(scatter beta_y15_compat qt_dt if tkr == "BK" & qt_dt >= yq(2013, 1) & qt_dt <= $charts_end, mcolor(red) mlabgap(*5) mlabsize(tiny) mlabel(tkr) mlabpos(11) mlabangle(rvertical) mlabcolor(black)) || ///
	(scatter beta_y15_compat qt_dt if tkr == "JPM" & qt_dt >= yq(2013, 1) & qt_dt <= $charts_end, mcolor(red) mlabgap(*5) mlabsize(tiny) mlabel(tkr) mlabpos(11) mlabangle(rvertical) mlabcolor(black)) || ///
	(scatter beta_y15_compat qt_dt if tkr == "GS" & qt_dt >= yq(2013, 1) & qt_dt <= $charts_end, mcolor(red) mlabgap(*5) mlabsize(tiny) mlabel(tkr) mlabpos(11) mlabangle(rvertical) mlabcolor(black)) || ///
	(scatter beta qt_dt if tkr == "JPM" & qt_dt >= yq(2013, 1) & qt_dt <= $charts_end, mcolor(green) mlabgap(*5) mlabsize(tiny) mlabel(tkr) mlabpos(11) mlabangle(rvertical) mlabcolor(black)) || ///
	(scatter beta qt_dt if tkr == "GS" & qt_dt >= yq(2013, 1) & qt_dt <= $charts_end, mcolor(green)  mlabgap(*5) mlabsize(tiny) mlabel(tkr) mlabpos(11) mlabangle(rvertical) mlabcolor(black)) || ///
	(scatter beta qt_dt if tkr == "STT" & qt_dt >= yq(2013, 1) & qt_dt <= $charts_end, mcolor(green)  mlabgap(*5) mlabsize(tiny) mlabel(tkr) mlabpos(11) mlabangle(rvertical) mlabcolor(black)) || ///
	(scatter beta qt_dt if tkr == "BK" & qt_dt >= yq(2013, 1) & qt_dt <= $charts_end, mcolor(green)  mlabgap(*5) mlabsize(tiny) mlabel(tkr) mlabpos(11) mlabangle(rvertical) mlabcolor(black)), ///	
	xtitle(" ", size(medsmall)) xlabel(, angle(45) labsize(small)) ylabel(0 .25 .5 .75 1,labsize(small)) ///
	ytitle(" ",  height(6)) subtitle("Inside Liabilities/Total Liabilities", position(11) size(small)) ///
	legend(order(1 "Y-15" 5 "Y-9C" ) r(4) region(style(none)) size(small)) ///
	title(" ", pos(11) size(medlarge) color(black)) graphregion(color(white) style(none)) name(scatter_betas_y15, replace)
graph export ../output/scatter_connect_y15_compat.pdf, as(pdf) replace
graph export ../output/scatter_connect_y15_compat.png, as(png) replace

twoway (scatter beta_max_y15_1 qt_dt if  max_bank_y15 == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end, mcolor(blue) mlabgap(*5) mlabsize(vsmall) mlabel(tkr) mlabpos(11) mlabangle(rvertical)) || ///
	(scatter beta_max_y15_2 qt_dt if  max_bank_y15 == 2 & qt_dt >= $charts_start & qt_dt <= $charts_end, mcolor(red) mlabgap(*5) mlabsize(vsmall) mlabel(tkr) mlabpos(11) mlabangle(rvertical)) || ///
	(scatter beta_max_y15_3 qt_dt if  max_bank_y15 == 3 & qt_dt >= $charts_start & qt_dt <= $charts_end, mcolor(green) mlabgap(*5) mlabsize(vsmall) mlabel(tkr) mlabpos(6) mlabangle(rvertical)) || ///
	(scatter beta_max_y15_4 qt_dt if  max_bank_y15 == 4 & qt_dt >= $charts_start & qt_dt <= $charts_end, mcolor(maroon) mlabgap(*5) mlabsize(vsmall) mlabel(tkr) mlabpos(6) mlabangle(rvertical)), ///
	ytitle("",  height(6)) xtitle("", size(medsmall)) ylabel(.4 .5 .6 .7 .8 .9) ///
	title(" ", pos(11) size(medlarge)) graphregion(color(white) style(none)) /// 
	xlabel($charts_start(4)$charts_end, angle(45) labsize(medsmall)) ///
	legend(order(1 "Highest {&beta}" 2 "Second Highest {&beta}" 3 "Third Highest {&beta}" 4 "Fourth Highest {&beta}") r(1) region(style(none)) size(small)) ///
	subtitle("Inside Liabilities/Total Liabilities", position(11) size(small)) name(beta_first, replace)
graph export ../output/betas_y15_test.pdf, as(pdf) replace



*New dealer nodes stuff
replace delta = delta * 100
replace delta_dealer = delta_dealer * 100

twoway (line delta qt_dt if tkr == "BRO10" & qt_dt >= $charts_start & qt_dt <= $charts_end) || ///
	(line delta qt_dt if tkr == "BRO25" & qt_dt >= $charts_start & qt_dt <= $charts_end) || ///
	(line delta_dealer qt_dt if qt_dt >= $charts_start & qt_dt <= $charts_end ), ///
	xtitle(" ", size(medsmall)) xlabel($charts_start(4)$charts_end, angle(45) labsize(small)) ylabel(,labsize(small)) ///
	subtitle("% per year", position(11) size(small)) ///
	ytitle(" ",  height(6)) ///
	legend(order(1 "Top 10 Dealers" 2 "Top 11-25 Dealers" 3 "All Dealers") r(2) region(style(none)) size(small)) ///
	title(" ", pos(11) size(medlarge)) graphregion(color(white) style(none))
graph export ../output/delta_externals_dealernodes.pdf, as(pdf) replace
//graph export ../output/delta_externals.png, as(png) replace


twoway (tsline nvi_benchmark_test if max_bank_dealertest == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end) || ///
	(tsline nvi_benchmark if max_bank == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end), ///
	ytitle("",  height(6)) xtitle("", size(medsmall)) ///
	subtitle("%", position(11) size(small)) ///
	title(" ", pos(11) size(medlarge)) graphregion(color(white) style(none)) /// 
	xlabel($charts_start(4)$charts_end, angle(45) labsize(small)) ///
	legend(order(1 "With Dealer Nodes" 2 "Original") r(1) region(style(none)) size(small)) ///
	ylabel(0 5 10 15 20 25 30, labsize(small))
graph export ../output/NVI_benchmark_dealernodes.pdf, as(pdf) replace
*graph export ../output/NVI_benchmark.png, as(png) replace

twoway (scatter beta_max_1 qt_dt if max_bank == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end, mcolor(blue) mlabgap(*5) mlabsize(vsmall) mlabel(tkr) mlabpos(11) mlabangle(rvertical)) || ///
	(scatter beta_max_1_deal qt_dt if max_bank_dealertest == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end, mcolor(red) mlabgap(*5) mlabsize(vsmall) mlabel(tkr) mlabpos(11) mlabangle(rvertical)), ///
	ytitle("",  height(6)) xtitle("", size(medsmall)) ylabel(.4 .5 .6 .7 .8 .9) ///
	title(" ", pos(11) size(medlarge)) graphregion(color(white) style(none)) /// 
	xlabel($charts_start(4)$charts_end, angle(45) labsize(medsmall)) ///
	legend(order(1 "Original {&beta}{superscript:+}" 2 "With Dealer Nodes") r(1) region(style(none)) size(small)) ///
	subtitle("Inside Liabilities/Total Liabilities", position(11) size(small)) name(beta_dealer_test, replace)
graph export ../output/connectivities_dealernodes.pdf, as(pdf) replace



*Testing how our call report matching went. Plotting different beta^+ for entities that didn't match that first time
levelsof tkr if matchedon == 2, local(tickers_nocall)
local num: word count `tickers_nocall'
local i=1
forv j=1(1)`=ceil(`num'/16)'{
	local to_combine_`j'
}
foreach tkr in `tickers_nocall'{
	local to_combine_`=ceil(`i'/16)' `to_combine_`=ceil(`i'/16)'' `tkr'
	local first_
	
	twoway (scatter beta qt_dt if tkr == "`tkr'" & qt_dt >= $charts_start & qt_dt <= $charts_end, mcolor(blue)) || ///
		(scatter liabs_in_nocall_frac qt_dt if tkr == "`tkr'" & qt_dt >= $charts_start & qt_dt <= $charts_end, mcolor(red)), ///
		ytitle("",  height(6)) xtitle("", size(medsmall)) ylabel(0 .25 .5 .75 1,angle(0) labsize(small)) ///
		title("`tkr'", size(medlarge)) graphregion(color(white) style(none)) /// 
		xlabel($charts_start(16)$charts_end, angle(45) labsize(medsmall)) ///
		legend(order(1 "With Hierarchical Match" 2 "With RSSD-Only Match") r(1) region(style(none)) size(small)) ///
		subtitle("Inside Liabilities/Total Liabilities", position(11) size(small)) name(`tkr', replace)

	local i = `i' + 1
		}
	graph close
forv j=1(1)`=ceil(`num'/16)'{
	if `j' == `=ceil(`num'/16)'{
		local rows 3
	}
	else{
		local rows 4
	}
	grc1leg `to_combine_`j'', graphregion(color(white) style(none)) rows(`rows') plotregion(margin(none))
	graph export ../output/call_match_compare`j'.pdf, replace
}








foreach var of varlist assets_kmv_* *_orig{
	replace `var' = `var'/1000000
}
tsline assets_kmv_dealer assets_kmv_reit assets_kmv_insurance assets_kmv_other if max_bank == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end, ///
	xtitle(" ", size(medsmall)) xlabel($charts_start(4)$charts_end, angle(45) labsize(small)) ylabel(,labsize(small)) ///
	ytitle(" ",  height(6)) subtitle("Sector Sample Assets", position(11) size(small)) ///
	legend(order(1 "Broker Dealers" 2 "REITs" 3 "Insurance" 4 "Other") r(2) region(style(none)) size(small)) ///
	title(" ", pos(11) size(medlarge)) graphregion(color(white) style(none))
graph export ../output/kmv_assets.pdf, as(pdf) replace
graph export ../output/kmv_assets.png, as(png) replace

tsline total_dealer_orig total_reit_orig total_insurance_orig total_other_orig if max_bank == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end, ///
	xtitle(" ", size(medsmall)) xlabel($charts_start(4)$charts_end, angle(45) labsize(small)) ylabel(,labsize(small)) ///
	ytitle(" ",  height(6)) subtitle("Total Flow of Funds Assets", position(11) size(small)) ///
	legend(order(1 "Broker Dealers" 2 "REITs" 3 "Insurance" 4 "Other") r(2) region(style(none)) size(small)) ///
	title(" ", pos(11) size(medlarge)) graphregion(color(white) style(none))
graph export ../output/ffunds_assets.pdf, as(pdf) replace
graph export ../output/ffunds_assets.png, as(png) replace

tsline coverage_insurance coverage_bhc coverage_reit coverage_other if max_bank_bhc == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end, ///
	xtitle(" ", size(medsmall)) xlabel($charts_start(4)$charts_end, angle(45) labsize(small)) ylabel(,labsize(small)) ///
	ytitle(" ",  height(6)) subtitle("Sector Sample Assets / Flow of Funds Total Assets", position(11) size(small)) ///
	legend(order(1 "Insurance" 2 "BHCs" 3 "REITs" 4 "Other") r(2) region(style(none)) size(small)) ///
	title(" ", pos(11) size(medlarge)) graphregion(color(white) style(none))
graph export ../output/coverage.pdf, as(pdf) replace
graph export ../output/coverage.png, as(png) replace



*Figure 4: Reasons for CDealer Coverage Drop During and After Financial Crisis

preserve
local axis_label = yq(2013, 1)
local line1 = yq(2008, 3)
local line2 = yq(2009, 1) 
replace BHCK2170 = BHCK2170/1000000000
replace total_dealer = total_dealer/1000000
twoway (line coverage_dealer qt_dt if qt_dt >= $charts_start & qt_dt <= $charts_end, yaxis(2)) || ///
	(line BHCK2170 qt_dt if tkr == "BAC" & qt_dt >= $charts_start & qt_dt <= $charts_end) || ///
	(line BHCK2170 qt_dt if tkr == "GS" & qt_dt >= $charts_start & qt_dt <= $charts_end) || ///
	(line BHCK2170 qt_dt if tkr == "MS" & qt_dt >= $charts_start & qt_dt <= $charts_end) || ///
	(line total_dealer qt_dt if max_bank == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end), ///
	xtitle(" ", size(medsmall)) xlabel($charts_start(4)$charts_end, angle(45) labsize(small)) ylabel(0 1 2 3 4 5,labsize(small)) ylabel(0 .25 .5 .75 1 ,axis(2) labsize(small)) ///
	ytitle(" ", axis(1) height(6)) ytitle(" ", axis(2) height(6)) subtitle("Broker-Dealer Coverage", position(11) size(small)) ///
	legend(order(1 "Broker-Dealer Coverage" 2 "BAC Assets" 3 "GS Assets" 4 "MS Assets" 5 "Dealer Assets, Outside FR-Y9C") r(3) region(style(none)) size(small)) ///
	title("", pos(11) size(medlarge)) graphregion(color(white) style(none)) ///
	text(5.38 `axis_label' "Total Assets (Y9C), Trillions", size(small) placement(east))
graph export ../output/dealer_coverage_details.pdf, as(pdf) replace
graph export ../output/dealer_coverage_details.png, as(png) replace
restore



twoway (scatter beta_max_1 qt_dt if max_bank == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end, mlabgap(*5) mlabsize(tiny) mlabel(tkr) mlabpos(11) mlabangle(rvertical)) || ///
	(scatter beta_max_all_1 qt_dt if beta_max_all_1<1 & max_bank_all == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end, mlabgap(*5) mlabsize(tiny) mlabel(tkr) mlabpos(6) mlabangle(rvertical)) || ///
	(scatter beta_max_all_nocall_1 qt_dt if beta_max_all_nocall_1<1 & max_bank_nocall == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end, mlabgap(*5) mlabsize(tiny) mlabel(tkr) mlabpos(6) mlabangle(rvertical)), ///
	xtitle(" ", size(medsmall)) xlabel($charts_start(4)$charts_end, angle(45) labsize(small)) ylabel(.5 .6 .7 .8 .9 1,labsize(small)) ///
	ytitle(" ",  height(6)) subtitle("Inside Liabilities/Total Liabilities", position(11) size(small)) ///
	legend(order(1 "Benchmark BHC Connectivity ({&beta}{superscript:+})" 2 "Highest BHC Connectivity, from Entire Sample ({&beta}{superscript:+})" 3 "Highest BHC Connectivity, from Entire Sample with only RSSD matching ({&beta}{superscript:+})" ) r(4) region(style(none)) size(small)) ///
	title("(b)", pos(11) size(medlarge) color(black)) graphregion(color(white) style(none)) name(scatter_betas_matching_call, replace)
graph export ../output/scatter_betanum_matching.pdf, as(pdf) replace

twoway (scatter beta qt_dt if max_bank_all == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end, mlabgap(*5) mlabsize(tiny) mlabel(tkr) mlabpos(11) mlabangle(rvertical)) || ///
	(scatter beta qt_dt if max_bank_all == 2 & qt_dt >= $charts_start & qt_dt <= $charts_end, mlabgap(*5) mlabsize(tiny) mlabel(tkr) mlabpos(6) mlabangle(rvertical)) || ///
	(scatter beta qt_dt if max_bank_all == 3 & qt_dt >= $charts_start & qt_dt <= $charts_end, mlabgap(*5) mlabsize(tiny) mlabel(tkr) mlabpos(6) mlabangle(rvertical)), ///
	xtitle(" ", size(medsmall)) xlabel($charts_start(4)$charts_end, angle(45) labsize(small)) ylabel(.5 .6 .7 .8 .9 1,labsize(small)) ///
	ytitle(" ",  height(6)) subtitle("Inside Liabilities/Total Liabilities", position(11) size(small)) ///
	legend(order(1 "Highest ({&beta}{superscript:+})" 2 "Second Highest ({&beta}{superscript:+})" 3 "Third Highest ({&beta}{superscript:+})" ) r(4) region(style(none)) size(small)) ///
	title("", pos(11) size(medlarge) color(black)) graphregion(color(white) style(none)) name(scatter_betanum_wholesample, replace)
graph export ../output/scatter_betanum_wholesample.pdf, as(pdf) replace

*Different Tiers for Max Beta


twoway (tsline nvi_benchmark if max_bank == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end) || ///
	(tsline index_w_ff_beta_all1 if max_bank_all == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end) || ///
	(tsline nvi_benchmark_all_nocall if max_bank_nocall == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end), ///
	xtitle(" ", size(medsmall)) xlabel($charts_start(4)$charts_end, angle(45) labsize(small)) ylabel(,labsize(small)) ///
	ytitle(" ",  height(6)) subtitle("%", position(11) size(small)) ylabel(,labsize(small)) ///
	legend(order(1 "Benchmark NVI" 2 "Highest BHC Connectivity, from Entire Sample" 3 "Highest BHC Connectivity, from Entire Sample with only RSSD Matching" ) r(4) region(style(none)) size(small)) ///
	title("(a)", pos(11) size(medlarge)) graphregion(color(white) style(none)) name(nvi_betanum_callmatch, replace)
//graph play ytitle_offset
graph export ../output/nvi_betanum_matchingcall.pdf, as(pdf) replace


graph combine nvi_betanum_callmatch scatter_betas_matching_call, title(" ", pos(11) size(medlarge)) graphregion(color(white) style(none))
graph export ../output/robustness_beta_selection_callmatch.pdf, as(pdf) replace












/*
*Playing around with default probs around bankruptcies
stop
clear
odbc load, exec("SELECT mkmv_id,annualized_edf_01yr_pct, risk_neutral_annualized_qdf_01yr_pct, dt_day,mkmvind_name_desc, total_assets_current_amt, kmv_country_id FROM vw_kmv_v9_dfs_v8map WHERE mkmv_id IN ('N13561', 'N05616', '524908', '882149', 'N11320', 'N23734')") ///
	connectionstring("DRIVER={SQL Server};SERVER=m1-wdmsql01.frb.gov;DATABASE=WDM;")
save ../temp/kmvs_select_bankruptcies, replace

use ../temp/CUSIP_MKMVID.dta, clear
keep mkmv_id company_name tkr
duplicates drop mkmv_id, force
merge 1:m mkmv_id using ../temp/kmvs_select_bankruptcies
gen date = date(dt_day, "YMD")
format date %tdCCYY
rename annualized_edf_01yr_pct edf01
keep if _merge == 3
drop _merge


keep if date >= dofq($charts_start) & date <= dofq($charts_end)
merge m:1 mkmv_id using ../temp/bankruptcies.dta
keep if _merge == 3
drop _merge
levelsof mkmv_id, clean local(firms)
sort mkmv_id date
local all	
foreach id in `firms'{
local all `all' probs_`id'
	preserve 
	keep if mkmv_id == "`id'"
	replace bankruptcy_datetime2 = td(01jan2017) if missing(bankruptcy_datetime2)
	replace resolution_datetime2 = td(01jan2017) if missing(resolution_datetime2)
	local bank = bankruptcy_datetime1[1]
	local bank_res = resolution_datetime1[1]
	local bank2 = bankruptcy_datetime2[1]
	local bank_res2 = resolution_datetime2[1]
	local comp_name = tkr[1]
	line edf01 date if date >= dofq($charts_start) & date <= dofq($charts_end), ///
	xtitle(" ", size(medsmall)) xlabel(,angle(45) labsize(small)) ylabel(,labsize(small)) ///
	ytitle(" ",  height(6)) subtitle("", position(11) size(small)) ///
	legend(r(1) region(style(none)) size(small)) ///
	title("`comp_name'", pos(11) size(medlarge)) graphregion(color(white) style(none)) name(probs_`id', replace) ///
	xline(`bank' `bank2', lcolor(green) lpattern(dash)) ///
	xline(`bank_res' `bank_res2', lcolor(red) lpattern(dash))
	//graph export ../output/bankruptcy_prob_`id'.pdf, as(pdf) replace
	restore
}
graph combine `all'
graph export ../output/some_bankruptcy_probs.pdf, as(pdf) replace

levelsof mkmv_id, clean local(y9c_sample)
global y9c_sample `y9c_sample'


use ../temp/CUSIP_MKMVID.dta, clear
keep mkmv_id company_name tkr
duplicates drop mkmv_id, force
merge 1:m mkmv_id using ../temp/high_edf_drops
keep if _merge == 3
collapse (firstnm) mkmvind_name_desc annualized_edf_01yr_pct date0 (mean) total_assets_current_amt, by(mkmv_id)
gsort mkmvind_name_desc -total_assets_current_amt
by mkmvind_name_desc: gen rank = _n
keep if rank <= 4

levelsof mkmv_id, local(ids)
local query
foreach id in `ids'{
	local query `query' '`id'', 
	}

local query = substr("`query'", 1, length("`query'")-1)
disp "`query'"
clear
odbc load, exec("SELECT mkmv_id,annualized_edf_01yr_pct, risk_neutral_annualized_qdf_01yr_pct, dt_day,mkmvind_name_desc, total_assets_current_amt, kmv_country_id FROM vw_kmv_v9_dfs_v8map WHERE mkmv_id IN (`query')") ///
	connectionstring("DRIVER={SQL Server};SERVER=m1-wdmsql01.frb.gov;DATABASE=WDM;")
save ../temp/kmvs_select_bankruptcies, replace


use ../temp/CUSIP_MKMVID.dta, clear
keep mkmv_id company_name tkr
duplicates drop mkmv_id, force
merge 1:m mkmv_id using ../temp/kmvs_select_bankruptcies
gen date = date(dt_day, "YMD")
format date %tdCCYY
rename annualized_edf_01yr_pct edf01
keep if _merge == 3
drop _merge


keep if date >= dofq($charts_start) & date <= dofq($charts_end)
merge m:1 mkmv_id using ../temp/bankruptcies.dta
keep if _merge == 3 | _merge == 1
drop _merge
sort mkmv_id date
levelsof mkmvind_name_desc, local(sectors)
disp "1"
foreach sector in `sectors'{
local name = subinstr("`sector'", " ", "_", 10)
local name = subinstr("`name'", "-", "", 10)
local name = subinstr("`name'", "/", "_", 10)
disp "`name'"
disp "2"
local all	
levelsof mkmv_id if mkmvind_name_desc == "`sector'", clean local(ids)
foreach id in `ids'{
local all `all' probs_`id'
	preserve 
	keep if mkmv_id == "`id'"
	replace bankruptcy_datetime2 = td(01jan2020) if missing(bankruptcy_datetime2)
	replace resolution_datetime2 = td(01jan2020) if missing(resolution_datetime2)
	replace resolution_datetime1 = td(01jan2020) if missing(resolution_datetime1)
	replace bankruptcy_datetime1 = td(01jan2020) if missing(bankruptcy_datetime1)
	local bank = bankruptcy_datetime1[1]
	local bank_res = resolution_datetime1[1]
	local bank2 = bankruptcy_datetime2[1]
	local bank_res2 = resolution_datetime2[1]
	local comp_name = company_name[1]
	line edf01 date if date >= dofq($charts_start) & date <= dofq($charts_end), ///
		xtitle(" ", size(medsmall)) xlabel(,angle(45) labsize(small)) ylabel(,labsize(small)) ///
		ytitle(" ",  height(6)) subtitle("", position(11) size(small)) ///
		legend(r(1) region(style(none)) size(small)) ///
		title("`comp_name'", pos(11) size(medlarge)) graphregion(color(white) style(none)) name(probs_`id', replace) ///
		xline(`bank', lcolor(green) lpattern(dash)) ///
		xline(`bank_res', lcolor(red) lpattern(dash))
		restore

}
graph combine `all'
disp "3"
graph export ../output/some_bankruptcy_probs_`name'.pdf, as(pdf) replace
}

use ../temp/CUSIP_MKMVID.dta, clear
keep mkmv_id company_name tkr
duplicates drop mkmv_id, force
merge 1:m mkmv_id using ../temp/DFS
keep if _merge == 3

keep if inlist(mkmvind_name_desc, "FINANCE COMPANIES", "FINANCE NEC", "INVESTMENT MANAGEMENT")
keep if kmv_country_id == "USA"
drop if missing(edf01)
drop if qt_dt < $charts_start
preserve

collapse (mean) total_assets_current_amt (lastnm) mkmvind_name_desc company_name, by(mkmv_id)
gsort mkmvind_name_desc -total_assets_current_amt
keep if !strpos(" $y9c_sample ", mkmv_id)
by mkmvind_name_desc: gen rank = _n
//keep if rank <= 4
rename company_name cname
rename mkmvind_name_desc sector_temp
keep mkmv_id cname sector_temp rank
save ../temp/big_mkmvid_not_in, replace
restore

preserve
keep if !strpos(" $y9c_sample ", mkmv_id)
format qt_dt %tq
drop if mkmv_id == "313400" | mkmv_id == "313586"
collapse (sum) total_assets_current_amt, by(mkmvind_name_desc qt_dt)
sort qt_dt
twoway (line total_assets_current_amt qt_dt if mkmvind_name_desc == "INVESTMENT MANAGEMENT") || ///
	(line total_assets_current_amt qt_dt if mkmvind_name_desc == "FINANCE NEC") || ///
	(line total_assets_current_amt qt_dt if mkmvind_name_desc == "FINANCE COMPANIES"), ///
	xtitle(" ", size(medsmall)) xlabel(,angle(45) labsize(small)) ylabel(,labsize(small)) ///
	ytitle(" ",  height(6)) subtitle("", position(11) size(small)) ///
	legend(order(1 "Investment Management" 2 "NEC" 3 "Finance Companies") r(1) region(style(none)) size(small)) ///
	title("`comp_name'", pos(11) size(medlarge)) graphregion(color(white) style(none))
graph export ../output/assets_other_KMV_cats.pdf, as(pdf) replace
restore


drop _merge
merge m:1 mkmv_id using ../temp/big_mkmvid_not_in
keep if _merge == 3
keep if rank <= 4

sort mkmv_id qt_dt
format qt_dt %tq
levelsof mkmvind_name_desc, local(sectors)
foreach sector in `sectors'{
local all
local name = subinstr("`sector'", " ", "_", 10)
local name = subinstr("`name'", "-", "", 10)
local name = subinstr("`name'", "/", "_", 10)
levelsof mkmv_id if mkmvind_name_desc == "`sector'", clean local(ids)
foreach id in `ids'{
local all `all' probs_`id'
	preserve 
	keep if mkmv_id == "`id'"
	local comp_name = company_name[1]
	line edf01 qt_dt if qt_dt >= $charts_start & qt_dt <= $charts_end, ///
		xtitle(" ", size(medsmall)) xlabel(,angle(45) labsize(small)) ylabel(,labsize(small)) ///
		ytitle(" ",  height(6)) subtitle("", position(11) size(small)) ///
		legend(r(1) region(style(none)) size(small)) ///
		title("`comp_name'", pos(11) size(medlarge)) graphregion(color(white) style(none)) name(probs_`id', replace)
		restore

}
graph combine `all'
//graph export ../output/large_firms_notin_`name'.pdf, as(pdf) replace

}
*/
