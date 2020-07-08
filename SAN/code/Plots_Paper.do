*********************************************************************************
*
*	Network Chart Production, Paper
*	Collin Jones 
*	Purpose: Output Tables and Charts for use in paper. Separate do-file for 
*		appendix (and other) charts 
*
*********************************************************************************

use ../output/Contagion_Data_select, clear

*Figure 1: A simple network
*		Written in Tikz

*Figure 2: A Simple Disconnected Network
*		Written in Tikz

*Figure 3: Total  Financial  Assets  for  each  Network  Subsector,  as  a  Percentage  of  Total Network Assets
local to_add bhc_ffunds total_insurance_orig total_dealer_orig total_reit_orig total_other_orig
local passed
local to_graph
local i : word count `to_add'
local legend
foreach var of varlist `to_add'{
	gen `var'_perc = 100*`var'/total_financial_ours `passed'
	local passed + `var'_perc
	disp "`passed'"
	local to_graph `var'_perc `to_graph'
	local name = strproper(subinstr(subinstr(subinstr("`var'", "_orig", "", 1), "total_", "", 1)), "_ffunds", "", 1)
	local name = subinstr("`name'", "Bhc", "BHC", 1)
	local name = subinstr("`name'", "Reit", "REIT", 1)
	local name = subinstr("`name'", "Dealer", "Broker-Dealer", 1)	
	local legend `"`legend' `i' "`name'""'
	local i = `i' - 1
}
twoway area `to_graph' qt_dt if max_bank == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end, ///
	xtitle(" ", size(medsmall)) xlabel($charts_start(4)$charts_end, angle(45) labsize(small)) ylabel(0 25 50 75 100,labsize(small)) ///
	ytitle(" ",  height(6)) subtitle("% of Total Network Assets", position(11) size(small)) ///
	legend(order(`legend') r(2) region(style(none)) size(small)) ///
	title(" ", pos(11) size(medlarge)) graphregion(color(white) style(none))
graph export ../output/ffunds_assets_area.pdf, as(pdf) replace
graph export ../output/ffunds_assets_area.png, as(png) replace

*Figure 4: Overall FFunds Coverage
tsline coverage_included coverage_possible if max_bank == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end, ///
	xtitle(" ", size(medsmall)) xlabel($charts_start(4)$charts_end, angle(45) labsize(small)) ylabel(0 .25 .5 .75 1,labsize(small)) ///
	ytitle(" ",  height(6)) subtitle("Sample Assets / Flow of Funds Total Assets", position(11) size(small)) ///
	legend(order(1 "Assets of Individual Firms Included" 2 "Assets of Approximated Nodes Included") r(3) region(style(none)) size(small)) ///
	title(" ", pos(11) size(medlarge)) graphregion(color(white) style(none))
graph export ../output/coverage_broad.png, as(png) replace
graph export ../output/coverage_broad.pdf, as(pdf) replace

*Figure 5: Network Vulnerability Index
twoway (tsline nvi_benchmark if max_bank == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end)  , ///
	ytitle("",  height(6)) xtitle("", size(medsmall)) ///
	subtitle("%", position(11) size(small)) ///
	title(" ", pos(11) size(medlarge)) graphregion(color(white) style(none)) /// 
	xlabel($charts_start(4)$charts_end, angle(45) labsize(small)) ///
	legend(order(1 "Amplification of Expected Losses Due to Network Spillover Effects") r(1) region(style(none)) size(small) on) ///
	ylabel(0 5 10 15 20 25 30 35 40, labsize(small))
graph export ../output/NVI_benchmark.pdf, as(pdf) replace
graph export ../output/NVI_benchmark.png, as(png) replace

*Figure 6: Multiplicative Components of NVI
local axis_label = yq(2015, 3)
twoway (tsline agg_loss if max_bank == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end) ///
	(tsline connectivity_comp if max_bank == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end, yaxis(2))  , ///
	ytitle("",  height(6)) xtitle("", size(medsmall)) ytitle("",  height(6) axis(2)) ///
	subtitle("% per year", position(11) size(small)) ///
	title(" ", pos(11) size(medlarge)) graphregion(color(white) style(none)) /// 
	xlabel($charts_start(4)$charts_end, angle(45) labsize(small)) ///
	legend(order(1 "Asset-Weighted Default Probability (left axis)" 2 "Connectivity Multiplier (right axis)") r(2) region(style(none)) size(small)) ///
	ylabel(0 5 10, axis(1) labsize(small)) ylabel(, axis(2) labsize(small))
graph export ../output/NVI_components.pdf, as(pdf) replace
graph export ../output/NVI_components.png, as(png) replace

*Figure 7: Maximum Liability Connectivity Among Large BHCs (and other high connectivities)
scatter beta_max_1 qt_dt if max_bank == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end, ///
	mcolor(blue) mlabgap(*5) mlabsize(small) mlabel(tkr) mlabpos(11) mlabangle(rvertical) ///
	ytitle("",  height(6)) xtitle("", size(medsmall)) ylabel(.5 .6 .7 .8) ///
	title(" ", pos(11) size(medlarge)) graphregion(color(white) style(none)) /// 
	xlabel($charts_start(4)$charts_end, angle(45) labsize(medsmall)) ///
	subtitle("Inside Liabilities/Total Liabilities", position(11) size(small))	
graph export ../output/connectivities_1.pdf, as(pdf) replace
graph export ../output/connectivities_1.png, as(png) replace

*Figure 8: Graph of average probabilities of default in FFunds sectors
foreach var of varlist delta_* avg_delta_assets{
	replace `var' = 100 * `var'
}
replace delta = 100*delta
twoway (line delta_insurance avg_delta_assets delta_reit delta_other qt_dt if max_bank_bhc == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end) || ///
	(line delta qt_dt if tkr == "BRO10" & qt_dt >= $charts_start & qt_dt <= $charts_end) || ///
	(line delta qt_dt if tkr == "BRO25" & qt_dt >= $charts_start & qt_dt <= $charts_end), ///
	xtitle(" ", size(medsmall)) xlabel($charts_start(4)$charts_end, angle(45) labsize(small)) ylabel(,labsize(small)) ///
	subtitle("% per year", position(11) size(small)) ///
	ytitle(" ",  height(6)) ///
	legend(order(1 "Insurance" 2 "Bank Holding Companies" 3 "REITs" 4 "Other" 5 "Top 10 Dealers" 6 "Top 11-25 Dealers") r(2) region(style(none)) size(small)) ///
	title(" ", pos(11) size(medlarge)) graphregion(color(white) style(none))
graph export ../output/delta_externals.pdf, as(pdf) replace
graph export ../output/delta_externals.png, as(png) replace

*Figure 9: Several Firm-Specific Variables for Select Large BHCs

preserve
collapse (mean) contag_index, by(tkr)
gsort -contag_index
keep if _n <= 4
levelsof(tkr), local(large_contagions) clean
restore

preserve
local to_combine
replace c = c/1000
replace w = w/1000
replace contag_index = contag_index/1000

foreach tic in `large_contagions'{
	if "`tic'" == "BAC"{
		local title = "Bank of America"
	}
	else if "`tic'" == "JPM"{
		local title = "JP Morgan Chase"
	}
	else if "`tic'" == "C"{
		local title = "Citigroup"
	}
	else if "`tic'" == "WFC"{
		local title = "Wells Fargo"
	}
	local axis_label = yq(2016, 1)
	twoway(tsline contag_index c w, yaxis(1)) || (tsline beta, yaxis(2)) if tkr == "`tic'", name(index_`tic', replace) ///
	subtitle("Billions", position(11) size(small) span) ///
	ytitle(" ", axis(1)) ytitle(" ", axis(2)) ///
	ylabel( 0 .25 .5 .75 1, axis(2) labsize(small) angle(0)) ///
	ylabel(, axis(1) labsize(small) angle(0)) ///
	xtitle(" ") title("`title'", size(small) c(black) ) ///
	xlabel($charts_start(8)$charts_end, angle(45) labsize(small)) ///
	graphregion(color(white) style(none)) ///
	legend(order(1 "Contagion Index (left axis)" 2 "Outside Assets (left axis)" 3 "Net Worth (left axis)" 4 "Financial Connectivity (right axis)") r(2) region(style(none)) symxsize(6) size(vsmall)) ///
	note("Inside/Total Liabilities", position(1) size(small) span ring(6))
	//graph play ytitle_offset_2axes
	local to_combine `to_combine' index_`tic'
}
graph close _all
grc1leg `to_combine', graphregion(color(white) style(none))
restore
graph export ../output/info_large_contagions.pdf, as(pdf) replace
graph export ../output/info_large_contagions.png, as(png) replace

*Figure 10: Network Vulnerability with Additional Costs of Bankruptcy
twoway (tsline nvi_benchmark index_w_ff_gam1 index_w_ff_gam5 index_w_ff_gam10 index_w_ff_gam15 if max_bank == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end),  ///
	xtitle(" ", size(medsmall)) xlabel($charts_start(4)$charts_end, angle(45) labsize(small)) ylabel(0 10 20 30 40 50 60, axis(1) labsize(small)) ///
	ytitle(" ",  height(6)) subtitle("%", position(11) size(small)) ///
	legend(order(1 "Benchmark NVI ({&gamma} = 0%)" 2 "{&gamma} = 1%" 3 "{&gamma} = 5%" 4 "{&gamma} = 10%" 5 "{&gamma} = 15%") r(2) region(style(none)) size(small)) ///
	title(" ", pos(11) size(medlarge)) graphregion(color(white) style(none))
graph export ../output/robustness_gamma.pdf, as(pdf) replace
graph export ../output/robustness_gamma.png, as(png) replace

*Figure 11: Network Vulnerability Under Different Balance Sheet Inside vs Outside Classificatoins
label values index_wffunds_unc_30perc big
replace index_wffunds_unc_70perc = 50 if index_wffunds_unc_70perc >= 50 & !missing(index_wffunds_unc_70perc)
twoway (tsline index_wffunds_unc_30perc index_wffunds_unc_50perc index_wffunds_unc_70perc if max_bank == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end),  ///
	xtitle(" ", size(medsmall)) xlabel($charts_start(4)$charts_end, angle(45) labsize(small)) ylabel(0 10 20 30 40, valuelabel labsize(small)) ///
	ytitle(" ",  height(6)) subtitle("%", position(11) size(small)) ///
	title(" ", pos(11) size(medlarge) color(black)) graphregion(color(white) style(none)) ///
	legend(order(1 "30% Inside Network" 2 "Benchmark NVI (50% Inside Network)" 3 "70% Inside Network") r(3) region(style(none)) size(small)) name(uncertain_robust, replace)
graph export ../output/robustness_perc_unc.pdf, as(pdf) replace

*Figure 12: Network Vulnerability and Maximum Connectivity, Different Selection Criteria
twoway (scatter beta_max_1 qt_dt if max_bank == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end, mlabgap(*5) mlabsize(tiny) mlabel(tkr) mlabpos(11) mlabangle(rvertical)) || ///
	(scatter beta_max_2 qt_dt if max_bank == 2 & qt_dt >= $charts_start & qt_dt <= $charts_end) || ///
	(scatter beta_max_3 qt_dt if max_bank == 3 & qt_dt >= $charts_start & qt_dt <= $charts_end) || ///
	(scatter beta_max_all_1 qt_dt if beta_max_all_1<1 & max_bank_all == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end, mlabgap(*5) mlabsize(tiny) mlabel(tkr) mlabpos(6) mlabangle(rvertical)), ///
	xtitle(" ", size(medsmall)) xlabel($charts_start(4)$charts_end, angle(45) labsize(small)) ylabel(.5 .6 .7 .8 .9 1,labsize(small)) ///
	ytitle(" ",  height(6)) subtitle("Inside Liabilities/Total Liabilities", position(11) size(small)) ///
	legend(order(1 "Benchmark BHC Connectivity ({&beta}{superscript:+})" 2 "Second Highest BHC Connectivity ({&beta}{superscript:+})" 3 "Third Highest BHC Connectivity ({&beta}{superscript:+})" 4 "Highest {&beta}, from Entire Sample ({&beta}{superscript:+})") r(4) region(style(none)) size(small)) ///
	title("(b)", pos(11) size(medlarge) color(black)) graphregion(color(white) style(none)) name(scatter_betas_top3, replace)
graph export ../output/scatter_betanum.pdf, as(pdf) replace

twoway (tsline nvi_benchmark if max_bank == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end) || ///
	(tsline index_w_ff_beta2  if max_bank == 2 & qt_dt >= $charts_start & qt_dt <= $charts_end) || ///
	(tsline index_w_ff_beta3 if max_bank == 3 & qt_dt >= $charts_start & qt_dt <= $charts_end) || ///
	(tsline index_w_ff_beta_all1 if max_bank_all == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end), ///
	xtitle(" ", size(medsmall)) xlabel($charts_start(4)$charts_end, angle(45) labsize(small)) ylabel(,labsize(small)) ///
	ytitle(" ",  height(6)) subtitle("%", position(11) size(small)) ylabel(0 10 20 30 40,labsize(small)) ///
	legend(order(1 "Benchmark NVI" 2 "Second Highest BHC Connectivity" 3 "Third Highest BHC Connectivity" 4 "Highest {&beta}, from Entire Sample") r(4) region(style(none)) size(small)) ///
	title("(a)", pos(11) size(medlarge)) graphregion(color(white) style(none)) name(nvi_betanum, replace)
graph export ../output/nvi_betanum.pdf, as(pdf) replace

graph combine nvi_betanum scatter_betas_top3, title(" ", pos(11) size(medlarge)) graphregion(color(white) style(none))
graph export ../output/robustness_beta_selection.pdf, as(pdf) replace

*Figure 13: Network Vulnerability Calculated with FR-Y15 Data
local y15_start = yq(2013, 1)
local y15_end = yq(2016, 1)
twoway (tsline nvi_benchmark if max_bank == 1 & qt_dt >= `y15_start' & qt_dt <= `y15_end') || ///
	(tsline nvi_y15_only if max_bank == 1 & qt_dt >= `y15_start' & qt_dt <= `y15_end') || ///
	(tsline nvi_benchmark_y15cov if max_bank_y15 == 1 & qt_dt >= `y15_start' & qt_dt <= `y15_end') || ///
	(tsline nvi_y15_only_wconnect if max_bank_y15 == 1 & qt_dt >= `y15_start' & qt_dt <= `y15_end'), ///
	xtitle(" ", size(medsmall)) xlabel(`y15_start'(4)`y15_end', angle(45) labsize(small)) ylabel(0 .5 1 1.5 2,labsize(small)) ///
	ytitle(" ",  height(6)) subtitle("%", position(11) size(small)) ///
	legend(order(1 "Benchmark NVI" 2 "w/ Only Y-15 Firms" 3 "w/ Y-15 Connectivity" 4 "w/ Y-15 Connectivity and only Y-15 Firms") r(2) region(style(none)) size(small)) ///
	title(" ", pos(11) size(medlarge)) graphregion(color(white) style(none))
graph export ../output/robustness_nvi_y15.pdf, as(pdf) replace
graph export ../output/robustness_nvi_y15.png, as(png) replace

*Figure 14: Network Vulnerability with Alternate Percentages of Uninsured Deposits Inside the Network
twoway (tsline nvi_benchmark index_wffunds_unc_depos_20in if max_bank == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end) || ///
	(tsline nvi_benchmark_y15depos if max_bank_y15depos == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end),  ///
	xtitle(" ", size(medsmall)) xlabel($charts_start(4)$charts_end, angle(45) labsize(small)) ylabel(,labsize(small)) ///
	ytitle(" ",  height(6)) subtitle("%", position(11) size(small)) ///
	legend(order(1 "Benchmark NVI (100% of Uninsured Deposits are In-Network Liabilities)" 2 "20% of Uninsured Deposits are In-Network Liabilities" 3 "Uninsured Deposits from FR-Y15") r(3) region(style(none)) size(small)) ///
	title(" ", pos(11) size(medlarge)) graphregion(color(white) style(none))
graph export ../output/robustness_y15_depos.pdf, as(pdf) replace

*Figure 15: Network Vulnerability with Extrapolated Quantities for FR-Y15 Off-Balance Sheet Items
twoway (tsline nvi_benchmark if max_bank == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end) || ///
	(tsline nvi_y15_only if max_bank_y15_samp == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end) || ///
	(tsline nvi_offbalance if max_bank_offbalance == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end), ///
	xtitle(" ", size(medsmall)) xlabel($charts_start(4)$charts_end, angle(45) labsize(small)) ylabel(,labsize(small)) ///
	ytitle(" ",  height(6)) subtitle("%", position(11) size(small)) ///
	legend(order(1 "Benchmark NVI" 2 "w/ BHCs only" 3 "w/ BHCs only, Estimating Off-Balance Sheet Magnitudes") r(3) region(style(none)) size(small)) ///
	title(" ", pos(11) size(medlarge)) graphregion(color(white) style(none))
graph export ../output/robustness_offbalance.pdf, as(pdf) replace


/*
*Figure 12: Network Vulnerability and Maximum Connectivity Under Different Selection Schemes for Maximum Connectivity
*(A)
twoway (tsline nvi_benchmark if max_bank == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end) || ///
	(tsline index_w_ff_beta_all1  if max_bank_all == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end) || ///
	(tsline index_w_ff_beta_full1  if max_bank_fullsample == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end),  ///
	xtitle(" ", size(medsmall)) xlabel($charts_start(4)$charts_end, angle(45) labsize(small)) ylabel(0 10 20 30 40 50 60,labsize(small)) ///
	ytitle(" ",  height(6)) subtitle("%", position(11) size(small)) ///
	legend(order(1 "Benchmark" 2 "Max Connectivity from Entire Y9C Sample" 3 "Max Connectivity from Constant BHC Sample") r(3) region(style(none)) size(small)) ///
	title("(a)", pos(11) size(medlarge) color(black)) graphregion(color(white) style(none)) name(nvi_betas_selection, replace)
graph export ../output/robustness_beta_selection.pdf, as(pdf) replace

*(B)
twoway (scatter beta_max_1 qt_dt if max_bank == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end, mlabgap(*3) mlabsize(tiny) mlabel(tkr) mlabpos(6) mlabangle(rvertical)) || ///
	(scatter beta_max_all_1 qt_dt if beta_max_all_1<1 & max_bank_all == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end, mlabgap(*5) mlabsize(tiny) mlabel(tkr) mlabpos(6) mlabangle(rvertical)) || ///
	(scatter beta_max_fullsample_1 qt_dt if max_bank_fullsample == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end, mlabgap(*5) mlabsize(tiny) mlabel(tkr) mlabpos(11) mlabangle(rvertical)), ///
	xtitle(" ", size(medsmall)) xlabel($charts_start(4)$charts_end, angle(45) labsize(small)) ylabel(.5 .6 .7 .8 .9 1,labsize(small)) ///
	ytitle(" ",  height(6)) subtitle("Inside Liabilities/Total Liabilities", position(11) size(small)) ///
	legend(order(1 "Benchmark" 2 "From Entire Y9C Sample" 3 "From Constant BHC Sample") r(3) region(style(none)) size(small)) ///
	title("(b)", pos(11) size(medlarge) color(black)) graphregion(color(white) style(none)) name(scatter_betas_selection, replace)
graph export ../output/scatter_beta_selection.pdf, as(pdf) replace

graph combine nvi_betas_selection scatter_betas_selection, title(" ") graphregion(color(white) style(none))
graph export ../output/robustness_beta_selection.pdf, as(pdf) replace
graph export ../output/robustness_beta_selection.png, as(png) replace
*/

*Table 1: Select Firm-Specific Model Variables for Large BHCs, 2016Q4.

preserve
replace c = c/1000
replace w = w/1000
replace contag_index = contag_index/1000
keep if qt_dt == yq(2016, 4)
keep if keep_obs == 1

gsort -contag_index
replace name = regexr(name, "&", "")
replace nm_short = regexr(nm_short, "&", "")
replace nm_short = string(_n, "%02.0f" ) + "." + nm_short
estpost tabstat contag_index beta c w, by(nm_short) nototal
local lab_actual
foreach lab in `e(labels)'{
disp "`lab'"
	if regexm("`lab'", "[0-9]*\."){
		global temp = regexr("`lab'", "[0-9]*\.", "")
		local lab_actual `"`lab_actual' `"$temp"'"'
	}
	else{
		local lab_actual `"`lab_actual' `lab'"'
	}
}

esttab using ../output/large_instit_breakdowns.tex, ///
	nonum label replace title(" ") noobs booktabs ///
	cells("contag_index(fmt(2)) beta(fmt(2)) c(fmt(2)) w(fmt(2))") coeflabels(`lab_actual') collabels("Contagion Index" "Financial Connectivity" "Outside Assets" "Net Worth")
restore

*Table 2: See sum_stats_wholesample.do 

*Table 3: See sum_stats_wholesample.do 

*Table 4: See sum_stats_wholesample.do 

*Table 5: See sum_stats_wholesample.do 

*Table 6: Not automated

*Table 7: Not automated

*Table 8: Categorization portion not automated. For % of Sector Assets portion, see ffunds.do 

*Table 9: Not automated

*Other paper statistics

*% of domestic financial system covered
qui:sum coverage_included
display "Share of Domestic Financial System Directly Covered: `r(max)'"
//keep if coverage_included == `r(max)'
sum qt_dt if coverage_included == `r(max)'
//local max_qt = qt_dt[1]
display "Max achieved at `r(max)'"
//restore
preserve
qui:sum nvi_benchmark if qt_dt == yq(2009, 1) & max_bank == 1
//local nvi = nvi_benchmark[1]
display "First Quarter of 2009 NVI: `r(max)'"
qui:sum connectivity_comp if qt_dt == yq(2002, 1) | qt_dt == yq(2007, 4)
local perc_diff = (`r(max)' - `r(min)')/`r(min)'
display "% Change in Connectivity Multiplier over 2002Q1-2007Q4: `perc_diff'"
qui:sum nvi_benchmark
display "NVI Benchmark between `r(min)' and `r(max)'"
qui:sum nvi_benchmark if qt_dt == yq(2008, 4)
display "NVI Benchmark in 2008Q4: `r(max)'"

forval i=1(1)4{
	qui:sum nvi_benchmark if qt_dt == yq(2009, `i')
	display "NVI Benchmark in 2009Q`i': `r(max)'"
}

restore

*Exploratory charts not yet in the paper

/* *Average default probabilities, scaled by their contribution to system-wide NVI component */
/* twoway (line contribution_insurance contribution_bhc contribution_reit contribution_other contribution_dealer qt_dt if max_bank_bhc == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end), /// */
/* 	xtitle(" ", size(medsmall)) xlabel($charts_start(4)$charts_end, angle(45) labsize(small)) ylabel(,labsize(small)) /// */
/* 	subtitle("% per year", position(11) size(small)) /// */
/* 	ytitle(" ",  height(6)) /// */
/* 	legend(order(1 "Insurance" 2 "Bank Holding Companies" 3 "REITs" 4 "Other" 5 "Top 25 Dealers") r(2) region(style(none)) size(small)) /// */
/* 	title(" ", pos(11) size(medlarge)) graphregion(color(white) style(none)) */
/* graph export ../output/NVI_contributions.pdf, as(pdf) replace */
/* graph export ../output/NVI_contributions.png, as(png) replace */

/* *Average default probabilities, scaled by their contribution to system-wide NVI component */
/* twoway (line contribution_insurance contribution_bhc contribution_reit contribution_kmv_other contribution_nonkmv_other contribution_dealer qt_dt if max_bank_bhc == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end), /// */
/* 	xtitle(" ", size(medsmall)) xlabel($charts_start(4)$charts_end, angle(45) labsize(small)) ylabel(,labsize(small)) /// */
/* 	subtitle("% per year", position(11) size(small)) /// */
/* 	ytitle(" ",  height(6)) /// */
/* 	legend(order(1 "Insurance" 2 "Bank Holding Companies" 3 "REITs" 4 "KMV-Included Other" 5 "Not KMV-Included Other" 6 "Top 25 Dealers") r(3) region(style(none)) size(small)) /// */
/* 	title(" ", pos(11) size(medlarge)) graphregion(color(white) style(none)) */
/* graph export ../output/NVI_contributions_othersplit.pdf, as(pdf) replace */
/* graph export ../output/NVI_contributions_othersplit.png, as(png) replace */

/* *Scaling used above... */
/* twoway (line weight_insurance weight_bhc weight_reit weight_other weight_dealer qt_dt if max_bank_bhc == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end), /// */
/* 	xtitle(" ", size(medsmall)) xlabel($charts_start(4)$charts_end, angle(45) labsize(small)) ylabel(,labsize(small)) /// */
/* 	subtitle("Proportion", position(11) size(small)) /// */
/* 	ytitle(" ",  height(6)) /// */
/* 	legend(order(1 "Insurance" 2 "Bank Holding Companies" 3 "REITs" 4 "Other" 5 "Top 25 Dealers") r(2) region(style(none)) size(small)) /// */
/* 	title(" ", pos(11) size(medlarge)) graphregion(color(white) style(none)) */
/* graph export ../output/NVI_wts.pdf, as(pdf) replace */
/* graph export ../output/NVI_wts.png, as(png) replace */

/* *Time-varying scatter of betas */
/* bys qt_dt: egen beta_mean = mean(beta) */
/* bys qt_dt: egen beta_25 = pctile(beta), p(25) */
/* bys qt_dt: egen beta_75 = pctile(beta), p(75) */

/* twoway (scatter beta qt_dt if qt_dt >= $charts_start & qt_dt <= $charts_end, msize(vtiny)) || /// */
/* 	(line beta_mean beta_25 beta_75 qt_dt if qt_dt >= $charts_start & qt_dt <= $charts_end) || /// */
/* 	(line beta_max_1 qt_dt if max_bank == 1 & qt_dt >= $charts_start & qt_dt <= $charts_end), /// */
/* 	ytitle("",  height(6)) xtitle("", size(medsmall)) ylabel(0 .1 .2 .3 .4 .5 .6 .7 .8 .9) /// */
/* 	title(" ", pos(11) size(medlarge)) graphregion(color(white) style(none)) /// */
/* 	legend(order(1 "Actual {&beta}s" 2 "Mean of {&beta}" 3 "25th Percentile" 4 "75th Percentile" 5 "{&beta}{superscript:+}") r(3) region(style(none)) size(small)) /// */
/* 	xlabel($charts_start(4)$charts_end, angle(45) labsize(medsmall)) /// */
/* 	subtitle("Inside Liabilities/Total Liabilities", position(11) size(small)) */
/* graph export ../output/beta_dist.pdf, as(pdf) replace */
/* graph export ../output/beta_dist.png, as(png) replace */

