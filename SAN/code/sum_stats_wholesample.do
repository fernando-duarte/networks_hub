*********************************************************************************
*
*	Summary Statistics for Y9C Sample
*	Collin Jones
*	Purpose: Construct some summary statistics tables of BHC assets and liabilities.
*		This new code runs the major tables on the entire Y9C subsample, with double-checked
*		consistency with numbers from rest of the analysis.
*		Also creates tables for sector-wide default probability samples.
*		Note: This code is called by Model_series_processing at an appropriate stage of the merging process
*
*********************************************************************************


**********************************************************************************************************
*
*		Summary statistics of relative size of different Y9C asset/liabilities categories (Tables 2 and 3)
*
*********************************************************************************************************

*Liabilities -- Easy
gen help_liab_in_trading = BHCKG209 + BHCKG210 + BHCKG211 + BHCKF624
label var help_liab_in_trading "Trading Liabilities"

gen help_liab_in_deriv = BHCK3547
label var help_liab_in_deriv "Derivatives"

gen help_liab_in_repo_ffunds = BHDMB993 + BHCKB995
label var help_liab_in_repo_ffunds "Repos and Fed Funds"

gen help_liab_in_depos = tot_uninsur_depos_y9c
label var help_liab_in_depos "Uninsured Domestic Deposits"

gen help_liab_in_longterm_debt = 0.5*BHCK2333
label var help_liab_in_longterm_debt "Longer Term Debt"

gen help_liab_in_subor_debt = .5*BHCKC699
label var help_liab_in_subor_debt "Subordinated Debt"

gen help_liab_in_other = .5 * BHCKB557 + .5 * BHCKB984
label var help_liab_in_other "Other"

gen help_liab_in_short_debt = .5 * BHCK2309 + .5 * BHCK2332
label var help_liab_in_short_debt "Short Term Debt"

gen help_liab_in_time_depos = .5*BHOD6648 + .5*BHOD2604
label var help_liab_in_time_depos "Time Balance Deposits"

gen help_liab_in_transac_depos = .5*BHOD3187
label var help_liab_in_transac_depos "Transaction Account Deposits"

gen help_liab_in_noninter_depos = .5*BHOD3189
label var help_liab_in_noninter_depos "Noninterest Bearing Deposits"

gen help_liab_in_savings_depos = .5*BHOD2389
label var help_liab_in_savings_depos "Savings Deposits"

gen help_liab_out_foreign_depos = BHFN6631+BHFN6636
label var help_liab_out_foreign_depos "Foreign Deposits"

gen help_liab_out_depos_insur = tot_insur_depos
label var help_liab_out_depos_insur "Insured Domestic Deposits"

gen help_liab_out_longterm_debt = 0.5*BHCK2333
label var help_liab_out_longterm_debt "Longer Term Debt"

gen help_liab_out_subor_debt = BHCK4062 + .5*BHCKC699
label var help_liab_out_subor_debt "Subordinated Debt"

gen help_liab_out_time_depos = .5*BHOD6648 + .5*BHOD2604
label var help_liab_out_time_depos "Time Balance Deposits"

gen help_liab_out_transac_depos = .5*BHOD3187
label var help_liab_out_transac_depos "Transaction Account Deposits"

gen help_liab_out_noninter_depos = .5*BHOD3189
label var help_liab_out_noninter_depos "Noninterest Bearing Deposits"

gen help_liab_out_savings_depos = .5*BHOD2389
label var help_liab_out_savings_depos "Savings Deposits"

gen help_liab_out_other = BHCK3049 + .5 * BHCKB557 + .5 * BHCKB984
label var help_liab_out_other "Other"

gen help_liab_out_short_debt = .5 * BHCK2309 + .5 * BHCK2332
label var help_liab_out_short_debt "Short Term Debt"

*Assets -- Hard

gen help_asset_in_deriv = .5 * BHCM3543
label var help_asset_in_deriv "Derivatives"

gen help_asset_in_agencymbs = .5*BHCKG316 + .5*BHCKG319
label var help_asset_in_agencymbs "Agency MBS"

gen help_asset_in_depos_inter = BHCK0397 + .5 * BHCK0395
label var help_asset_in_depos_inter "Interest Bearing Deposits"

gen help_asset_in_loans = .5 * BHCKF618
label var help_asset_in_loans "Loans"

gen help_asset_in_other = BHCK2130 + .5 * (BHCK5507 + BHCK6438 + BHCKA519 + BHCKA520 + BHCKB026 + ///
	BHCKB556 + BHCKK201 + BHCKK202 + BHCKK270 + BHCK2168)

gen help_asset_in_misctrade = .5 * BHCM3541
label var help_asset_in_misctrade "Other Trading Assets"

gen help_asset_in_goodwill = .5 * BHCK3163
label var help_asset_in_goodwill "Goodwill"

label var help_asset_in_other "Other"

gen help_asset_in_othermbs = .5*(BHCKG308 + BHCKG311 + BHCKG381 + BHCKK146 + BHCKK149 + BHCKK154 + BHCKK157 + BHCKK198)
label var help_asset_in_othermbs "Other MBS"

gen help_asset_in_othersec = .5*BHCKA511
label var help_asset_in_othersec "Other Securities"

gen help_asset_in_priv_lab_abs = BHCKC026 + BHCKC027 + BHCKG336 + BHCKG339 + BHCKG340 + ///
	BHCKG343 + BHCKG344 + BHCKG347 + BHCKG383 + BHCKG384 + BHCKG385
label var help_asset_in_priv_lab_abs "Private Label ABS"

gen help_asset_in_priv_lab_mbs = BHCKG320 + BHCKG323
label var help_asset_in_priv_lab_mbs "Private Label MBS"

gen help_asset_in_repo_ffunds = BHCKB989 + BHDMB987
label var help_asset_in_repo_ffunds "Repos and Fed Funds"

gen help_asset_out_agencymbs = BHCKG300 + BHCKG303 + BHCKG304 + BHCKG307 + BHCKG312 + BHCKG315 + BHCKG379 + BHCKG380 + ///
	BHCKK142 + BHCKK145 + BHCKK150 + BHCKK153 + BHCKK197 + .5*BHCKG316 + .5*BHCKG319
label var help_asset_out_agencymbs "Agency MBS"

gen help_asset_out_deriv = .5 * BHCM3543
label var help_asset_out_deriv "Derivatives"

gen help_asset_out_depos_inter = .5 * BHCK0395
label var help_asset_out_depos_inter "Interest Bearing Deposits"

gen help_asset_out_loans = BHCK5369 + BHCKB529 + BHCKF610 + BHCKF614 + BHCKF615 + BHCKF616 + BHCKK199 + BHCKK210 + .5 * BHCKF618
label var help_asset_out_loans "Loans"

gen help_asset_out_noninters_depos = BHCK0081
label var help_asset_out_noninters_depos "Noninterest Bearing Deposits"

gen help_asset_out_other =  BHCK3164 + BHCK2145 + BHCK2148 + BHCK3656 + BHCK2150 + ///
	.5 * (BHCK5507 + BHCK6438 + BHCKA519 + BHCKA520 + BHCKB026 + ///
	BHCKB556 + BHCKK201 + BHCKK202 + BHCKK270 + BHCK2168)
label var help_asset_out_other "Other"

gen help_asset_out_misctrade = .5 * BHCM3541
label var help_asset_out_misctrade "Other Trading Assets"

gen help_asset_out_goodwill = .5 * BHCK3163
label var help_asset_out_goodwill "Goodwill"

gen help_asset_out_othermbs = .5*(BHCKG308 + BHCKG311 + BHCKG381 + BHCKK146 + BHCKK149 + BHCKK154 + BHCKK157 + BHCKK198)
label var help_asset_out_othermbs "Other MBS"

gen help_asset_out_othersec = BHCK1737 + BHCK1741 + BHCK1742 + BHCK1746 + BHCK1752 + BHCKG386 + .5*BHCKA511
label var help_asset_out_othersec "Other Securities"

gen help_asset_out_gov_debt = BHCK0211 + BHCK1287 + BHCK1289 + BHCK1293 + BHCK1294 + BHCK1298 + BHCK8496 + BHCK8499 + BHCM3531 + BHCM3532 + BHCM3533
label var help_asset_out_gov_debt "State, Treasury, and Agency Debt"

*Sanity-checking of some of the numbers. Making sure that totals for in, out, and all assets/liabilities match what they were in the rest of the analysis
egen liab_in_test = rowtotal(help_liab_in*)
egen liab_out_test = rowtotal(help_liab_out*)
gen liab_total_test = liab_in_test + liab_out_test

gen liab_in_actuals = liab_in + .5*liab_in_unc
gen liab_out_actual = BHCK2948 - liab_in_actuals

egen asset_in_test = rowtotal(help_asset_in*)
egen asset_out_test = rowtotal(help_asset_out*)
gen asset_total_test = asset_in_test + asset_out_test

gen asset_in_actual = asset_in + .5*asset_in_unc
gen asset_out_actual = BHCK2170 - asset_in_actual

*Creating some summary statistics tables
label var asset_in_test "% of BHC Assets"
label var asset_out_test "% of BHC Assets"

label var liab_in_test "% of BHC Liabilities"
label var liab_out_test "% of BHC Liabilities"

keep if qt_dt == $snapshot_date
*Want to save current labels for later...
foreach v of varlist *{
	capture local l`v': variable label `v'
}
collapse (sum) help_* *_test

*Re-attached labels after collapse
foreach v of varlist*{
	label var `v' "`l`v''"
}

*Generating stats for % attributable to each category generated above
foreach side in asset liab{
	foreach var of varlist help_`side'_in*{
		replace `var' = 100 * `var' / `side'_in_test
	}

	*If it's going to be display on the table with <1%, then loop it into Other
	foreach var of varlist help_`side'_in*{
		if `var'[1] < 1{
			replace help_`side'_in_other = help_`side'_in_other + `var'
			drop `var'
		}
	}

	*Same for `outside' liabs and assets
	foreach var of varlist help_`side'_out*{
		replace `var' = 100 * `var' / `side'_out_test
	}

	foreach var of varlist help_`side'_out*{
		if `var'[1] < 1{
			replace help_`side'_out_other = help_`side'_out_other + `var'
			drop `var'
		}
}

replace `side'_in_test = 100*`side'_in_test/`side'_total_test
replace `side'_out_test = 100*`side'_out_test/`side'_total_test

}

*Pivoting things for a better-formatted table
keep help* asset_in_test asset_out_test liab_in_test liab_out_test
xpose, clear varname
gen side = ""
replace side = "Assets" if strpos(_varname, "asset")
replace side = "Liabilities" if strpos(_varname, "liab")
gen in_out = ""
replace in_out = "In" if strpos(_varname, "_in_") & !strpos(_varname, "test")
replace in_out = "Out" if strpos(_varname, "_out_") & !strpos(_varname, "test")

*Making sure those labels we made above make it onto the table...
levelsof _varname, clean local(vars)
gen component = ""
foreach v in `vars'{
	replace component = "`l`v''" if _varname == "`v'"
}

*Ordering so that larger assets/liabilities are displayed first (with other always last)
gsort -v1
gen sorter = _n
replace sorter = 100 if component == "Other"
replace sorter = 1000 + _n if strpos(component, "%")
labmask sorter, values(component)

*Inside Liabilities
estpost tabstat v1 if (side == "Liabilities" & in_out == "In") | (_varname == "liab_in_test"), columns(statistics) by(sorter) nototal
esttab using ../output/liab_in.tex, coeflabels(`e(labels)') main(mean 2) booktabs nonum replace title(" ") mtitle("BHC Liabilities Inside Financial System (%)") noobs nonotes

*Outside Liabilities
estpost tabstat v1 if (side == "Liabilities" & in_out == "Out") | (_varname == "liab_out_test"), columns(statistics) by(sorter) nototal
esttab using ../output/liab_out.tex, coeflabels(`e(labels)') main(mean 2) booktabs nonum label replace title(" ") mtitle("BHC Liabilities Outside Financial System (%)") noobs nonotes

*Inside Assets
estpost tabstat v1 if (side == "Assets" & in_out == "In") | (_varname == "asset_in_test"), columns(statistics) by(sorter) nototal
esttab using ../output/asset_in.tex, coeflabels(`e(labels)') main(mean 2) booktabs nonum label replace title(" ") mtitle("BHC Assets Inside Financial System (%)") noobs nonotes

*Outside Assets
estpost tabstat v1 if (side == "Assets" & in_out == "Out") | (_varname == "asset_out_test"), columns(statistics) by(sorter) nototal
esttab using ../output/asset_out.tex, coeflabels(`e(labels)') main(mean 2) booktabs nonum label replace title(" ") mtitle("BHC Assets Outside Financial System(%) ") noobs nonotes


**********************************************************************************************************
*
*	Shows firms in sample for sector-wide default probabilities in our snapshot period. Also displays their asset weighting.
*		Tables 4 and 5
*
*********************************************************************************************************
use ../output/firm_list_final_wweights, clear

tempvar sorter
gsort sector -weight tkr, gen(`sorter')
replace company_name = subinstr(company_name, "&", "", 10)

drop nm_short
rename company_name nm_short
rename total_liabilities_amt p_bar
rename total_assets_current_amt assets
rename edf01 delta
rename edf01_mean delta_alt

gen w = assets - p_bar
foreach var in $simulation_data_vars{
	capture confirm variable `var', exact
	disp _rc
	if _rc != 0{
		gen `var' = .
		}
	}

levelsof sector, local(sectors)
levelsof qt_dt, local(quarters)
foreach sec in `sectors'{
	/* keep if sector == "`sec'" & weight > 0 */
	/* keep if weight > 0 */
	/* keep company_name tkr total_assets_current_amt edf01 edf01_mean total_liabilities_amt qt_dt */
	export excel $simulation_data_vars using ../temp/node_stats_forsimulation_all if sector == "`sec'" & weight > 0, ///
	  sheet("`sec'", replace) firstrow(variables)
}

keep if qt_dt == $snapshot_date
expand 3 if sector != sector[_n+1]
gsort sector -weight tkr

//nm_short tkr p_bar assets c b delta beta w

replace nm_short = "Number of Firms in Sample" if tkr == tkr[_n+1] & tkr != tkr[_n+2]
replace weight = count if nm_short == "Number of Firms in Sample"
replace nm_short = "Weighting from Rest of Sample" if sector != sector[_n+1]
replace weight = 0 if nm_short == "Weighting from Rest of Sample"
bys sector: egen temp = total(weight) if weight < .01 | nm_short == "Weighting from Rest of Sample"
replace weight = temp if nm_short == "Weighting from Rest of Sample"
keep if weight > .01 | inlist(sector, "Top 10 Dealers", "Top 11-25 Dealers")
gen sorter = _n
replace sorter = 1000 if nm_short == "Number of Firms in Sample"
replace sorter = 1001 if nm_short == "Weighting from Rest of Sample"
drop if nm_short == "Weighting from Rest of Sample" & inlist(sector, "Top 10 Dealers", "Top 11-25 Dealers")
labmask sorter, values(nm_short)

estpost tabstat weight if sector == "Insurance", columns(statistics) by(sorter) nototal
esttab using ../output/insurance_sample.tex, coeflabels(`e(labels)') main(mean 2) booktabs nonum replace title(" ") mtitle("Asset Weighting") noobs nonotes not

estpost tabstat weight if sector == "REITs", columns(statistics) by(sorter) nototal
esttab using ../output/reit_sample.tex, coeflabels(`e(labels)') main(mean 2) booktabs nonum replace title(" ") mtitle("Asset Weighting") noobs nonotes not

estpost tabstat weight if sector == "Other", columns(statistics) by(sorter) nototal
esttab using ../output/other_sample.tex, coeflabels(`e(labels)') main(mean 2) booktabs nonum replace title(" ") mtitle("Asset Weighting") noobs nonotes not

estpost tabstat weight if sector == "Top 10 Dealers", columns(statistics) by(sorter) nototal
esttab using ../output/dealers_top10_sample.tex, coeflabels(`e(labels)') main(mean 2) booktabs nonum replace title(" ") mtitle("Asset Weighting") noobs nonotes not

estpost tabstat weight if sector == "Top 11-25 Dealers", columns(statistics) by(sorter) nototal
esttab using ../output/dealers_top25_sample.tex, coeflabels(`e(labels)') main(mean 2) booktabs nonum replace title(" ") mtitle("Asset Weighting") noobs nonotes not
