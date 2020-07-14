
*********************************************************************************
*
*			Average Probabilities of Default Within Sectors: Firm Assets
*				Uses mkmv_id to do matching between Y9C and KMV data
*				Uses KMV's own classifications to get data from as many sector firms as possible
*
*********************************************************************************

*Saving mkmv_id of firms in Y9C. Don't want to include any firms from Y9C when computing default probs.
*		From what I've seen, big items this excludes (though not for whole sample): Goldman Sachs, Charles Schwab, E-Trade, AIG, Metlife
clear all
clear matrix
clear mata
set maxvar 10000

use ../temp/Contagion_preprocessing.dta, clear
drop if FGN_CALL_FAM_ID != 0 & entity != 3232316 & entity != 1132449
keep if !missing(prob_default_phys)
bys entity: egen count_q 	= count(qt_dt) if qt_dt <= $charts_end & qt_dt >= $charts_start

*All mkmv_id in any quarter. Not currently used.
qui:levelsof mkmv_id, local(y9c_sample) clean
global y9c_sample `y9c_sample'

*Separate locals for kmv_id samples in EACH quarter. So firms are included for quarters when not in Y9C.
qui:levelsof qt_dt, local(y9c_quarters)
foreach qt in `y9c_quarters'{
	qui:levelsof mkmv_id if qt_dt == `qt', local(q_`qt') clean
	qui:levelsof mkmv_id if qt_dt == `qt' & count_q == $charts_end - $charts_start + 1, local(q_fullsamp_`qt') clean
}

*Restricting ourselves to mkmvind_name_desc sector labels that fit what we're looking for
use ../temp/CUSIP_MKMVID.dta, clear
duplicates drop mkmv_id, force
merge 1:m mkmv_id using ../temp/DFS.dta

//merge 1:m mkmv_id using ../temp/DFS_nodrop.dta
keep if mkmvind_name_desc == "SECURITY BROKERS & DEALERS" | mkmvind_name_desc == "INVESTMENT MANAGEMENT" ///
  | mkmvind_name_desc == "INSURANCE - LIFE" | mkmvind_name_desc == "INSURANCE - PROP/CAS/HEALTH" | ///
  mkmvind_name_desc == "REAL ESTATE INVESTMENT TRUSTS" | mkmvind_name_desc == "FINANCE NEC" | mkmvind_name_desc == "FINANCE COMPANIES"
keep if _merge == 3

drop _merge
merge 1:1 mkmv_id qt_dt using ../temp/DFS8.dta
keep if _merge == 3 | _merge == 1
drop _merge
*Only concerned with American companies
keep if kmv_country_id == "USA" | kmv_country_id == "BRU" | kmv_country_id == "CYM"
keep if qt_dt >= $charts_start

*drop if negative net worth
drop if total_assets_current_amt - total_liabilities_amt < 0

save ../temp/DFS_merged, replace

bys mkmv_id (qt_dt): ipolate total_assets_current_amt qt_dt, gen(total_assets_current_amt_ip)
replace total_assets_current_amt = total_assets_current_amt_ip
drop total_assets_current_amt_ip

*Some manual hard-coding away of anomalies. These include:
replace total_assets_current_amt = 0 if mkmv_id == "N04946"
replace total_assets_current_amt = 0 if mkmv_id == "313586"
replace total_assets_current_amt = 0 if mkmv_id == "313400"
replace total_assets_current_amt = 0 if mkmv_id == "441815"
replace total_assets_current_amt = 0 if mkmv_id == "590188" & qt_dt >= yq(2009, 1)

*Holding onto the mkmv_id of our difference sector samples
levelsof mkmv_id if mkmvind_name_desc == "SECURITY BROKERS & DEALERS", local(dealer_kmv) clean
levelsof mkmv_id if mkmvind_name_desc == "INSURANCE - LIFE" | mkmvind_name_desc == "INSURANCE - PROP/CAS/HEALTH" | tkr == "AIG" | tkr == "ERIE", local(insurance_kmv) clean
levelsof mkmv_id if mkmvind_name_desc == "INVESTMENT MANAGEMENT", local(mutual_kmv) clean
levelsof mkmv_id if mkmvind_name_desc == "REAL ESTATE INVESTMENT TRUSTS", local(reit_kmv) clean
levelsof mkmv_id if mkmvind_name_desc == "INVESTMENT MANAGEMENT" | mkmvind_name_desc == "FINANCE NEC" | mkmvind_name_desc == "FINANCE COMPANIES", local(other_kmv) clean
//levelsof mkmv_id if mkmvind_name_desc == "BANKS AND S&LS", local(reit_kmv) clean

local tags insurance reit other
foreach tag in `tags'{
	global `tag'_vars
	global `tag'_vars_full
	global `tag'_kmv ``tag'_kmv'
	foreach tkr in $`tag'_kmv{
		global `tag'_vars $`tag'_vars weights`tkr'
		global `tag'_vars_full $`tag'_vars_full full_weights`tkr'
	}
}

keep if qt_dt >= $charts_start

*Sector labelling
gen sector = ""
replace sector = "Broker Dealers" if mkmvind_name_desc == "SECURITY BROKERS & DEALERS"
replace sector = "Insurance" if mkmvind_name_desc == "INSURANCE - LIFE" | mkmvind_name_desc == "INSURANCE - PROP/CAS/HEALTH" | tkr == "AIG" | tkr == "ERIE"
replace sector = "Mutual" if mkmvind_name_desc == "INVESTMENT MANAGEMENT"
replace sector = "REITs" if mkmvind_name_desc == "REAL ESTATE INVESTMENT TRUSTS"
replace sector = "Other" if mkmvind_name_desc == "INVESTMENT MANAGEMENT" | mkmvind_name_desc == "FINANCE NEC" | mkmvind_name_desc == "FINANCE COMPANIES"

*For troubleshooting
preserve
save ../temp/firm_list, replace
restore

*********************************************************************************
*
*			Top 1-10 and Top 11-25 Dealer Aggregate Nodes.
*				Import c, assets, beta from python script output
*				Compute average probability by assets for these firms
*
*********************************************************************************
preserve

import delimited ../input/all.csv, varnames(1) clear
gen date = date(v1, "YMD")
gen qt_dt = qofd(date)
format qt_dt %tq

save ../temp/focus_consolidated, replace

restore

frame create Y9C
frame change Y9C
use ../temp/Contagion_preprocessing, clear
drop if FGN_CALL_FAM_ID != 0 & entity != 3232316 & entity != 1132449
keep if !missing(prob_default_phys)
drop if no_call == 1
drop if tkr == "IFIN"
gen non_bd = BHCK2170 - BHCKC252

/* frame create Y9C_agg */
frame put if BHCKC252 > 0, into(Y9C_agg)
frame change Y9C_agg
collapse (sum) bd_subsid_all = BHCKC252, by(qt_dt)
frame change default

tempname bd_node_frame
/* local bd_node_frame temp */
frame put if sector == "Broker Dealers", into(`bd_node_frame')
frame `bd_node_frame'{
	/* keep if sector == "Broker Dealers" */
	gsort qt_dt -total_assets_current_amt
	by qt_dt: gen rank = _n
	keep if rank <= 25
	gen name			= ""
	replace name		= "Top 10 Dealers" if rank <= 10
    replace name		= "Top 11-25 Dealers" if rank >= 11 & rank <= 25

	*Average default prob by assets, using V9 and V8 of KMV EDF. Will sum these across quarters
    /* bys qt_dt name: egen assets_total	= sum(total_assets_current_amt) */
	/* gen weight							= (total_assets_current_amt/assets_total) */
	gen assets_prob = total_assets_current_amt

	*We're going to hold onto the total assets and liabilities of 'dealer' firms that
	*	show up in the FR-Y9C. Not doing anything with this now, but at some point we
	*	may want to remove these assets/liabilities
	gen to_reduce_ass	= 0
	gen to_reduce_liab	= 0

	frlink 1:1 qt_dt mkmv_id, frame(Y9C)
    frget BHCKC252 non_bd, from(Y9C)
	replace to_reduce_ass		= total_assets_current_amt if ~missing(Y9C)
	replace to_reduce_liab		= total_liabilities_amt if ~missing(Y9C)

	if inlist($bd_config, 1, 3){
    	/* replace weight = 0 if ~missing(Y9C) */
		replace assets_prob = 0 if ~missing(Y9C)
	}

	if $bd_config == 2{
		/* replace weight = weight + BHCKC252/1000 / assets_total if ~missing(Y9C) */
		replace assets_prob = assets_prob + BHCKC252/1000 if ~missing(Y9C)
	}

    /* tempvar weight_resum */
    bys qt_dt name: egen assets_total = sum(assets_prob)
    gen weight = assets_prob/assets_total

	gen prob_sum						= edf01 * weight
	gen prob_sum_8						= edf01_8 * weight

	/* foreach qt in `y9c_quarters'{ */
	/* 	replace to_reduce_ass = total_assets_current_amt if (qt_dt == `qt') & strpos(" `q_`qt'' ", " " + mkmv_id + " ") */
	/* 	*replace total_assets_current_amt = 0 if (qt_dt == `qt') & strpos(" `q_`qt'' ", " " + mkmv_id + " ") */
	/* 	replace to_reduce_liab = total_liabilities_amt if (qt_dt == `qt') & strpos(" `q_`qt'' ", " " + mkmv_id + " ") */
	/* 	*replace total_assets_current_amt = 0 if (qt_dt == `qt') & strpos(" `q_`qt'' ", " " + mkmv_id + " ") */
	/* } */

	replace sector = name
	//keep sector name true_weight company_name
	*Nice for troubleshooting
    save ../temp/broker_nodes_list_config$bd_config, replace
    save ../temp/broker_nodes_list_keep_config$bd_config, replace

	collapse (sum) BHCKC252 to_reduce_ass to_reduce_liab prob_sum prob_sum_8 non_bd (lastnm) assets_total, by(qt_dt name)
	frlink m:1 qt_dt, frame(Y9C_agg)
	frget bd_subsid_all, from(Y9C_agg)
	format qt_dt %tq

	merge 1:1 name qt_dt using ../temp/focus_consolidated
	keep if _merge == 3
	drop _merge

	*Making new data file compatible with other FR-Y9C stuff, so we can just append rows
	foreach var of varlist bhck*{
    	rename `var' `=upper("`var'")'
	}

	*Assigning some identifier information to these nodes, so they are processed correctly
	*	later on
	gen tkr				= "BRO10" if name == "Top 10 Dealers"
	replace tkr			= "BRO25" if name == "Top 11-25 Dealers"
	gen nm_short		= name
	gen entity			= 123 if name == "Top 10 Dealers"
	replace entity		= 456 if name == "Top 11-25 Dealers"
	rename prob_sum prob_default_phys
	rename prob_sum_8 prob_default_phys8
	replace assets_total = 1000*assets_total

	if $bd_config == 1{
		replace BHCK2170 = BHCK2170 - BHCKC252
	}
	else if $bd_config == 2{
		replace BHCK2170 = BHCK2170 + non_bd
	}
	else if $bd_config == 3{
		tempvar focus_all
		bys qt_dt (name): egen `focus_all' = sum(BHCK2170)
		replace BHCK2170 = BHCK2170 - (BHCK2170 / `focus_all') * bd_subsid_all
	}
	label var assets_total "KMV assets, EDF sample"
	save ../temp/broker_nodes_config$bd_config.dta, replace
}

drop if sector == "Broker Dealers"

keep qt_dt mkmv_id total_assets_current_amt total_liabilities_amt edf01 edf01_8 edf01_mean

*bit more intuitive name, now
rename total_assets_current_amt weights

reshape wide weights total_liabilities_amt edf01 edf01_8 edf01_mean, i(qt_dt) j(mkmv_id) str
foreach var of varlist weights*{
	gen full_`var' = `var'
}

*Will use FFunds data to generate some coverage statistics of different sectors
merge 1:1 qt_dt using ../temp/ffunds.dta
keep if _merge == 3
drop _merge
local dupes
format qt_dt %tq

foreach tag in `tags'{
	foreach tkr in $`tag'_kmv{
		foreach qt in `y9c_quarters'{
			if strpos(" `q_`qt'' ", " `tkr' "){
				replace weights`tkr' = 0 if qt_dt == `qt'
			}
			if strpos(" `q_fullsamp_`qt'' ", " `tkr' "){
				if !inlist("`tag'", "other"){
					replace total_`tag'_full = total_`tag'_full - full_weights`tkr' if qt_dt == `qt'
				}
				replace full_weights`tkr' = 0 if qt_dt == `qt'
			}
		}
}
egen assets_kmv_`tag' 	= rowtotal($`tag'_vars)
egen assets_kmv_`tag'_full 	= rowtotal($`tag'_vars_full)
}

/* egen total_compustat = rowtotal($dealer_vars) */
/* foreach var of varlist $dealer_vars{ */
/* //replace `var' = . if `var' == 0 */
/* gen test_perc_`var' = `var'/1000/total_dealer */
/* } */
/* gen cov_initial = total_compustat/1000/total_dealer */

save ../temp/firm_weightings.dta, replace

*And now 'weights_[tkr]' becomes a true weight

foreach tag in `tags'{
	foreach var in $`tag'_kmv{
		replace weights`var' = weights`var'/assets_kmv_`tag'
		replace full_weights`var' = full_weights`var'/assets_kmv_`tag'_full
}
}

*This is used in sum_stats_wholesample to show weights on some firms on the snapshot date
levelsof qt_dt if inrange(qt_dt, $charts_start, $charts_end) , local(quarters) clean
capture rm ../output/firm_list_final_wweights.dta
foreach quarter in `quarters'{

preserve
	disp %tq `quarter'

	keep if qt_dt == `quarter'
	keep weights* full*
	xpose, clear varname
	gen mkmv_id = regexs(1) if regexm(_varname, "weights(.*)")
	gen type = regexs(1) if regexm(_varname, "(.*weight)")
	drop _varname
	reshape wide v1, i(mkmv_id) j(type) string
	rename v1full_weight full_weight
	rename v1weight weight
	merge 1:m mkmv_id using ../temp/firm_list.dta
	keep if _merge == 3
	append using ../temp/broker_nodes_list_config$bd_config
	replace mkmv_id = mkmv_id + "_X" if inlist(sector, "Top 10 Dealers", "Top 11-25 Dealers")
	keep if qt_dt == `quarter'
	gen flag = 0
	replace flag = 1 if weight > 0
	bys sector: egen count = sum(flag)
	drop _merge
	tempfile intermed
	save `intermed', replace
	use ../temp/Contagion_preprocessing.dta, clear
	keep if qt_dt == `quarter'
	duplicates drop mkmv_id, force
	keep mkmv_id nm_short BHCK2170
	merge 1:1 mkmv_id using `intermed'
	keep if _merge == 3 | _merge == 2


	keep if qt_dt == `quarter'
	capture confirm file ../output/firm_list_final_wweights.dta
	if _rc != 0{
		save ../output/firm_list_final_wweights.dta, replace
	}
	else{
		tempfile temp
		save `temp'
		use ../output/firm_list_final_wweights.dta, clear
		append using `temp'
		save ../output/firm_list_final_wweights.dta, replace
	}


	*save ../output/firm_list_final_wweights.dta, replace

restore

}

drop edf01_mean*
*Each firm's contribution to the weighted average probability of its sector
foreach var of varlist edf*{
	local tkr = regexr("`var'", "edf01(_8)*", "")
	gen `var'_full = `var' * full_weights`tkr'
	replace `var' = `var' * weights`tkr'
}


gen ffunds_assets_in_samples = 0
foreach tag in `tags'{
	local `tag'_probs
	local `tag'_probs8
	local `tag'_probs_full
	foreach tkr in $`tag'_kmv{
		local `tag'_probs ``tag'_probs' edf01`tkr'
		local `tag'_probs8 ``tag'_probs8' edf01_8`tkr'
		local `tag'_probs_full ``tag'_probs_full'  edf01`tkr'_full
	}
	egen delta_`tag'		= rowtotal(``tag'_probs')
	egen delta_`tag'_full	= rowtotal(``tag'_probs_full')
	egen delta_`tag'8		= rowtotal(``tag'_probs8')

	gen coverage_`tag' 		= assets_kmv_`tag'/1000/total_`tag'

	replace ffunds_assets_in_samples = ffunds_assets_in_samples + assets_kmv_`tag'

}


format qt_dt %tq


foreach var of varlist total* dealer* insurance* reit* other* liab_total*{
	replace `var' = `var'*1000
}


*Output
drop edf* weights* full_weights*
save ../temp/deltas_external.dta, replace
