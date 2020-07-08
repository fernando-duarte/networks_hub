*********************************************************************************
*
*	Network output series construcion
*	Collin Jones and Rui Yu
*	Purpose: Construct variables needed to recreate fields in Glasserman(2014). Subsequent do-file
*		performs actual plotting. 
*
*********************************************************************************

clear

******************************************************************************************
*
*			Initial imports and definitions
*
******************************************************************************************

//First load and format Y15, so we don't need to load Y9C twice
import excel ../input/Y15.xlsx, firstrow
gen qt_dt 	= qofd(dofc(date))
rename ID_RSSD entity
save ../temp/Y15_formatted, replace

// Input prepared data

use ../temp/Contagion_preprocessing.dta, clear
drop _merge

merge 1:1 qt_dt entity using ../temp/Y15_formatted
sort qt_dt
drop if _merge == 2

gen liab_in_actual = liab_in + (1/2)*liab_in_unc

gen liab_in_possible = liab_in + liab_in_unc

*Generate variable to flag institutions that ever appear in the FR-Y15
levelsof mkmv_id if !missing(RISKM370), local(y15_sample)
global y15_sample = subinstr(`"`y15_sample'"', " ", ", ", 100)
gen y15_sample							= 0
replace y15_sample						= 1 if !missing(RISKM370)
bys mkmv_id: egen y15_sample_temp		= max(y15_sample)
replace y15_sample						= y15_sample_temp
drop y15_sample_temp

*Finding % of uninsured deposits (total - foreign + domestic) that Fr-Y15 says are inside system,
*	then extending that all back to the rest of the sample
gen y_15_intra_deposits = RISKM363	+	RISKM364
replace y_15_intra_deposits = . if tkr == "ALLY"
gen y9c_all_deposits_minusFDIC = BHCB2210 + BHCB3187 + BHCB2389 + BHCB6648 + BHCB2604 + BHOD3189 + BHOD3187 + BHOD2389 + BHOD2604 + BHOD6648 + BHFN6631 + BHFN6636 - tot_insur_depos
gen y9c_all_deposits_minusFDICfor = BHCB2210 + BHCB3187 + BHCB2389 + BHCB6648 + BHCB2604 + BHOD3189 + BHOD3187 + BHOD2389 + BHOD2604 + BHOD6648 - tot_insur_depos
gen y9c_all_deposits = BHCB2210 + BHCB3187 + BHCB2389 + BHCB6648 + BHCB2604 + BHOD3189 + BHOD3187 + BHOD2389 + BHOD2604 + BHOD6648 + BHFN6631 + BHFN6636
gen liab_in_y15deposits = y_15_intra_deposits +  BHDMB993 + BHCKB995 + BHCK3548 + (1/2)*(BHCK2309+BHCK2332+BHCK2333+BHCKC699+BHCKB557+BHCKB984)
gen perc_deposits_in_y15 = y_15_intra_deposits / y9c_all_deposits_minusFDIC
bys entity: egen avg_deposits_in = mean(perc_deposits_in_y15)
gen deposits_in_extrapolated = avg_deposits_in*(y9c_all_deposits_minusFDICfor)

*Regenerating some stats from Analysis_Y9C, taking this all into account
egen liab_in_y15 		= rowtotal(deposits_in_extrapolated BHDMB993 BHCKB995 BHCK3548)
egen liab_in_unc_y15 	= rowtotal(BHCK2309 BHCK2332 BHCK2333 BHCKC699 BHCKB557 BHCKB984)

*BK and STT likely have almost 100% of foreign deposits in the system
replace RISKM370 = RISKM370 - (BHFN6631 + BHFN6636) if inlist(mkmv_id, "857473", "064057")

/***************************************************************************
 Y15 off-balance sheet stuff
 
RISKJ458 = Assets, unused commited lines to other financial institutions
RISKM360 = Future exposure of derivatives

RISKM365 = liabilities, unused portion of commited lines from financial instutitions
RISKM368 = liabilities, potential future derivative exposure
****************************************************************************/

*Figure out what percentage of assets or liabilities tend to be made up of these off-balance-sheet things
gen perc_assets_offbalance = (RISKJ458 + RISKM360)/BHCK2170
gen perc_liab_offbalance = (RISKM365 + RISKM368)/BHCK2948
bys entity: egen avg_asset_perc_offbalance = mean(perc_assets_offbalance)
bys entity: egen avg_liab_perc_offbalance = mean(perc_liab_offbalance)
gen assets_offbalance_extrap = avg_asset_perc_offbalance * BHCK2170
gen liab_offbalance_extrap = avg_liab_perc_offbalance * BHCK2948

*Generating new beta variables using Y-15 data. Either using Y15 variable values themselves, or
*	the aggregate value of total intrafinancial liabilities
*	OR extrapolating the value of some off-balance sheet items back through time
gen y_15_in_compat 		= RISKM370 - RISKM368 - RISKM365
gen beta_y15 			= RISKM370/BHCK2948
gen beta_y15_compat 	= y_15_in_compat/BHCK2948
gen beta_y15_deposits 	= (liab_in_y15 + (1/2)*liab_in_unc_y15)/BHCK2948 if y15_sample == 1
gen beta_offbalance = ((liab_in + liab_offbalance_extrap) + (1/2)*liab_in_unc)/(BHCK2948 + liab_offbalance_extrap) if y15_sample == 1

*This line used to be in Analysis_Y9C. Now here instead, in case we want to hold onto Foreign BHCs. Uncomment to remove them. Note that most of them
*	don't have KMV probabilities, so this will sort of happen later.
*IMPORTANT: Dropping anything without a KMV default probability. 
drop if FGN_CALL_FAM_ID != 0 & entity != 3232316 & entity != 1132449
keep if !missing(prob_default_phys)
drop if no_call == 1
drop if tkr == "IFIN"

*Check for insured deposit double counting

replace id_rssd_top = entity if missing(id_rssd_top)
bys qt_dt id_rssd_top: gen count_top_matches = _N if matchedon != 3
gen double_count = 0
replace double_count = 1 if count_top_matches > 1 & matchedon == 2

******************************************************************************************
*
*			Dealer Aggregate Nodes
*
******************************************************************************************

*Format then merge list of nodes that appear in Top 25 broker dealers at each quarter
*	If they appear in those aggregated nodes, drop them from FR-Y9C sample
preserve
use ../temp/broker_nodes_list_config$bd_config, clear
keep qt_dt mkmv_id
gen broker_node = 1
save ../temp/broker_nodes_list_config$bd_config, replace
restore
drop _merge


/* if $output_simulation_data == 0{ */
*Comment out this next portion if we don't want the aggregate nodes
merge 1:1 qt_dt mkmv_id using ../temp/broker_nodes_list_config$bd_config
drop if _merge == 2
drop if broker_node == 1

*Append data for dealer aggregate nodes, from aggregated FOCUS SIFMA data
append using ../temp/broker_nodes_config$bd_config

******************************************************************************************
*
*			Some unit fixing, formatting, and initial identifiers
*
******************************************************************************************

gen keep_obs = 0
format qt_dt  %tq
drop rank

*This will give us a tag of the bank in the global macro above (which is the old sample). Helpful for setting beta_max later. 
foreach bank of numlist $bankSample {
	replace keep_obs = 1 if entity == `bank'
}
*Allow the new broker aggregated nodes to become the Beta+ for the NVI
replace keep_obs = 1 if tkr == "BRO10" | tkr == "BRO25"

drop _*
bys entity: egen count_q 	= count(qt_dt) if qt_dt <= $charts_end & qt_dt >= $charts_start
gen full_sample 			= 0
replace full_sample			= 1 if count_q == $charts_end - $charts_start + 1
gen crisis					= 0
replace crisis 				= 1 if inrange(qt_dt, yq(2008, 1), yq(2012, 2))

*Set time variable
xtset entity qt_dt, quarterly

******************************************************************************************
*
*			Defining Glasserman (and related) Fields
*
******************************************************************************************

*Beta = connectivity
gen beta 	= liabs_in_frac
gen beta_nocall = liabs_in_nocall_frac
replace beta_nocall = beta if missing(beta_nocall)

*w=net worth
gen w = BHCK2170*10^(-3) - BHCK2948*10^(-3)

* Assume that negative net worth is an already-defaulting firm that the bankrutpcy screens
* didn't pick up (seems to mostly be the case)
drop if w < 0

*c = outside assets
gen c 				= (BHCK2170*10^(-3) - (asset_in*10^(-3) + (1/2)*asset_in_unc*10^(-3))) // Convert to millions
gen c_y15 			= BHCK2170*10^(-3) - RISKM362*10^(-3)
gen c_offbalance 	= BHCK2170*10^(-3) + assets_offbalance_extrap*10^(-3)

gen assets 	= BHCK2170*10^(-3)
capture confirm variable assets_total
if _rc != 0{
	gen assets_total = .
}
replace assets_total = assets_total*10^(-3)
replace assets_total = assets if missing(assets_total)

*Ranking of firms
gsort qt_dt -assets_total
by qt_dt: gen rank_nodeal = _n
gen sim_nodes = 0
replace sim_nodes = 1 if rank_nodeal <= 10

*lambda = leverage of outside assets
gen lambda 	= c/w

*Probability of default from KMV (rename to Glasserman nomenclature)
rename prob_default_phys delta
rename prob_default_phys_alt delta_alt
rename prob_default_neut delta_neut
rename prob_default_phys8 delta8
rename prob_default_neut8 delta_neut8

*Contribution to weighted default probability
gen delta_c 				= delta*c
//gen delta_c_deal 				= delta*c
gen delta_c_y5				= delta *c_y15
gen delta_c_offbalance		= delta *c_offbalance
gen delta8_c				= delta8*c
gen delta_neut_c			= delta_neut*c

*Construct some helper variables
bys qt_dt: egen sum_c 					= sum(c)
bys qt_dt: egen sum_c_bhc 				= sum(c) if !inlist(tkr, "BRO10", "BRO25")
bys qt_dt: egen sum_c_top10				= sum(c) if sim_nodes == 1
bys qt_dt: egen sum_c_dealer 			= sum(c) if inlist(tkr, "BRO10", "BRO25")
bys qt_dt: egen sum_c_full				= sum(c) if full_sample == 1
bys qt_dt: egen sum_c_y15				= sum(c_y15) if y15_sample == 1
bys qt_dt: egen sum_c_y15_samp			= sum(c) if y15_sample == 1
bys qt_dt: egen sum_c_offbalance		= sum(c_offbalance) if y15_sample == 1
gen l_sum_c 							= ln(sum_c)
bys qt_dt: egen sum_w 					= sum(w)
bys qt_dt: egen sum_assets 				= sum(assets)
bys qt_dt: egen sum_assets_prob_covered	= sum(assets_total)
bys qt_dt: egen sum_assets_bhc 			= sum(assets) if !inlist(tkr, "BRO10", "BRO25")
bys qt_dt: egen sum_assets_full			= sum(assets) if full_sample == 1
bys qt_dt: egen sum_assets_y15			= sum(assets) if y15_sample == 1
bys qt_dt: egen sum_delta_c 			= sum(delta_c)
bys qt_dt: egen sum_delta_c_bhc			= sum(delta_c) if !inlist(tkr, "BRO10", "BRO25")
bys qt_dt: egen sum_delta_c_top10		= sum(delta_c) if sim_nodes == 1
bys qt_dt: egen sum_delta_c_dealer		= sum(delta_c) if inlist(tkr, "BRO10", "BRO25")
bys qt_dt: egen sum_delta_c_full 		= sum(delta_c) if full_sample == 1
bys qt_dt: egen sum_delta_c_y15 		= sum(delta_c_y5) if y15_sample == 1
bys qt_dt: egen sum_delta_c_offbalance 	= sum(delta_c_offbalance) if y15_sample == 1
bys qt_dt: egen sum_delta_c_y15_samp 	= sum(delta_c) if y15_sample == 1
bys qt_dt: egen sum_delta8_c 			= sum(delta8_c)
bys qt_dt: egen sum_delta_neut_c 		= sum(delta_neut_c)

******************************************************************************************
*
*			Defining Beta+ and second, third, fourth, etc Betas
*						
******************************************************************************************

*Produce ranked fractions of finacial liabilities. Benchmark only takes Beta+ from top BHCs.
*Currently goes out to fifth highest
foreach v in max_bank max_bank_all max_bank_fullsample max_bank_y15 max_bank_y15_compat max_bank_y15depos max_bank_top10 ///
	max_bank_offbalance max_bank_bhc max_bank_nocall max_bank_y15_samp{
	gen `v' = .
}

gsort qt_dt -BHCK2170 
by qt_dt: gen rank = _n
drop if tot_insur_depos < 0
foreach num of numlist 1/5 {
	bys qt_dt: egen beta_max_`num'	= max(beta) if max_bank == . & keep_obs == 1
	replace max_bank = `num' if beta_max_`num'== beta
	
	bys qt_dt: egen beta_max_bhc_`num'	= max(beta) if max_bank_bhc == . & (keep_obs == 1 | !inlist(tkr, "BRO10", "BRO25")) 
	replace max_bank_bhc = `num' if beta_max_bhc_`num'== beta & missing(max_bank_bhc)
	
	bys qt_dt: egen beta_max_top10_`num'	= max(beta) if max_bank_top10 == . & (sim_nodes == 1) 
	replace max_bank_top10= `num' if beta_max_top10_`num'== beta & missing(max_bank_top10)
	
	bys qt_dt: egen beta_max_y15_depos_`num' = max(beta_y15_deposits) if max_bank_y15depos == . & y15_sample == 1 & keep_obs == 1
	replace max_bank_y15depos = `num' if beta_max_y15_depos_`num'== beta_y15_deposits
	
	bys qt_dt: egen beta_max_offbalance_`num' = max(beta_offbalance) if max_bank_offbalance == . & y15_sample == 1 & keep_obs == 1
	replace max_bank_offbalance = `num' if beta_max_offbalance_`num'== beta_offbalance

	bys qt_dt: egen beta_max_all_`num'	= max(beta) if max_bank_all == .
	replace max_bank_all = `num' if beta_max_all_`num'== beta
	
	bys qt_dt: egen beta_max_all_nocall_`num'	= max(beta_nocall) if max_bank_nocall == .
	replace max_bank_nocall = `num' if beta_max_all_nocall_`num'== beta_nocall
	
	bys qt_dt: egen beta_max_fullsample_`num' = max(beta) if max_bank_fullsample == . & keep_obs == 1 & full_sample == 1
	replace max_bank_fullsample = `num' if beta_max_fullsample_`num'== beta

	bys qt_dt: egen beta_max_y15_`num' = max(beta_y15) if max_bank_y15 == . & keep_obs == 1 & y15_sample == 1
	replace max_bank_y15 = `num' if beta_max_y15_`num'== beta_y15
	
	bys qt_dt: egen beta_max_y15_samp_`num' = max(beta) if max_bank_y15_samp == . & keep_obs == 1 & y15_sample == 1
	replace max_bank_y15_samp = `num' if beta_max_y15_samp_`num'== beta

	bys qt_dt: egen beta_max_y15_compat_`num' = max(beta_y15_compat) if max_bank_y15_compat == . & keep_obs == 1 & y15_sample == 1
	replace max_bank_y15_compat = `num' if beta_max_y15_compat_`num'== beta_y15_compat

}


*Create contagion index
gen contag_index_glasserman = beta * (c - w) //Glasserman-Young Contagion index = w*beta*(lambda-1)
gen contag_index = beta*c //Duarte Jones Contagion index = w*beta*lambda 

gen gamma_max = (1/beta_max_1)-1

*Run sum stats now, before things get more complicated 
/* if $output_simulation_data == 0{ */
preserve
drop if inlist(tkr, "BRO10", "BRO25")
do ../code/sum_stats_wholesample.do
restore
/* } */

***************************************************************************************************
*
*			Merging in Flow of Funds work, building new indices with that data...
*
***************************************************************************************************

merge m:1 qt_dt using ../temp/deltas_external.dta, nogen

local denom sum_c
local num sum_delta_c
*Create some more relevant fields for the subsectors being included in aggregate. This is no longer
*	necessary for dealers, as those are included as their own rows.
foreach tag in insurance reit other{
	egen temp_`tag'_in = rowtotal(`tag'_in*)
	egen temp_`tag'_unc = rowtotal(`tag'_unc*)
	gen perc_in_`tag' = (temp_`tag'_in + 0.5 * temp_`tag'_unc)/total_`tag'_orig
	local num `num' + (1-perc_in_`tag')*total_`tag'*delta_`tag'
	local denom `denom' + (1-perc_in_`tag')*total_`tag'
	gen perc_in_`tag'_temp = .
	
	local assets_to_include `assets_to_include' + assets_kmv_`tag'
	local assets_possible `assets_possible' + total_`tag'
	local assets_possible_orig `assets_possible_orig' + total_`tag'_orig
}
foreach tag in insurance reit other{
	tempvar contrib
	gen `contrib' = ((1-perc_in_`tag')*total_`tag'*delta_`tag') / ((1-(1+$gamma_benchmark/100)*beta_max_1)*(`denom'))
	bys qt_dt: egen contribution_`tag' = max(`contrib')
	gen weight_`tag' = ((1-perc_in_`tag')*total_`tag')/(`denom')
	
	replace `contrib' = ((1-perc_in_`tag')*(assets_kmv_`tag')*delta_`tag') / ((1-(1+$gamma_benchmark/100)*beta_max_1)*(`denom'))
	bys qt_dt: egen contribution_kmv_`tag' = max(`contrib')
	
	replace `contrib' = ((1-perc_in_`tag')*(total_`tag' - assets_kmv_`tag')*delta_`tag') / ((1-(1+$gamma_benchmark/100)*beta_max_1)*(`denom'))
	bys qt_dt: egen contribution_nonkmv_`tag' = max(`contrib')
}
tempvar temp_bhc temp_dealer temp_bhc_wt temp_dealer_wt
gen contribution_bhc 		= sum_delta_c_bhc/((1-(1+$gamma_benchmark/100)*beta_max_1)*(`denom'))
gen contribution_dealer		= sum_delta_c_dealer/((1-(1+$gamma_benchmark/100)*beta_max_1)*(`denom'))
gen weight_bhc = sum_c_bhc/(`denom')
gen weight_dealer = sum_c_dealer/(`denom')
bys qt_dt: egen `temp_bhc' 	= max(contribution_bhc)
bys qt_dt: egen `temp_dealer' 	= max(contribution_dealer)
bys qt_dt: egen `temp_bhc_wt' 	= max(weight_bhc)
bys qt_dt: egen `temp_dealer_wt'	= max(weight_dealer)
replace contribution_bhc = `temp_bhc'
replace contribution_dealer = `temp_dealer'
replace weight_bhc = `temp_bhc_wt'
replace weight_dealer = `temp_dealer_wt'

* Indices with differnet subsamples included. Note that, in this iteration, dealers are included using ROWS for Top 10 and Top 11-25. All others are included in aggregate from Flow of funds, with
*	columns (constant within quarters) that give the necessary fields
gen index_bhc_only = 100*(sum_delta_c_bhc)/((1-(1+$gamma_benchmark/100 )*beta_max_bhc_1)*(sum_c_bhc))
gen index_top10_only = 100*(sum_delta_c_top10)/((1-(1+$gamma_benchmark/100 )*beta_max_top10_1)*(sum_c_top10))
gen index_w_insurance = 100*(sum_delta_c + (1-perc_in_insurance)*total_insurance*delta_insurance)/((1-(1+$gamma_benchmark/100)*beta_max_1)*(sum_c+(1-perc_in_insurance)*total_insurance))
gen index_w_dealer = 100*(sum_delta_c)/((1-(1+$gamma_benchmark/100)*beta_max_1)*(sum_c))
gen index_w_other = 100*(sum_delta_c + (1-perc_in_other)*total_other*delta_other)/((1-(1+$gamma_benchmark/100)*beta_max_1)*(sum_c + (1-perc_in_other)*total_other))
gen index_w_reit = 100*(sum_delta_c + (1-perc_in_reit)*total_reit*delta_reit)/((1-(1+$gamma_benchmark/100)*beta_max_1)*(sum_c + (1-perc_in_reit)*total_reit))

gen nvi_benchmark = 100*(`num')/((1-(1+$gamma_benchmark/100)*beta_max_1)*(`denom'))

*Breaking down benchmark into multiplicative components
gen agg_loss = 100*(`num')/(`denom')
gen sum_delta_c_all = (100*(`num'))/1000000
gen sum_c_all = (`denom')/1000000
gen connectivity_comp = 1/(1-(1+$gamma_benchmark/100)*beta_max_1)

******************************************************************************************
*
*			Some misc stats on the Y9C sample (need Flow of Funds data to calc)
*						
******************************************************************************************

*Coverage series for Y9C, BHC sample
gen bhc_ffunds = total_us_depository + total_holding_companies
gen coverage_bhc = sum_assets_bhc / bhc_ffunds if max_bank_bhc ==1
format qt_dt %tqCCYY

tempvar focus
gen `focus' = assets if inlist(tkr, "BRO10", "BRO25")
bys qt_dt: egen dealer_FOCUS = total(`focus')

gen total_assets_in_all = sum_assets_prob_covered `assets_to_include'
gen total_assets_cov_all = sum_assets `assets_possible'
gen total_assets_sector = bhc_ffunds + dealer_FOCUS `assets_possible'

gen coverage_included = total_assets_in_all/total_financial_ours
gen coverage_possible = total_assets_cov_all/total_financial_ours
gen coverage_sectors = total_assets_sector/total_financial_ours

*"Average Delta" series for BHCs
bys qt_dt: gen weight_assets = assets/sum_assets_bhc if !inlist(tkr, "BRO10", "BRO25")
bys qt_dt: gen weight_assets_full = assets/sum_assets_full if full_sample == 1
//bys qt_dt: gen weight_outside_assets = (1-asset_in_frac)*assets/sum_c
gen prob_contribution_assets = delta*weight_assets
gen prob_contribution_assets_full = delta*weight_assets_full
gen prob_contribution_assets8 = delta8*weight_assets
//gen prob_contribution_outside_assets = delta*weight_outside_assets
bys qt_dt: egen avg_delta_assets = sum(prob_contribution_assets)
//bys qt_dt: egen avg_delta_assets_full = sum(prob_contribution_assets_full) if full_sample == 1
bys qt_dt: egen avg_delta8_assets = sum(prob_contribution_assets8)
//bys qt_dt: egen avg_delta_outside_assets = sum(prob_contribution_outside_assets)


***************************************************************************************************
*
*			A number of different index configurations, calculated for robustness purposes
*
***************************************************************************************************

*Indices with a range of %assets inside or outside the system (for broker dealer and insurance companies)
local per_out_grades

*Benchmark index with second, third, etc highest betas
foreach val of numlist 2/5{
	gen index_w_ff_beta`val' = 100*(`num')/((1-(1+$gamma_benchmark/100)*beta_max_`val')*(`denom'))
}

*Benchmark index selecting beta+ (and second, third, etc) from whole Y9C sample
foreach val of numlist 1/5{
	gen index_w_ff_beta_all`val'= 100*(`num')/((1-(1+$gamma_benchmark/100)*beta_max_all_`val')*(`denom'))
	
	gen index_w_ff_beta_full`val'= 100*(`num')/((1-(1+$gamma_benchmark/100)*beta_max_fullsample_`val')*(`denom'))
}

*Index with bankruptcy costs. Global defined in NtwkCntgn.do
foreach val in $gammas 0{
	gen index_w_ff_gam`val' = 100*(`num')/((1-(1+`val'/100)*beta_max_1)*(`denom'))
	replace index_w_ff_gam`val' = . if `val'/100 >= gamma_max

}

*Index with some different Y9C "unclear" classifications
	*Note these are percentages IN the system
forval perc=10(10)100{
	
	//replace perc_in_dealer_temp = (temp_dealer_in + (`perc'/100) * temp_dealer_unc)/total_dealer_orig
	replace perc_in_insurance_temp = (temp_insurance_in + (`perc'/100) * temp_insurance_unc)/total_insurance_orig
	replace perc_in_reit_temp = (temp_reit_in + (`perc'/100) * temp_reit_unc)/total_reit_orig
	replace perc_in_other_temp = (temp_other_in + (`perc'/100) * temp_other_unc)/total_other_orig

	gen liab_in_frac_temp`perc' =  (liab_in + (`perc'/100)*liab_in_unc)/BHCK2948
	gen c_temp				= (BHCK2170*10^(-3) - (asset_in*10^(-3) + (1 - `perc'/100)*asset_in_unc*10^(-3)))
	bys qt_dt: egen beta_max_1_perc`perc' = max(liab_in_frac_temp`perc') if keep_obs == 1
	bys qt_dt: egen beta_max_1_perc`perc'_actual = max(beta_max_1_perc`perc')

	gen delta_c_temp					= delta*c_temp
	bys qt_dt: egen sum_delta_c_temp 	= sum(delta_c_temp)
	bys qt_dt: egen sum_c_temp 			= sum(c_temp)
	
	gen index_wffunds_unc_`perc'perc = 100*(sum_delta_c_temp + (1-perc_in_insurance_temp)*total_insurance*delta_insurance + ///
		(1-perc_in_reit_temp)*total_reit*delta_reit+ (1-perc_in_other_temp)*total_other*delta_other)/((1-(1+$gamma_benchmark/100)*beta_max_1_perc`perc'_actual)*(sum_c_temp+(1-perc_in_insurance_temp)*total_insurance + ///
		(1-perc_in_reit_temp)*total_reit+ (1-perc_in_other_temp)*total_other))
	
	drop delta_c_temp sum_delta_c_temp sum_c_temp c_temp
}

*When you do this, you also have to re-calculated betas. This is done below.
gen beta_temp = .
forval perc=10(10)100{
	replace beta_temp = ((1/2)*liab_in_unc + (`perc')/100 * tot_uninsur_depos_y9c + BHDMB993 + BHCKB995 + BHCK3548)/BHCK2948
	bys qt_dt: egen beta_max_temp = max(beta_temp) if keep_obs == 1
	
	gen index_wffunds_unc_depos_`perc'in = 100*(`num')/((1-(1+$gamma_benchmark/100)*beta_max_temp)*(`denom'))

	drop beta_max_temp
}

gen nvi_benchmark_all_nocall = 100*(`num')/((1-(1+$gamma_benchmark/100)*beta_max_all_nocall_1)*(`denom'))
gen nvi_benchmark_full = 100*(sum_delta_c_full + (1-perc_in_other)*total_other_full*delta_other_full+ (1-perc_in_insurance)*total_insurance_full*delta_insurance_full + (1-perc_in_reit)*total_reit_full*delta_reit_full)/((1-(1+$gamma_benchmark/100)*beta_max_fullsample_1)*(sum_c_full+(1-perc_in_insurance)*total_insurance_full + (1-perc_in_reit)*total_reit_full + (1-perc_in_other)*total_other_full))
gen nvi_benchmark8 = 100*(sum_delta8_c + (1-perc_in_other)*total_other*delta_other8 + (1-perc_in_insurance)*total_insurance*delta_insurance8 + (1-perc_in_reit)*total_reit*delta_reit8)/((1-(1+$gamma_benchmark/100)*beta_max_1)*(sum_c+(1-perc_in_insurance)*total_insurance + (1-perc_in_reit)*total_reit + (1-perc_in_other)*total_other))
gen nvi_benchmark_deltafix = 100*$delta_fixed*(sum_c + (1-perc_in_other)*total_other+ (1-perc_in_insurance)*total_insurance + (1-perc_in_reit)*total_reit)/((1-(1+$gamma_benchmark/100)*beta_max_1)*(sum_c+(1-perc_in_insurance)*total_insurance + (1-perc_in_reit)*total_reit + (1-perc_in_other)*total_other))
gen nvi_benchmark_y15cov = 100*(`num')/((1-(1+$gamma_benchmark/100)*beta_max_y15_1)*(`denom'))
gen nvi_benchmark_y15depos = 100*(`num')/((1-(1+$gamma_benchmark/100)*beta_max_y15_depos_1)*(`denom'))
gen nvi_y15_only = 100*(sum_delta_c_y15_samp)/((1-(1+$gamma_benchmark/100)*beta_max_y15_samp_1)*(sum_c_y15_samp))
gen nvi_offbalance = 100*(sum_delta_c_offbalance)/((1-(1+$gamma_benchmark/100)*beta_max_offbalance_1)*(sum_c_offbalance))
gen nvi_y15_only_wconnect = 100*(sum_delta_c_y15)/((1-(1+$gamma_benchmark/100)*beta_max_y15_1)*(sum_c_y15))

***************************************************************************************************
*
*					Y15 Exploratory Work
*
***************************************************************************************************

format qt_dt %tqCCYY

gen liab_in_depos15 = deposits_in_extrapolated + BHDMB993 + BHCKB995 + BHCK3548 + BHCK4062 + BHCK3049 
gen beta_depos15 = (liab_in_depos15 + (1/2)*liab_in_unc)/BHCK2948


local perc_plot
foreach var of varlist y_15_intra_deposits BHDMB993 BHCKB995 BHCK3548 BHCK4062 BHCK3049{
	gen perc_`var'_y9c = 100*`var'/liab_in_y15deposits
	label var perc_`var'_y9c "`var'"
	local perc_plot `perc_plot' `var'
}
foreach var of varlist BHCK2309 BHCK2332 BHCK2333 BHCKC699 BHCKB557 BHCKB984{
	gen perc_`var'_y9c = 100*(1/2)*`var'/liab_in_y15deposits
	label var perc_`var'_y9c "`var'"
	local perc_plot `perc_plot' `var'
}

local perc_plot_y15
foreach var of varlist y_15_intra_deposits RISKY833 RISKM365 RISKM366 RISKM367 RISKM368{
	gen perc_`var'_y15 = 100*`var'/RISKM370
	label var perc_`var'_y15 "`var'"
	local perc_plot_y15 `perc_plot_y15' `var'
}

/* if $output_simulation_data == 0{ */
	*Output large dataset. Will be used for Plots_Appendix, and general exploratory work.
	save ../output/Contagion_Data, replace


	*Output smaller dataset with only those variables need for Plots_Paper
	preserve

	keep qt_dt name nm_short tkr entity keep_obs beta w c assets lambda delta delta_neut delta_c ///
		sum_c sum_assets sum_delta_c max_bank beta_max_1 beta_max_all_1 beta_max_fullsample_1 ///
		max_bank_fullsample max_bank_all contag_index nvi_benchmark perc_in_insurance total_insurance ///
		delta_insurance total_dealer delta_dealer coverage_insurance max_bank_bhc ///
		coverage_dealer coverage_bhc BHCK2170 agg_loss connectivity_comp index_w_ff_gam1 ///
		index_w_ff_gam5 index_w_ff_gam10 index_w_ff_gam15 index_w_ff_gam30 index_wffunds_unc_30perc ///
		index_wffunds_unc_50perc index_wffunds_unc_70perc index_wffunds_unc_depos_20in index_wffunds_unc_depos_40in ///
		index_wffunds_unc_depos_60in index_wffunds_unc_depos_80in index_w_ff_beta_all1 index_w_ff_beta_full1 ///
		delta_insurance delta_dealer avg_delta_assets nvi_benchmark_full coverage_reit delta_reit coverage_included coverage_possible coverage_sectors ///
		coverage_other delta_other *_orig assets_kmv_* bhc_ffunds total_financial_ours nvi_y15_only nvi_benchmark_y15cov nvi_y15_only_wconnect max_bank_y15 ///
		index_wffunds_unc_depos_20in nvi_benchmark_y15depos max_bank_y15depos nvi_offbalance max_bank_offbalance index_w_ff_beta2 index_w_ff_beta3 index_w_ff_beta_all1 ///
		beta_max_2 beta_max_3 max_bank_y15_samp contribution* index_top10_only beta_max_top10_1 max_bank_top10 weight_*
	save ../output/Contagion_Data_select, replace

	restore
/* } */
***************************************************************************************************
*
*			Putting into excel file for simulations
*
***************************************************************************************************

/* if $output_simulation_data == 1{ */
preserve
format qt_dt %8.0g
tempvar temp
bys qt_dt: egen `temp' = max(nvi_benchmark)
replace nvi_benchmark = `temp'
drop if inlist(tkr, "BRO10", "BRO25")
gen p_bar = BHCK2948*(10^(-3))
gen b = (BHCK2948 - (liab_in + (1/2)*liab_in_unc))*(10^(-3))

/* keep nm_short tkr p_bar assets c b delta delta_alt beta w qt_dt nvi_benchmark */
foreach var in $simulation_data_vars{
	capture confirm variable `var', exact
	disp _rc
	if _rc != 0{
		gen `var' = .
		}
}
export excel $simulation_data_vars using ../temp/node_stats_forsimulation_all.xls, sheet("BHCs", replace) firstrow(variables)
restore

gen p_bar = .
gen b = .
local insurance_name "Insurance Aggregate"
local reit_name "REIT Aggregate"
local other_name "Other Aggregate"

local insurance_tkr "INSUR"
local reit_tkr "REIT"
local other_tkr "OTHER"
tempfile temp

levelsof qt_dt, local(qt) clean
foreach qt in `qt'{
preserve
keep if qt_dt == `qt'
keep if _n == 1
expand 4 if _n == 1
local row  2
foreach sec in insurance reit other{
	replace nm_short = "``sec'_name'" if _n == `row'
	replace tkr = "``sec'_tkr'" if _n == `row'
	replace c = (1-perc_in_`sec')*total_`sec' if _n == `row'
	replace p_bar = liab_total_`sec' if _n == `row'
	replace assets = total_`sec' if _n == `row'
	replace delta = delta_`sec' if _n == `row'
	replace w = total_`sec' - liab_total_`sec' if _n == `row'
	replace beta = . if _n == `row'
	replace b = . if _n == `row'
	local row = `row' + 1	
}
	drop if _n == 1
	/* keep qt_dt nm_short tkr c p_bar assets delta w beta b */
	capture append using `temp'
	save `temp', replace
restore
}
use `temp', clear
format qt_dt %8.0g
foreach var in $simulation_data_vars{
	capture confirm variable `var', exact
	disp _rc
	if _rc != 0{
		gen `var' = .
		}
}

export excel $simulation_data_vars using ../temp/node_stats_forsimulation_all.xls, sheet("Approx Aggregates", replace) firstrow(variables)
/* } */
