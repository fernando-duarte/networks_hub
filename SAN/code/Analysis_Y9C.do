*********************************************************************************
*
*	Y9C Analysis
*	Authors: Collin Jones and Rui Yu
*	Updated(8/9/18): Updated to pull data up to 2018Q1 - Francisco Ruela
*	Purpose: Execute variable categorizations defined in codebook to get numbers for each BHC's 
*		assets in and out of system. Also: clean and merge with KMV data.
*
*********************************************************************************

clear
set more off

*********************************************************************************
*
*	First, some KMV data cleaning.
*		Considering moving these to a separate file. 
*
*********************************************************************************


foreach ver in RawDFS RawDFS8{

local filename = regexr("`ver'", "Raw", "")
capture confirm file ../temp/`filename'.dta

if _rc != 0{
	use ../temp/`ver'.dta, clear
	*Date clean-up
	*Note: might need to use this first line if using KMV V8

	if "`ver'" == "RawDFS8"{
		gen date0 	= dofc(dt_day)
	}
	else{
	capture gen date0 	= date(dt_day, "YMD")
	if _rc != 0{
		gen date0 = dt_day
	}
	}

	merge m:1 mkmv_id using ../temp/bankruptcies.dta
	drop _merge

	foreach k in 1 2 3{
		drop if date0 >= bankruptcy_datetime`k' & !missing(bankruptcy_datetime`k') & bankruptcy_datetime`k' >= dofq($charts_start)
	}

	*Uncomment if dropping based on high EDF scores
	/*
	if "`ver'" == "RawDFS"{
	preserve
		keep if kmv_country_id == "USA" 
		keep if date0 >= $charts_start
		keep if inlist(mkmvind_name_desc, "BANKS AND S&LS", "FINANCE COMPANIES", "FINANCE NEC", "INSURANCE - LIFE", ///
			"INSURANCE - PROP/CAS/HEALTH", "INVESTMENT MANAGEMENT", "REAL ESTATE INVESTMENT TRUSTS", "SECURITY BROKERS & DEALERS")
		drop if mkmv_id == "N04946"
		keep if annualized_edf_01yr_pct >= 34
		save ../temp/high_edf_drops, replace
	restore
	}
	drop if annualized_edf_01yr_pct >= 34
	*/

	gen qt_dt 	= qofd(date0)
	gen year 	= year(date0)
		
	*Will aggregate both by last obs and by mean. The one named edf01 will be used in actual series.
	rename annualized_edf_01yr_pct edf01
	rename annualized_qdf_01yr_pct qdf01
	gen edf01_mean = edf01
	collapse (mean) edf01_mean (lastnm) edf01 qdf01 (lastnm) mkmvind_name_desc total_assets_current_amt total_liabilities_amt kmv_country_id, by(mkmv_id qt_dt)
	replace 				edf01 = edf01/100
	replace 				edf01_mean = edf01_mean/100
	replace 				qdf01 = qdf01/100

	bys mkmv_id (qt_dt): replace mkmvind_name_desc = mkmvind_name_desc[_n-1] if mkmvind_name_desc == ""
	bys mkmv_id (qt_dt): replace mkmvind_name_desc = mkmvind_name_desc[_N]

	if "`ver'" == "RawDFS8"{
		keep mkmv_id qt_dt edf01 qdf01 edf01_mean
		rename edf01 edf01_8
		rename qdf01 qdf01_8
		}
	else{
		*KMV industry category sometimes missing. We forward fill on the industry tag,
		*	then sets all industry tags within a company to the last value for that company.
		*	A problemmatic if a firm legitimately switches industries mid-sample. Haven't spotted
		*	any major issues, though.
		keep mkmv_id qt_dt edf01 qdf01 mkmvind_name_desc total_assets_current_amt total_liabilities_amt kmv_country_id edf01_mean
		bys mkmv_id (qt_dt): replace mkmvind_name_desc = mkmvind_name_desc[_n-1] if mkmvind_name_desc == ""
		bys mkmv_id (qt_dt): replace mkmvind_name_desc = mkmvind_name_desc[_N]
	}

	save ../temp/`filename'.dta, replace
}

}

*********************************************************************************
*
*					Preliminary Y9C clean-up 
*
*********************************************************************************

use ../temp/RawY9C.dta, clear
rename entity entity_type2
rename ID_RSSD entity

tostring dt, replace
gen year 	= substr(dt, 1, 4)	
gen month 	= substr(dt, 5, 2)
destring year, replace
destring month, replace

*Define quarters
gen quarter 	= 1 if month >= 1 & month <= 3
replace quarter = 2 if month >= 4 & month <= 6
replace quarter = 3 if month >= 7 & month <= 9
replace quarter = 4 if month >= 10 & month <= 12

*Subset to date range
gen qt_dt 	= yq(year, quarter)

keep if qt_dt >= $charts_start

*Use link from Match_RSSID_MKMVID to match Y9C BHCs to their KMVID
joinby entity using ../temp/RSSID_MKMVID.dta, unmatched(master)
drop _merge

*Merge in KMV default probs. Drop unmatched KMV rows (i.e. not in Y9C)
merge m:1 mkmv_id qt_dt using ../temp/DFS.dta
keep if _merge == 3 | _merge == 1
drop _merge
merge m:1 mkmv_id qt_dt using ../temp/DFS8.dta
keep if _merge == 3 | _merge == 1
drop _merge

*Average over edf01 duplicates.
bys entity qt_dt: egen prob_default_phys 	= mean(edf01)
bys entity qt_dt: egen prob_default_phys_alt 	= mean(edf01_mean)
bys entity qt_dt: egen prob_default_neut 	= mean(qdf01)
bys entity qt_dt: egen prob_default_phys8 	= mean(edf01_8)
bys entity qt_dt: egen prob_default_neut8 	= mean(qdf01_8)
duplicates drop entity qt_dt, force

*Drop if a BHC is only there for one period
bys entity: gen count 		= _n
bys entity: egen count_max 	= max(count)
drop if count_max == 1

drop count count_max

*Year/quarter convention. Q1 = YYYY00, Q2 = YYYY25, Q3 = YYYY50, Q4 = YYYY75
gen yyyyq 		= 100*year if quarter == 1
replace yyyyq 	= 100*year + 25 if quarter == 2
replace yyyyq 	= 100*year + 50 if quarter == 3
replace yyyyq 	= 100*year + 75 if quarter == 4

gsort qt_dt -BHCK2170
bys qt_dt: gen rank = _n


******************************************************************************************
*
*					Asset classifications
*						- Using codebook's classifications to put assets into "in" or "out"
*						- Our convention is to put 50% in, 50% out when we're not sure
*						- frac_[var] fields try to pull variables out from more aggregated 
*							fields in earliy Y9C formats. Usually the % of a variable over a year's time.
*
******************************************************************************************

foreach var of varlist BH*{
	capture confirm numeric variable `var'
	if !_rc{
		replace `var' = 0 if missing(`var')
	}
}	

*Component of in-financial assets unchanging from 3/2002
egen assets_in_unchanged 	= rowtotal(BHCK0397	BHDMB987 BHCKB989 BHCK2130)
egen assets_unc_unchanged 	= rowtotal(BHCK0395	BHCKA511 BHCK3163 BHCKA519 BHCK6438	BHCKB026 BHCK5507 BHCKB556 BHCKA520 BHCK2168)
	***BHCKA519, BHCKA520 end on 2018-03-31

*Non Agency Hold-To-Maturity MBS
*In
	*BHCKG320 decomposition
	gen frac_BHCKG320 					= BHCKG320/(BHCKG308+BHCKG320+BHCKG324+BHCKG328) if yyyyq >= 200925 & yyyyq < 201025
	replace frac_BHCKG320 				= 0 if frac_BHCKG320 == . & yyyyq >= 200925 & yyyyq < 201025
	bys entity: egen avg_frac_BHCKG320 	= mean(frac_BHCKG320) if yyyyq >= 200925 & yyyyq < 201025
	bys entity: egen temp 				= max(avg_frac_BHCKG320)
	replace avg_frac_BHCKG320 			= temp if avg_frac_BHCKG320 == .
	qui: sum frac_BHCKG320 if yyyyq >= 200925 & yyyyq < 201025
	replace avg_frac_BHCKG320 = `r(mean)' if missing(avg_frac_BHCKG320)
	drop temp frac_BHCKG320

	*define
	gen assets_htm_mbs 		= BHCKG320 if yyyyq >=200925
	replace assets_htm_mbs 	= (BHCK1709+BHCK1733)*avg_frac_BHCKG320 if yyyyq >=200200 & yyyyq < 200925

*Unclear
	gen assets_in_htm_mbs_unc 		= (BHCKG308+BHCKK146+BHCKK154) if yyyyq >=201100
	replace assets_in_htm_mbs_unc 	= BHCKG308+BHCKG324+BHCKG328 if yyyyq >= 200925 & yyyyq < 201100
	replace assets_in_htm_mbs_unc	= (1-avg_frac_BHCKG320)*(BHCK1709+BHCK1733) if yyyyq >= 200200 & yyyyq < 200925
	
*Non Agency Available-For-Sale MBS
*In
	*BHCKG323 decomposition	
	gen frac_BHCKG323 					= BHCKG323/(BHCKG311+BHCKG323+BHCKG327+BHCKG331) if yyyyq >= 200925 & yyyyq < 201025
	replace frac_BHCKG323 				= 0 if frac_BHCKG323 == . & yyyyq >= 200925 & yyyyq < 201025						
	bys entity: egen avg_frac_BHCKG323 	= mean(frac_BHCKG323) if yyyyq >= 200925 & yyyyq < 201025
	bys entity: egen temp 				= max(avg_frac_BHCKG323)
	replace avg_frac_BHCKG323 			= temp if avg_frac_BHCKG323 == .
	qui: sum frac_BHCKG323 if yyyyq >= 200925 & yyyyq < 201025
	replace avg_frac_BHCKG323 = `r(mean)' if missing(avg_frac_BHCKG323)
	drop temp frac_BHCKG323

	*Define 
	gen assets_in_afs_mbs			= BHCKG323 if yyyyq >=200925
	replace assets_in_afs_mbs 		= (BHCK1713+BHCK1736)*avg_frac_BHCKG323 if yyyyq >=200200 & yyyyq < 200925

*Unclear
	gen assets_in_afs_mbs_unc 		= (BHCKG311+BHCKK157+BHCKK149) if yyyyq >=201100
	replace assets_in_afs_mbs_unc 	= (BHCKG311+BHCKG327+BHCKG331) if yyyyq >= 200925 & yyyyq < 201100
	replace assets_in_afs_mbs_unc	= (1-avg_frac_BHCKG323)*(BHCK1713+BHCK1736) if yyyyq >= 200200 & yyyyq < 200925

*Non Agency MBS Trading
	gen assets_in_trad_mbs_unc  = BHCK3536 + BHCM3536 + BHCKG381 + BHCKG382 + BHCKK198
	***BHCK3536 and BHCM3536 end 2009Q1, BHCKG382 ends 2010Q4 - make sure they don't get reassigned

*Other Debt Securities
*In
	*BHCKG383+BHCKG384+BHCKG385 decomposition
	gen frac_oth_debt 					= (BHCKG383+BHCKG384+BHCKG385)/(BHCKG383+BHCKG384+BHCKG385+BHCKG386) if yyyyq >= 200925 & yyyyq < 201025
	replace frac_oth_debt 				= 0 if frac_oth_debt == . & yyyyq >= 200925 & yyyyq < 201025		
	bys entity: egen avg_frac_oth_debt 	= mean(frac_oth_debt) if yyyyq >= 200925 & yyyyq < 201025
	bys entity: egen temp 				= max(avg_frac_oth_debt)
	replace avg_frac_oth_debt 			= temp if avg_frac_oth_debt == .
	qui: sum frac_oth_debt if yyyyq >= 200925 & yyyyq < 201025
	replace avg_frac_oth_debt = `r(mean)' if missing(avg_frac_oth_debt)
	drop temp frac_oth_debt

	*Define
	gen assets_in_oth_debt			= BHCKG383+BHCKG384+BHCKG385 if yyyyq >=200925
		***ALL THREE OF THESE ASSETS, BHCKG383-385 end in 2018Q1
	replace assets_in_oth_debt 		= BHCM3537*avg_frac_oth_debt if yyyyq >=200800 & yyyyq < 200925
	replace assets_in_oth_debt 		= BHCK3537*avg_frac_oth_debt if yyyyq >=199500 & yyyyq < 200825
	*Note: G386 is out

*Non Agency Available-For-Sale CMO
*Unclear
	*BHCKG315 decomposition
	gen frac_afs_cmo 					= BHCKG319/(BHCKG315+BHCKG319) if yyyyq >= 200925 & yyyyq < 201025
	replace frac_afs_cmo 				= 0 if frac_afs_cmo == . & yyyyq >= 200925 & yyyyq < 201025				
	bys entity: egen avg_frac_afs_cmo 	= mean(frac_afs_cmo) if yyyyq >= 200925 & yyyyq < 201025
	bys entity: egen temp 				= max(avg_frac_afs_cmo)
	replace avg_frac_afs_cmo 			= temp if avg_frac_afs_cmo == .
	qui: sum frac_afs_cmo if yyyyq >= 200925 & yyyyq < 201025
	replace avg_frac_afs_cmo = `r(mean)' if missing(avg_frac_afs_cmo)
	drop temp frac_afs_cmo

	*Define 
	gen assets_in_afs_cmo_unc		= BHCKG319 if yyyyq >=200925
	replace assets_in_afs_cmo_unc 	= (BHCK1717+BHCK1732)*avg_frac_afs_cmo if yyyyq >=200100 & yyyyq < 200925

* Non Agency Hold-To-Maturity CMO
* Unclear
	*BHCKG316 decomposition
	gen frac_htm_cmo 					= BHCKG316/(BHCKG312+BHCKG316) if yyyyq >= 200925 & yyyyq <= 201025
	replace frac_htm_cmo 				= 0 if frac_htm_cmo == . & yyyyq >= 200925 & yyyyq <= 201025						
	bys entity: egen avg_frac_htm_cmo 	= mean(frac_htm_cmo) if yyyyq >= 200925 & yyyyq <= 201025
	bys entity: egen temp 				= max(avg_frac_htm_cmo)
	replace avg_frac_htm_cmo 			= temp if avg_frac_htm_cmo == .
	qui: sum frac_htm_cmo if yyyyq >= 200925 & yyyyq < 201025
	replace avg_frac_htm_cmo = `r(mean)' if missing(avg_frac_htm_cmo)
	drop temp frac_htm_cmo

	*Define 
	gen assets_in_htm_cmo_unc		= BHCKG316 if yyyyq >=200925
	replace assets_in_htm_cmo_unc 	= (BHCK1714+BHCK1718)*avg_frac_htm_cmo if yyyyq >=200100 & yyyyq < 200925

* Asset-Backed Securities Hold-To-Maturity
* In 
	gen assets_in_htm_abs 		= BHCKC026+BHCKG336+BHCKG340+BHCKG344 if yyyyq >= 200925
		***G336, G340, G344 end 2018Q1 - if updating past 2018Q1 may want to reconsider definition
	replace assets_in_htm_abs 	= BHCKC026 if yyyyq >= 200600 & yyyyq < 200925
	replace assets_in_htm_abs 	= BHCKB838+BHCKB842+BHCKB846+BHCKB850+BHCKB854+BHCKB858 if yyyyq >= 200100 & yyyyq < 200600

*Asset-Backed Securities Available-For-Sale
	gen assets_in_afs_abs 		= BHCKC027+BHCKG339+BHCKG343+BHCKG347 if yyyyq >= 200925
		***G339, G343, G347 end 2018Q1 - if updating past 2018Q1 may want to reconsider definition
	replace assets_in_afs_abs 	= BHCKC027 if yyyyq >= 200600 & yyyyq < 200925
	replace assets_in_afs_abs 	= BHCKB841+BHCKB845+BHCKB849+BHCKB853+BHCKB857+BHCKB861 if yyyyq >= 200100 & yyyyq < 200600


*Other Trading Assets
	gen assets_in_trad_asst_unc 	= BHCM3541 if yyyyq >= 200800
	replace assets_in_trad_asst_unc = BHCK3541 if yyyyq >= 199500 & yyyyq < 200800

*Derivatives 
	gen assets_in_deriv_unc 	= BHCM3543 if yyyyq >= 200800
	replace assets_in_deriv_unc = BHCK3543 + BHFN3543 if yyyyq >= 199500 & yyyyq < 200800

*Life insurance assets
	gen assets_in_insur_unc 	= BHCKK201 + BHCKK202 + BHCKK270
	replace assets_in_insur_unc = 0 if yyyyq < 201100

*Other loans
	gen assets_in_loans_unc 	= BHCKF618
	replace assets_in_loans_unc = 0 if yyyyq < 200800


*TOTALING
gen asset_in 		= assets_in_unchanged + assets_htm_mbs + assets_in_afs_mbs + assets_in_oth_debt + assets_in_htm_abs + assets_in_afs_abs
gen asset_in_unc 	= assets_unc_unchanged + assets_in_htm_mbs_unc + assets_in_afs_mbs_unc + assets_in_trad_mbs_unc + assets_in_afs_cmo_unc + assets_in_htm_cmo_unc + assets_in_trad_asst_unc + assets_in_deriv_unc + assets_in_insur_unc + assets_in_loans_unc	
gen asset_in_frac 	= (asset_in + (1/2)*asset_in_unc)/BHCK2170


******************************************************************************************
*
*					Call Report Merging
*						Want to collect FDIC-insured deposits from each BHC
*						Issue: Sometimes commercial bank will not match to FR-Y9C firm, but WILL
*							connect to that firm's parent. In that case, we match to the FR-Y9C firm,
*							and flag. If we DON'T do this, then those firms
*							will have vastly overestimated connectivities. 
*
******************************************************************************************

*Use FI's master parent matching to see if this BHC has another parent
gen rssd9001 = entity
merge 1:1 qt_dt rssd9001 using ../temp/rssd_hh_match_proced
drop if _merge == 2
gen entity_call = entity
drop _merge

*Try matching on FR-Y9C RSSID
preserve
merge 1:1 qt_dt entity_call using ../temp/InsuredDeposits.dta
keep if _merge == 3
gen matchedon = 1
save ../temp/y9c_callmatch, replace
restore

*For BHCs where that didn't work, try matching based on the BHC's parent's RSSID
preserve 
replace entity_call = id_rssd_top
merge m:1 qt_dt entity_call using ../temp/InsuredDeposits.dta
keep if _merge != 2
gen matchedon = 2
replace matchedon = 3 if _merge == 1
append using ../temp/y9c_callmatch
save ../temp/y9c_callmatch, replace
restore

use ../temp/y9c_callmatch, clear
sort qt_dt entity matchedon

by qt_dt entity: keep if _n == 1

*Add up deposit variables, subtract out our Call report estimates of INSURED deposits
gen tot_depos 	= BHCB2210 + BHCB3187 + BHCB2389 + BHCBHK29 + BHCBJ474
replace tot_depos = BHCB2210 + BHCB3187 + BHCB2389 + BHCB6648 + BHCB2604 if yyyyq < 201700
	/*BHCB6648 and BHCB2604 ended 2016Q4 - these variables were time deposits of more (less) than 100k, 
	  so these two can just be treated as total time deposits. After 2016Q4 time deposits of more(less)
	  than 250k are available and therefore ~should~ perfectly replace the variables					*/

sort entity qt_dt

*Tally up FDIC-insured deposits. Categorize the rest of FR-Y9C deposits from BHC as 
*	uninsured. If strange things happen, make sure uninsured deposits between 0 and sum of
*	all deposits.
gen no_call = 0
replace no_call = 1 if (tot_insur_depos < 0 | missing(tot_insur_depos)) & (entity != 1132449)
replace tot_insur_depos = 0 if (tot_insur_depos < 0 | missing(tot_insur_depos)) & (entity != 1132449)
gen tot_uninsur_depos_y9c = tot_depos - tot_insur_depos

preserve
keep if tot_uninsur_depos_y9c < 0
save ../temp/depos_mismatches.dta, replace
restore

replace tot_uninsur_depos_y9c = 0 if tot_uninsur_depos_y9c < 0

*******************************************************************************************
*
*					Liabilities classifications
*						This one is easier. Many fewer variables, more consistency. 
*
******************************************************************************************

*Fix <- unsure what this comment refers to
egen liab_in 		= rowtotal(tot_uninsur_depos_y9c BHDMB993 BHCKB995 BHCK3548)
egen liab_in_nocall	= rowtotal(tot_depos BHDMB993 BHCKB995 BHCK3548)
egen liab_in_unc 	= rowtotal(BHCK2309 BHCK2332 BHCK2333 BHCKC699 BHCKB557 BHCKB984 BHOD3189 BHOD3187 BHOD2389 BHOD2604 BHOD6648 BHODHK29 BHODJ474)
	/* I throw a wonky assumption here. The sum BHOD2604 and BHOD6648 can be replaced by BHODHK29 and BHCOD474 and their times are mutually exclusive
	 the thought here is that since egen treats missing as 0, there shouldn't be an issue using the sum like this */
gen liabs_in_frac 	= (liab_in + (1/2)*liab_in_unc)/BHCK2948 
gen liabs_in_nocall_frac 	= (liab_in_nocall + (1/2)*liab_in_unc)/BHCK2948 if matchedon == 2

keep if qt_dt <= $charts_end

save ../temp/Contagion_preprocessing, replace
