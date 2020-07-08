*********************************************************************************
*
*	KMV Data Cleaning
*	Authors: Collin Jones and Rui Yu
*	Updated(6/16/20): Separated into separate file from Analysis_Y9C
*	Purpose: Execute variable categorizations defined in codebook to get numbers for each BHC's 
*		assets in and out of system. Also: clean and merge with KMV data.
*
*********************************************************************************
clear
set more off

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
