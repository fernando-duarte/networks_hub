*********************************************************************************
*
*	Network data updating
*	Collin Jones and Rui Yu
*	Purpose: Pull data necessary for further analysis. Include permco-rssid matches
*		from FI, probabilites of default from KMV, and KMV company identifiers
*

*	Most data pulls will only refresh data when the file is missing in its directory.
*		So delete those files before running, if you want to update everything
*
*   NOTE Raw9YC.stcmd needs to be updated to save Raw9YC to YOUR temp folder
*
* 	Files Updated:
* 		RawDFS 				
*		RawDFS8
*		CUSIP_MKMVID
*		CUSIP_MKMVID8
*		RawY9C
*		permco_cusip
*********************************************************************************

clear

local pwd: pwd

*permco_cusip, RSSID_MKMVID, RSSID_PERMCO are temporarily not being removed in the
*remove temp option. The RSSID data are hand modified by Fernando, the permco_cusip
*is held because WRDS isn't complying

if "`c(os)'" == "Unix" {
	set odbcmgr unixodbc
	local connect `""DSN=WDM""'
	if $wipe_temp == 1{
		cd ../temp/
		shell mv ffunds.dta ../
		shell mv permco_cusip.csv ../
		shell mv RSSID_MKMVID.dta ../
		shell mv RSSID_PERMCO.dta ../
		shell rm *.xls
		shell rm *.dta
		shell mv ../ffunds.dta ./
		shell mv ../permco_cusip.csv ./
		shell mv ../RSSID_MKMVID.dta ./
		shell mv ../RSSID_PERMCO.dta ./
		cd ../code
	}

}

else if "`c(os)'" == "Windows"{
	local connect `"DRIVER={SQL Server};SERVER=m1-wdmsql01.frb.gov;DATABASE=WDM;"'
	if $wipe_temp == 1{
		cd ../temp/
		shell erase *.xls
		shell erase *.dta
		cd ../code
	}
	
}

 //shell /data/apps/Anaconda2-5.0.1/bin/python ../code/format_regs.py
// Install STCMD Stat-transfer Stata utility
//cap which stcmd
//if _rc ssc install stcmd, replace

*Pull FI Web entity-permco Link - Now done in Match_RSSID.do
//copy https://www.newyorkfed.org//research/banking_research/datasets/crsp_$entity_permco_date.csv ///
//	../temp/entity-permco_$entity_permco_date.csv, replace

*Execute SQL pull from KMV database
capture confirm file ../temp/RawDFS.dta
if _rc != 0{
	clear 
	odbc load, exec("SELECT mkmv_id,annualized_edf_01yr_pct, risk_neutral_annualized_qdf_01yr_pct, dt_day,mkmvind_name_desc, total_assets_current_amt, kmv_country_id, total_liabilities_amt FROM vw_kmv_v9_dfs_v8map") ///
		connectionstring(`connect')
	rename var3 annualized_qdf_01yr_pct
	save ../temp/RawDFS.dta, replace
}

*In case we want to test the V9 vs V8 difference...
capture confirm file ../temp/RawDFS8.dta
if _rc != 0{
	clear
	odbc load, exec("SELECT mkmv_id,annualized_edf_01yr_pct, risk_neutral_annualized_qdf_01yr_pct, dt_day,mkmvind_name_desc, total_assets_current_amt, total_liabilities_amt, kmv_country_id FROM vw_kmv_v8_dfs") ///
		connectionstring(`connect')
	rename var3 annualized_qdf_01yr_pct
	save ../temp/RawDFS8.dta, replace
}

*Pull MKMV-CUSIP match (also incl. stock market tickers, which will be useful later)
clear
capture confirm file ../temp/CUSIP_MKMVID.dta
if _rc != 0{
	odbc load, exec("SELECT * FROM vw_kmv_v9_cmf_v8map") connectionstring(`connect')
	drop if issuer_cusip_id == "" | mkmv_id == ""
	duplicates drop issuer_cusip_id mkmv_id, force
	save ../temp/CUSIP_MKMVID.dta, replace
}

*In case we want to test the V9 vs V8 difference...
capture confirm file ../temp/CUSIP_MKMVID8.dta
if _rc != 0{
	clear
	odbc load, exec("SELECT mkmv_id,issuer_cusip_id,tkr,company_name FROM vw_kmv_v8_cmf") connectionstring(`connect')
	drop if issuer_cusip_id == "" | mkmv_id == ""
	duplicates drop issuer_cusip_id mkmv_id, force
	save ../temp/CUSIP_MKMVID8.dta, replace
}

*Pull Y9C Data, if missing. Send batch job to head node if working on the san, run with local stattransfer if not.
capture confirm file ../temp/RawY9C.dta
if _rc != 0{
	if "`c(os)'" == "Windows"{
		stcmd $regData/fr_y9c/cuv_bhcf_nc.sas7bdat $san/RDS/Work/cmf/b1far01/Fernando/Network_Contagion/input/RawY9C.dta
	}
	else if "`c(os)'" == "Unix"{
		
		shell ssh b1rsflx23.re.ny.frb.org stattransfer13-batch-withemail 20 `pwd'/RawY9C.stcmd
		local wait_time = 1000*60*8
		sleep `wait_time'
	}
}

*Update permco-cusip links from WRDS
capture confirm file ../temp/permco_cusip.csv
if _rc != 0{
	if "`c(os)'" == "Unix"{
		if $changed_WRDS_files == 1{
			shell /data/apps/Anaconda2-5.0.1/bin/python WRDS_meta.py 
		}
		else if $changed_WRDS_files == 0{
			disp "Won't continue, set global changed_WRDS_files to 1 if you really want to do this"
			assert $changed_WRDS_files == 1
		}
		
	}
}


use ../temp/RawY9C, clear

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
keep if qt_dt <= $charts_end
keep entity qt_dt
rename entity entity_test
save ../temp/valid_entities.dta, replace
