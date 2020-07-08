*********************************************************************************
*
*	RSSID-MKMVID Matching
*	Collin Jones and Rui Yu
*	Purpose: Link Y9C's RSSID field to KMV's ID values
*
*********************************************************************************

clear
set more off

*Merge with FI reg data to obtain RSSID - CUSIP match
	*NOTE: when running, FI was updating reg data - commented out data should be right generally
*use $regData/fr_y9c/qbhc_nc_clean.dta, clear

if `"$MY_ENV"' ~= "RAN"{
	use $regData/fr_y9c/archive/1803/qbhc_nc_clean.dta, clear
}
else{
	use ../temp/qbhc_nc_all.dta, clear
}

gen qt_dt 	= yq(year, quarter)
keep entity qt_dt id_cusip name
drop if id_cusip == "0"
duplicates drop entity id_cusip, force
save ../temp/FI_Link.dta, replace

*Match PERMCO to CUSIP
import delimited "https://www.newyorkfed.org/medialibrary/media/research/banking_research/data/crsp_$entity_permco_date.csv?la=en", clear
	*Troubleshooting replacement HARDCODED
	replace dt_end = 20171231 if dt_end==20161231

drop notice
drop in 1/2

*New version of FI's permco-cusip doesn't seem to have this. Don't think it impacts anything.
//keep if inst_type == "Bank Holding Company"

tostring dt_start, replace
tostring dt_end, replace
gen start_yr 	= substr(dt_start, 1, 4)
gen end_yr 		= substr(dt_end, 1, 4)
gen start_month = substr(dt_start, 5, 2)
gen end_month 	= substr(dt_end, 5, 2)
destring start_yr, replace
destring end_yr, replace
destring start_month, replace
destring end_month, replace
gen start_qt 		= 1 if start_month >= 1 & start_month <= 3
replace start_qt 	= 2 if start_month >= 4 & start_month <= 6
replace start_qt 	= 3 if start_month >= 7 & start_month <= 9
replace start_qt 	= 4 if start_month >= 10 & start_month <= 12
gen end_qt 			= 1 if end_month >= 1 & end_month <= 3
replace end_qt 		= 2 if end_month >= 4 & end_month <= 6
replace end_qt 		= 3 if end_month >= 7 & end_month <= 9
replace end_qt 		= 4 if end_month >= 10 & end_month <= 12
gen start_qtdt 		= yq(start_yr, start_qt)
gen end_qtdt 		= yq(end_yr, end_qt)

//drop if end_qtdt < $charts_start

keep name entity permco start_qtdt end_qtdt

gen quarters = end_qtdt - start_qtdt + 1
expand quarters
bys entity permco start_qtdt end_qtdt: gen count = _n
replace count = count - 1
gen qt_dt = start_qtdt + count

keep name entity permco qt_dt
duplicates drop entity permco qt_dt, force
save ../temp/RSSID_PERMCO.dta, replace

*Input PERMCO-CUSIP link
insheet using ../temp/permco_cusip.csv, comma clear
duplicates drop permco cusip, force
keep permco cusip

*merge with Entity-PERMCO
joinby permco using ../temp/RSSID_PERMCO.dta

keep permco cusip name entity qt_dt
order entity permco cusip name 
duplicates drop entity cusip qt_dt, force

*Combine two Entity-CUSIP matches

append using ../temp/FI_Link.dta

*First preference for cusip from RSSID-PERMCO-CUSIP match
gen cusip_6 = substr(cusip, 1, 6)

//Fill missing cusip from FI reg data match
replace cusip_6 = substr(id_cusip, 1, 6) if cusip_6 == "" | cusip_6 == "0" 
drop permco cusip id_cusip 
rename cusip_6 cusip
order entity cusip name
duplicates drop entity cusip qt_dt, force

*BONY Mellon seems to break in this process. That's the only one I've spotted, though.
replace cusip = "064058" if cusip == "640581"

save ../temp/RSSID_CUSIP.dta, replace

*Match CUSIP-MKMVID
use ../temp/CUSIP_MKMVID.dta, clear
rename issuer_cusip_id cusip
joinby cusip using ../temp/RSSID_CUSIP.dta, unmatched(none)
duplicates drop entity mkmv_id, force
keep mkmv_id entity name tkr
save ../temp/RSSID_MKMVID.dta, replace
