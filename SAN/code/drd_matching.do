*********************************************************************************
*
*	KMV DRD Database Bankruptcy Search
*	Collin Jones
*	Purpose: Use Moody's Analytics' default database to find BHC bankruptcies, so that
*		we can remove from the final sample firms who have already filed for bankruptcy.
*
*********************************************************************************


clear

if "`c(os)'" == "Unix" {
	set odbcmgr unixodbc
	local connect `""DSN=WDM""'
}

else if "`c(os)'" == "Windows"{
	local connect `"DRIVER={SQL Server};SERVER=m1-wdmsql01.frb.gov;DATABASE=WDM;"'
	
}
odbc load, exec("SELECT * FROM  vw_moodys_drd_issr_ids") connectionstring(`connect')
save ../temp/drd_matches.dta, replace

use ../temp/CUSIP_MKMVID.dta, clear
duplicates drop mkmv_id, force

collapse (lastnm) issuer_cusip_id tkr company_name, by(mkmv_id)
drop if missing(company_name)
save ../temp/matcher_cusipmkmvid, replace

preserve
	rename mkmv_id issr_id
	merge 1:m issr_id using ../temp/drd_matches.dta
	rename issr_id mkmv_id
	keep if _merge == 3
	gen type = 1
	save ../temp/matched, replace
restore

preserve

rename issuer_cusip_id issr_id
merge m:m issr_id using ../temp/drd_matches.dta
rename issr_id issuer_cusip_id
keep if _merge == 3
gen type = 2
append using ../temp/matched
save ../temp/matched, replace
restore

preserve
rename tkr issr_id
drop if mkmv_id == "N11845"
merge m:m issr_id using ../temp/drd_matches.dta
rename issr_id tkr
keep if _merge == 3
gen type = 3
append using ../temp/matched

save ../temp/matched, replace
restore

clear

clear
odbc load, exec("SELECT * FROM vw_moodys_drd_mast_issr") connectionstring(`connect')
save ../temp/mast_issr, replace

clear
odbc load, exec("SELECT * FROM vw_moodys_drd_family_structure") connectionstring(`connect')
save ../temp/family, replace

clear
odbc load, exec("SELECT * FROM  vw_moodys_drd_mast_dflt ") connectionstring(`connect')
merge m:1 mast_issr_num using ../temp/mast_issr
keep if _merge == 3
drop _merge


sort bankruptcy_datetime
keep if dofc(bankruptcy_datetime) > td(01jan2002)
keep if bankruptcy_typ_cd != ""

gen company_name = subinstr(issuer_nam, ".", "", 10)
replace company_name = upper(company_name)
replace company_name = strtrim(company_name)

merge m:m company_name using ../temp/matcher_cusipmkmvid.dta
keep if _merge == 3
gen type = 4
keep mkmv_id tkr company_name mast_issr_num type
append using ../temp/matched
sort type
drop _merge
collapse (firstnm) company_name tkr type issuer_cusip_id issr_id_num , by(mkmv_id mast_issr_num)
drop if mkmv_id == "792860"
drop if missing(mast_issr_num)

save ../temp/matched, replace

clear

odbc load, exec("SELECT * FROM  vw_moodys_drd_mast_dflt ") connectionstring(`connect')
save ../temp/kmv_bankruptcies, replace
keep if bankruptcy_typ_cd != ""
duplicates drop bankruptcy_datetime mast_issr_num resolution_datetime, force
keep mast_issr_num bankruptcy_datetime resolution_datetime blurb

merge m:m mast_issr_num using ../temp/matched
replace bankruptcy_datetime = dofc(bankruptcy_datetime)
replace resolution_datetime = dofc(resolution_datetime)
format bankruptcy_datetime %td
format resolution_datetime %td
replace resolution_datetime = td(01jan2050) if missing(resolution_datetime)

keep if _merge == 3

keep mkmv_id bankruptcy_datetime resolution_datetime
keep if bankruptcy_datetime >= td(01jan2002) & bankruptcy_datetime <=td(01jan2017)
duplicates drop mkmv_id bankruptcy_datetime, force

bys mkmv_id: gen count = _n
reshape wide bankruptcy_datetime resolution_datetime, i(mkmv_id) j(count)

save ../temp/bankruptcies.dta, replace
