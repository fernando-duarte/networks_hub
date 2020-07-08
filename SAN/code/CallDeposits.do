
*********************************************************************************
*
*	Call Report pull and analysis
*	Collin Jones and Rui Yu
*	Purpose: Find the FDIC-insured deposits of the commercial bank entities within each BHC
*
*********************************************************************************

clear
clear matrix
set more off
set maxvar 9000

*Import Call Reports
******use $regData/call/qbank_nc_clean.dta, clear
	*NOTE: when running, FI was updating reg data - commented out data should be right generally
*use $regData/call/qbank_nc_all.dta, clear
if `"$MY_ENV"' ~= "RAN"{
	use $regData/call/archive/1803/qbank_nc_all.dta, clear	
	}
else{
	use ../temp/qbank_nc_all.dta, clear
	}

*Separate dates
tostring date, replace
gen qt_dt 	= yq(year, quarter)
keep if qt_dt <= $charts_end
keep if qt_dt >= $charts_start
destring date, replace

******************************************************************
*
*	So what's going on here?
*		The fields used below correspond to the following:
*			fdic_uninsur_depos_1 = Call report variable rconf051
*				Amount in accounts > $250k
*			fdic_uninsur_depos_2 = Call report variable rconf052
*				# of accounts > $250k
*			fdic_uninsur_depos_3 = Call report variable rconf047
*				Amount in retirement accounts > $250k
*			fdic_uninsur_depos_4 = Call report variable rconf048
*				# of retirement accts > $250k
*			fdic_uninsur_depos 	= Call report variable rcon5597
*				Estimated total amount of uninsured deposits
*
*	If we have fdic_uninsur_depos, use that. Otherwise, try to 
*		construct it using the variables that we do have.
*
******************************************************************


******************************************************************
*
*	Using FI's master list of umbrella BHCs
*
******************************************************************

preserve
*use ../input/rssd_hh_match_all_extended, clear
use ../input/rssd_hh_match_all2, clear
sort id_rssd date_q date_d_end
collapse (lastnm) id_rssd_top, by(id_rssd date_q)
rename id_rssd rssd9001
rename date_q qt_dt
save ../temp/rssd_hh_match_proced, replace
restore

merge 1:1 qt_dt rssd9001 using ../temp/rssd_hh_match_proced
keep if _merge == 3
gen entity_call = id_rssd_top

*Define uninsured deposits
gen uninsur_depos 		= fdic_uninsur_depos

******************************************************************
*
*	Section only relevent if we're missing fdic_uninsur_depos...
*
******************************************************************

*Calculate uninsured 2008Q4 and before insurance limit change
replace uninsur_depos 	= fdic_unins_depos_3 - 100*(fdic_unins_depos_4 + fdic_unins_depos_2) + fdic_unins_depos_1 if uninsur_depos ==. & date <= 20081231

*Calculate uninsured from 2009Q1 to 2009Q3 after limit change but before form change
replace uninsur_depos 	= fdic_unins_depos_3 - 100*(fdic_unins_depos_4 + fdic_unins_depos_2) + fdic_unins_depos_1 if uninsur_depos ==. & date >= 20090331 & date <= 20090630

*Calculate uninsured after 2009Q3 after limit and form change
replace uninsur_depos 	= fdic_unins_depos_3 - 250*(fdic_unins_depos_4 + fdic_unins_depos_2) + fdic_unins_depos_1 if uninsur_depos ==. & date >= 20090930

*Fill missing uninsured with 0's
replace uninsur_depos 	= 0 if uninsur_depos == .

******************************************************************
*
*	Some final clean-up
*
******************************************************************

*In case multiple commercial banks in each BHC 
bys date entity_call: egen tot_uninsur_depos 	= sum(uninsur_depos)

*Fernando: This is where we implement our domestic insured deposits fixed. Makes sure
*	the insured deposits number is internally consistent with call report numbers.
bys date entity_call: egen tot_insur_depos = sum(domestic_dep - uninsur_depos)

duplicates drop date entity_call, force

keep qt_dt entity_call tot_uninsur_depos tot_insur_depos

save ../temp/InsuredDeposits.dta, replace
