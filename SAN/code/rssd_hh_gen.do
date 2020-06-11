********************************************************************************
*
*	rssd_hh_match_all generation
*	Francisco Ruela
*	Purpose: Convert reg_data raw_toph.dta into quarterly RSSID matches
*	Input: reg_data -> raw_toph.dta
*	Output: input/rssd_hh_match_all
*
*
********************************************************************************

clear all
set more off

if "`c(os)'" == "Unix" {
   global san "/san"
}
else if "`c(os)'" == "Windows" {
   global san "//rb.win.frb.org/B1/NYRESAN"
}

*open FIs raw data
use "$san/RDS/Derived/reg_data_nc/nic/output/raw_toph.dta"

****************************************************
***		Apply FI's Filters		 ***
****************************************************

keep if holder_type=="REG"
keep if UNIQUE_TIER_IND =="Y"
*keep if ctrl_flg=="Y"

keep id_rssd d_dt_start d_dt_end id_rssd_top

****************************************************
***		  Convert to Long		 ***
****************************************************

gen q1 = yq(year(d_dt_start), quarter(d_dt_start))
	format q1 %tq
gen q2 = yq(year(d_dt_end), quarter(d_dt_end))
	replace q2 = yq(2018,3) if q2 == yq(9999,4)
	format q2 %tq
gen foo = q2-q1+1

expand foo

bys id_rssd id_rssd_top q1 q2 foo: gen date_q=q1-1+_n 
	replace date_q = q2 if foo==1
	format date_q %tq
	
****************************************************
***		Format Appropriately		 ***
****************************************************

drop foo q1 q2
*duplicates drop id_rssd id_rssd_top date_q, force

rename d_dt_start  date_d_start
rename d_dt_end  date_d_end
order id_rssd id_rssd_top 

save "$san/RDS/Work/cmf/b1far01/Fernando/Network_Contagion/input/rssd_hh_match_all2.dta", replace
*merge m:m id_rssd id_rssd_top date_q using "$san/RDS/Work/cmf/b1far01/Fernando/Network_Contagion/input/rssd_hh_match_all_original.dta"
