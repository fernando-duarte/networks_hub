/*
* This code artificially extends RSSD_hh_match_all_original to 2018
* Run Locally, should be quick
*/
clear all
set more off


use "\\rb.win.frb.org\B1\NYRESAN\RDS\Work\cmf\b1far01\Fernando\Network_Contagion\input\rssd_hh_match_all_original.dta", replace 

gen foo = 2 if date_q==230
expand foo
sort id_rssd id_rssd_top date_d_start date_d_end date_q
replace date_q=date_q + 1 if foo==foo[_n-1] & id_rssd_top==id_rssd_top[_n-1] & id_rssd==id_rssd[_n-1] & date_q==date_q[_n-1]
drop foo

save "\\rb.win.frb.org\B1\NYRESAN\RDS\Work\cmf\b1far01\Fernando\Network_Contagion\input\rssd_hh_match_all_extended.dta"
