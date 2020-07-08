/************************************************************

Aggregate Sector Composition 
Purpose: This code counts the number of each type of 
	firm in each aggregate node for each quarter

This code is extremely short but I'm not sure if it makes sense
	to include in another file

*****************************************************************/


*UNCOMMENT BLOCK BELOW IF RUNNING THIS CODE SEGMENT INDEPENDENTLY
/****************************************************************

	clear all
	
	if "`c(os)'" == "Unix" {
	   global san "/san"
	}
	else if "`c(os)'" == "Windows" {
	   global san "\\rb.win.frb.org\B1\NYRESAN"
	}
	
	global charts_start = yq(2002, 1)
	global charts_end	= yq(2017, 4)

****************************************************************/
	
/****************************************************************
Broker Node Info
****************************************************************/

/* use $san/RDS/Work/cmf/b1far01/Fernando/Network_Contagion/temp/broker_nodes.dta, clear  //This file finalized in Model_series_processing.do */
use ../temp/broker_nodes_config$bd_config.dta, clear
	
	gen beta 	= liabs_in_frac
	gen w 		= BHCK2170*10^(-3) - BHCK2948*10^(-3)
	gen c 		= (BHCK2170*10^(-3) - (asset_in*10^(-3) + (1/2)*asset_in_unc*10^(-3))) // 
	gen p_bar	= BHCK2948*10^(-3)
	gen b 		= (BHCK2948 - (liab_in + (1/2)*liab_in_unc))*(10^(-3))
	rename prob_default_phys delta
	capture confirm variable assets_total
	gen lambda 	= c/w

	gen assets 	= BHCK2170*10^(-3)
	if _rc != 0{
		gen assets_total = .
	}
	replace assets_total = assets_total*10^(-3)
	replace assets_total = assets if missing(assets_total)
	
	format qt_dt %9.0g

	keep nm_short tkr p_bar assets c b delta beta w qt_dt  
	order nm_short tkr qt_dt delta assets p_bar b c w beta   
    /* export excel $san/RDS/Work/cmf/b1far01/Fernando/Network_Contagion/temp/node_stats_forsimulation_all.xls, sheet("Dealer Aggregates") firstrow(variables) sheetreplace */
	export excel ../temp/node_stats_forsimulation_all.xls, sheet("Dealer Aggregates") firstrow(variables) sheetreplace


*open firm list generated by compile_agg_sector_data
/* use $san/RDS/Work/cmf/b1far01/Fernando/Network_Contagion/temp/firm_list, clear */
use ../temp/firm_list, clear

gen Dealer 		=(sector=="Broker Dealers")
gen Insurance 	=(sector=="Insurance")
gen Other 		=(sector=="Other")
gen REITs 		=(sector=="REITs")

collapse (sum) Dealer Insurance Other REITs, by(qt_dt)

/* save $san/RDS/Work/cmf/b1far01/Fernando/Network_Contagion/temp/aggregate_composition_count.dta, replace */
/* export excel $san/RDS/Work/cmf/b1far01/Fernando/Network_Contagion/temp/node_stats_forsimulation_all.xls, sheet("Aggregates Composition") firstrow(variables) sheetreplace */

save ../temp/aggregate_composition_count.dta, replace
export excel ../temp/node_stats_forsimulation_all.xls, sheet("Aggregates Composition") firstrow(variables) sheetreplace





/* Robustness checks that aren't necessary to run *//*

*Same idea but directly using the DFS data rather than the processed firm list

use ../temp/DFS.dta, clear

duplicates drop qt_dt mkmv_id mkmvind_name_desc, force

keep if mkmvind_name_desc == "SECURITY BROKERS & DEALERS" | mkmvind_name_desc == "INVESTMENT MANAGEMENT" | mkmvind_name_desc == "INSURANCE - LIFE" | /// 
 		mkmvind_name_desc == "INSURANCE - PROP/CAS/HEALTH" | mkmvind_name_desc == "REAL ESTATE INVESTMENT TRUSTS" | mkmvind_name_desc == "FINANCE NEC" | mkmvind_name_desc == "FINANCE COMPANIES"

keep if kmv_country_id == "USA" | kmv_country_id == "BRU" | kmv_country_id == "CYM"
keep if qt_dt >= $charts_start

*Sector labelling
gen sector = ""
replace sector = "Broker Dealers" if mkmvind_name_desc == "SECURITY BROKERS & DEALERS"
replace sector = "Insurance" if mkmvind_name_desc == "INSURANCE - LIFE" | mkmvind_name_desc == "INSURANCE - PROP/CAS/HEALTH"
replace sector = "Mutual" if mkmvind_name_desc == "INVESTMENT MANAGEMENT"
replace sector = "REITs" if mkmvind_name_desc == "REAL ESTATE INVESTMENT TRUSTS"
replace sector = "Other" if mkmvind_name_desc == "INVESTMENT MANAGEMENT" | mkmvind_name_desc == "FINANCE NEC" | mkmvind_name_desc == "FINANCE COMPANIES"

*collapse
gen Dealer=(sector=="Broker Dealers")
gen Insurance=(sector=="Insurance")
gen Other=(sector=="Other")
gen REITs=(sector=="REITs")

collapse (sum) Dealer Insurance Other REITs, by(qt_dt)

save $san/RDS/Work/cmf/b1far01/Fernando/Network_Contagion/temp/aggregate_composition_count2.dta, replace

*Why are these two identical???
merge 1:1 qt_date using $san/RDS/Work/cmf/b1far01/Fernando/Network_Contagion/temp/aggregate_composition_count.dta

*/
