*********************************************************************************
*
*	Network Contagion MAIN
*	Collin Jones and Rui Yu
*	Purpose: Execute all NVI code sequentially
*
*	Note: If you are about to run this for the first time 
*		1) update WRDS_meta and WRDS_query
*		2) Update RawY9C to extract to the correct folder
*		3) Make sure all the packages are installed correctly
*		4) I would recomend running the code once with wipe_temp=0 to make sure the code runs
*	Note: If you are trying to extend the data
*		1) Update the charts_end macro
*		2) Update Analysis_Y9C.do (non-trivial)
*		3) Update process_FOCUS to pull most recent SIFMA data 
*
*********************************************************************************

clear all
capture cd code/
local pwd: pwd
capture log close
set maxvar 8000

global runthrough_tag rerun_FAR

log using ../code/log_$runthrough_tag.txt, replace
log close

***********************************************************
*
*	Manual steps, to be done occassionally:
*		(1) occassionally run ffunds.do on your local machine. Make sure to verify that Haver is correctly pathed.
*		(2) Periodically (quarterly, Rui seemed to think) pull a PERMCO-to-CUSIP 
*			match from WRDS. Can be done using SQL query select cusip, permco from CRSPM.STOCKNAMES (or equivalent)
*			Must be added manually, it will NOT be cleared from wipe_temp
*
***********************************************************

*Packages/commands you'll need:
* help labmask 
* help grc1leg
* ssc install egenmore
* ssc install stcmd
* ssc instal estout

*Set to 1 if you want to update EVERYTHING in temp folder (will drastically increase runtime)
global wipe_temp 0

*IMPORTANT: The contents of WRDS_meta.py will automatically pull the permco-cusip match for you. 
*	BUT if the usernames and passwords in that file are incorrect, then you could
*	lock the floor's remote WRDS access through too many authentication errors. 
*	If you're definitely sure your usernames and passwords are right, set this to 1.
global changed_WRDS_files 1

*Are we doing a full update and chart regeneration, or just outputting data for simulations?
*0 for full update, 1 for just outputting data
global output_simulation_data 1

* Standardized variable list for simulation data output
global simulation_data_vars nm_short qt_dt tkr delta delta_alt beta w c assets nvi_benchmark p_bar b

*Which different bankruptcy costs do we want to try in the NVI? (in 5)
global gammas 1 5 10 15 20 30

*What level of bankruptcy costs are the "benchmark"? (in %)
global gamma_benchmark 0

*What are the RSSIDs of the BHCs we want to select the Beta+ from? Currently,
*	the top 10, top 11-25 aggregate broker dealer nodes are also appended to 
*	this list
global bankSample 1039502 1068025 1069778 1070345 1073757 1074156 ///
	1111435 1119794 1120754 1131787 1275216 1562859 1951350 2162966 ///
	2277860 2380443 3242838 3587146 3232316 1132449

*What probability of default do you want for the fixed-delta robustness exercise? (percent in decimal form)
global delta_fixed 0.06

*Needed for installing packages on san
if "`c(os)'"~="Windows" & `"$MY_ENV"' != "RAN"{
	set httpproxyhost "p1web1.frb.org"
	set httpproxyport 8080
	set httpproxy on

	set odbcmgr unixodbc
}

global  MY_ENV: env MY_ENV
if `"$MY_ENV"' == ""{
    global MY_ENV SAN
}

capture ssc install labutil
capture ssc install egenmore
capture net get dm88_1.pkg 
capture net install grc1leg, from(http://www.stata.com/users/vwiggins)

*For charts
global charts_start 	= yq(2002, 1)
global charts_end		= yq(2017, 4)

*For quarter-specific sum stats
global snapshot_date 	= yq(2017, 4)

if "`c(os)'" == "Unix" {
   global san "/san"
}
else if "`c(os)'" == "Windows" {
   global san "//rb.win.frb.org/B1/NYRESAN"
}

*Location of FI's reg data folder
global regData "$san/RDS/Derived/reg_data_nc"
global entity_permco_date 	= "20181231"

/* How we deal with broker dealer asset double counting: */
/* 0 -- Nothing. GS, MS, and the like in both Y9C sample and Top 10 Dealer node. */
/* 1 -- Maximum data detail. Switch GS, MS assets to Y9C node. Remove BHCKC252 values of KMVID overlap from Top 10 Node. */
/* 2 -- Maximum consistency. Keep GS, MS in Top 10 dealer node. ADD BHCK2170 - BHCKC252 to Top 10 Node. */
/* 3 -- Maximum data detail and double counting avoidance. Switch GS, MS assets to Y9C node. Remove BHCKC252 values of ALL Y9C sample BHCs. */

global bd_config 3

*Location of Haver Analytics data folder
global haverData "K:/DLX/data" 

* Run this do file locally, once a quarter (others should be run on san)
//do ../code/ffunds.do

if $output_simulation_data == 0 {
	do ../code/Update_Data.do
    do ../code/Match_RSSID_MKMVID.do
    do ../code/CallDeposits.do
    if `"$MY_ENV"' ~= "RAN"{
        do ../code/drd_matching.do
    }
    do ../code/KMV_clean.do
    do ../code/Analysis_Y9C.do
	shell /data/apps/Anaconda2-5.0.1/bin/python `code'/process_FOCUS.py
}
*unsure what the code macro does above - FR
do ../code/compile_agg_sector_data.do
do ../code/Model_series_processing.do
do ../code/agg_sector_composition.do

if $output_simulation_data == 0{
	do ../code/Plots_Paper.do
}
//do ../code/Plots_Appendix.do

shell /data/apps/Anaconda2-5.0.1/bin/python clean_sumstats.py
