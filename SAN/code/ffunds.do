*********************************************************************************
*
*	Flow of Funds Analysis
*	Collin Jones
*	Purpose: Import flow of funds data from Haver, classify aggregate sector nodes
*		as inside or outside financial system, save results for use 
*		in compile_agg_sector_data.dta
*

*   Note: This is the only do-file that MUST be run locally (import haver does not
*	work otherwise). 
*   Note 2: Make sure you have all packages listed in NtwkCntagn.do (this may be 
*	an incomplete list - if you encounter an additional package please add to this
* 	list and the list in NtwkCntagn.do):
*		labmask
*		egenmore
*		grc1leg
*		stcmd
*		xpose

*********************************************************************************

* set folder where haver data is - note globals are copied from NtwkCntagn.do
global haverData "K:/DLX/data" 
set haverdir "$haverData"

*For charts
global charts_start 	= yq(2002, 1)
global charts_end	= yq(2018, 1)



*Insurance Companies
*Unclear
clear
import haver oa54cof5@ffunds oa54cog5@ffunds oa54dpp5@ffunds oa54mfs5@ffunds ///
	oa54mms5@ffunds oa54otl5@ffunds oa54una5@ffunds oa51cof3@ffunds ///
	oa51cog5@ffunds oa51dpp3@ffunds oa51mfs3@ffunds oa51mms3@ffunds oa51una5@ffunds

renvars oa*, map("insurance_unc_" + "@")
rename time qt_dt
save ../temp/ffunds.dta, replace

*In
clear
import haver oa51ahy3@ffunds oa51ccd3@ffunds oa54ahv3@ffunds oa54ccd5@ffunds
renvars oa*, map("insurance_in_" + "@") 
rename time qt_dt
merge 1:1 qt_dt using ../temp/ffunds.dta, nogen
save ../temp/ffunds.dta, replace

*Brokers and Dealers
*Unclear
clear
import haver OA66COF5@FFUNDS OA66COG3@FFUNDS OA66DPP3@FFUNDS OA66MFT5@FFUNDS OA66OTL5@FFUNDS
renvars oa*, map("dealer_unc_" + "@")
rename time qt_dt
merge 1:1 qt_dt using ../temp/ffunds.dta, nogen
save ../temp/ffunds.dta, replace

*In
clear
import haver OA66CCD3@FFUNDS OA66AHY3@FFUNDS 
renvars oa*, map("dealer_in_" + "@")
rename time qt_dt
merge 1:1 qt_dt using ../temp/ffunds.dta, nogen
save ../temp/ffunds.dta, replace

*Mutual Funds
*Unclear
clear
import haver oa65cof5@ffunds oa65cog0@ffunds oa65dpp0@ffunds oa65una5@ffunds
renvars oa*, map("mutual_unc_" + "@")
rename time qt_dt
merge 1:1 qt_dt using ../temp/ffunds.dta, nogen
save ../temp/ffunds.dta, replace

*In
clear
import haver oa65ahy3@ffunds
renvars oa*, map("mutual_in_" + "@")
rename time qt_dt
merge 1:1 qt_dt using ../temp/ffunds.dta, nogen
save ../temp/ffunds.dta, replace

*Depository Institutions
*Unclear
clear
import haver OA70CPA5@FFUNDS OA70COF5@FFUNDS OA70OTL5@FFUNDS OA70COG5@FFUNDS OA70MFS5@FFUNDS OA70MFT5@FFUNDS 
renvars oa*, map("depository_unc_" + "@")
rename time qt_dt
merge 1:1 qt_dt using ../temp/ffunds.dta, nogen
save ../temp/ffunds.dta, replace

*In
clear
import haver OA70FFS5@FFUNDS OA70BLN5@FFUNDS
renvars oa*, map("depository_in_" + "@")
rename time qt_dt
merge 1:1 qt_dt using ../temp/ffunds.dta, nogen
save ../temp/ffunds.dta, replace

*REITs
*Unclear
clear
import haver OA64COF5@FFUNDS OA64MFT5@FFUNDS
renvars oa*, map("reit_unc_" + "@")
rename time qt_dt
merge 1:1 qt_dt using ../temp/ffunds.dta, nogen
save ../temp/ffunds.dta, replace

*In
clear
import haver OA64CCD5@FFUNDS 
renvars oa*, map("reit_in_" + "@")
rename time qt_dt
merge 1:1 qt_dt using ../temp/ffunds.dta, nogen
save ../temp/ffunds.dta, replace
	
*Other
clear
import haver OA47SLM3@FFUNDS OA47COF5@FFUNDS OA47MFT5@FFUNDS OA61COF3@FFUNDS OA61FLB0@FFUNDS ///
	OA61MFT5@FFUNDS OA50DPP5@FFUNDS OA50COF5@FFUNDS OA50CMB5@FFUNDS OA50COG5@FFUNDS OA67OTL5@FFUNDS OA67AIW3@FFUNDS
renvars oa*, map("other_unc_" + "@")
rename time qt_dt
merge 1:1 qt_dt using ../temp/ffunds.dta, nogen
save ../temp/ffunds.dta, replace

*In
clear
import haver  OA47SRP3@FFUNDS OA61CCD3@FFUNDS OA61TID3@FFUNDS OA50AHY3@FFUNDS OA50JUM5@FFUNDS
renvars oa*, map("other_in_" + "@")
rename time qt_dt
merge 1:1 qt_dt using ../temp/ffunds.dta, nogen
save ../temp/ffunds.dta, replace

	
*Temporary output	
save ../temp/ffunds.dta, replace

*********************************************************************************
*		Create summary stats of % of sector-wide assets attributable
*			to each of some broad categories. 
*********************************************************************************

clear
import haver OA51CCD3@FFUNDS OA51MMS3@FFUNDS OA51AHY3@FFUNDS OA51DPP3@FFUNDS ///
	OA51TRE3@FFUNDS OA51AGI3@FFUNDS OA51STO3@FFUNDS OA51COF3@FFUNDS OA51CMT3@FFUNDS ///
	OA51COG5@FFUNDS OA51MFS3@FFUNDS OA51TRD3@FFUNDS OA51EQG3@FFUNDS OA51MFA3@FFUNDS ///
	OA51UNA5@FFUNDS OA54TAO5@FFUNDS OA54CCD5@FFUNDS OA54MMS5@FFUNDS OA54AIJ5@FFUNDS ///
	OA54DPP5@FFUNDS OA54TRE5@FFUNDS OA54AGI5@FFUNDS OA54STO5@FFUNDS OA54COF5@FFUNDS ///
	OA54OTL5@FFUNDS OA54MOR5@FFUNDS OA54COG5@FFUNDS OA54MFS5@FFUNDS OA54AFX3@FFUNDS ///
	OA54AFY3@FFUNDS OA54AFV3@FFUNDS OA54UNA5@FFUNDS OA54AHV3@FFUNDS OA51TAO5@FFUNDS ///
	OA66TAO5@FFUNDS OA66CCD3@FFUNDS OA66AHY3@FFUNDS OA66DPP3@FFUNDS OA66TRE5@FFUNDS ///
	OA66AGI3@FFUNDS OA66STO3@FFUNDS OA66COF5@FFUNDS OA66OTL5@FFUNDS OA66COG3@FFUNDS ///
	OA66EQG3@FFUNDS OA66MFT5@FFUNDS OA66TAO5@FFUNDS OA66CCD3@FFUNDS OA66AHY3@FFUNDS OA66DPP3@FFUNDS ///
	OA66TRE5@FFUNDS OA66AGI3@FFUNDS OA66STO3@FFUNDS OA66COF5@FFUNDS OA66OTL5@FFUNDS ///
	OA66COG3@FFUNDS OA66EQG3@FFUNDS OA66MFT5@FFUNDS OA64CCD5@FFUNDS OA64SAS3@FFUNDS ///
	OA64COF5@FFUNDS OA64MOR5@FFUNDS OA64MFT5@FFUNDS OA64TAO5@FFUNDS FLBZB@FFUNDS FLAZB@FFUNDS ///
	FLIBC@FFUNDS  oa47mbr5@ffunds oa47srp3@ffunds oa47slm3@ffunds oa47tre5@ffunds oa47agi5@ffunds ///
	oa47sto5@ffunds oa47cof5@ffunds oa47abn0@ffunds oa47mfs5@ffunds oa47mft5@ffunds oa47tao5@ffunds ///
	oa61ccd3@ffunds oa61tid3@ffunds oa61cof3@ffunds oa61flb0@ffunds oa61mor0@ffunds oa61cnc5@ffunds ///
	oa61eqg3@ffunds oa61mft5@ffunds oa61tao5@ffunds oa50mms3@ffunds oa50ahy3@ffunds oa50dpp5@ffunds ///
	oa50cof5@ffunds oa50cmb5@ffunds oa50cog5@ffunds ol75jum3@ffunds ol66jum5@ffunds oa50tao5@ffunds ///
	oa67tre3@ffunds oa67agi3@ffunds oa67otl5@ffunds oa67mor5@ffunds oa67cnc0@ffunds oa67trd3@ffunds ///
	oa67aiw3@ffunds oa67tao5@ffunds

*Insurance Companies
gen insurance_corp_bond = oa51cof3_ffunds + oa54cof5_ffunds
label var insurance_corp_bond "Corporate and Foreign Bonds"
gen insurance_corp_equities = oa51cog5_ffunds + oa54cog5_ffunds
label var insurance_corp_equities "Corporate Equities" 
gen insurance_mark_paper = oa51dpp3_ffunds + oa54dpp5_ffunds
label var insurance_mark_paper "Market Paper"
gen insurance_mutual = oa51mfs3_ffunds + oa54mfs5_ffunds
label var insurance_mutual "Mutual Fund Shares"
gen insurance_money_mark = oa51mms3_ffunds + oa54mms5_ffunds
label var insurance_money_mark "Money Market Mutual Fund Shares"
gen insurance_other_loans = oa54otl5_ffunds
label var insurance_other_loans "Other Loans and Advances"
gen insurance_other = oa51una5_ffunds + oa54una5_ffunds + oa51trd3_ffunds
label var insurance_other "Other/Unallocated Claims"
gen insurance_repos = oa51ahy3_ffunds + oa54ahv3_ffunds
label var insurance_repos "Security Repurchase Agreements"
gen insurance_checkable = oa51ccd3_ffunds + oa54ccd5_ffunds
label var insurance_checkable "Checkable Deposits and Currency"
gen insurance_agency = oa51agi3_ffunds + oa54agi5_ffunds
label var insurance_agency "Agency and GSE-backed securities"
gen insurance_munis = oa54sto5_ffunds + oa51sto3_ffunds
label var insurance_munis "Municipal Securities"
gen insurance_mortages = oa51cmt3_ffunds + oa54mor5_ffunds
label var insurance_mortages "Mortgages"
gen insurance_usdi = oa54afx3_ffunds + oa51eqg3_ffunds
label var insurance_usdi "US direct investment"
gen insurance_fhlb = oa51mfa3_ffunds + oa54afy3_ffunds
label var insurance_fhlb "Equity in FHLB"
gen insurance_premiums = oa54afv3_ffunds
label var insurance_premiums "Deferred and Unpaid Life ins Premiums"
gen insurance_treasuries = oa54tre5_ffunds + oa51tre3_ffunds
label var insurance_treasuries "Treasury Securities"

gen all_insurance = oa54tao5_ffunds + oa51tao5_ffunds

*Dealers
gen dealer_corp_bonds = oa66cof5_ffunds
label var dealer_corp_bonds "Corporate and Foreign Bonds"
gen dealer_equities = oa66cog3_ffunds
label var dealer_equities "Corporate Equities"
gen dealer_paper = oa66dpp3_ffunds
label var dealer_paper "Market Paper"
gen dealer_other = oa66mft5_ffunds
label var dealer_other "Other"
gen dealer_other_loans = oa66otl5_ffunds
label var dealer_other_loans "Other Loans"
gen dealer_repo = oa66ahy3_ffunds
label var dealer_repo "Security Repurchase Agreements"
gen dealer_checkable = oa66ccd3_ffunds
label var dealer_checkable "Checkable Deposits and Currency"
gen dealer_treas = oa66tre5_ffunds
label var dealer_treas "Treasury Securities"
gen dealer_agency = oa66agi3_ffunds
label var dealer_agency "Agency and GSE-Backed Securities"
gen dealer_muni = oa66sto3_ffunds
label var dealer_muni "Municipal Securities"
gen dealer_usdi = oa66eqg3_ffunds
label var dealer_usdi "US Direct Investment"

gen all_dealer = oa66tao5_ffunds

*REITs
gen reit_deposits = oa64ccd5_ffunds
label var reit_deposits "Checkable Deposits and Currency"
gen reit_agency = oa64sas3_ffunds
label var reit_agency "Agency and GSE-backed securities"
gen reit_corporates= oa64cof5_ffunds
label var reit_corporates "Corporate and Foreign Bonds"
gen reit_loans = oa64mor5_ffunds
label var reit_loans "Loans"
gen reit_other = oa64mft5_ffunds
label var reit_other "Other"
gen reit_nonfinancial = (flbzb_ffunds + flazb_ffunds)/1000 - oa64tao5_ffunds
label var reit_nonfinancial "Nonfinancial Assets"

gen all_reit = oa64tao5_ffunds + reit_nonfinancial

*Credit unions
gen credit_union_reserves = oa47mbr5_ffunds
label var credit_union_reserves "Reserves"
gen credit_union_ffunds = oa47srp3_ffunds
label var credit_union_ffunds "Fed Funds and Repos"
gen credit_union_paper = oa47slm3_ffunds
label var credit_union_paper "Open Market Paper"
gen credit_union_treas = oa47tre5_ffunds
label var credit_union_treas "Treasury Securities"
gen credit_union_agency = oa47agi5_ffunds
label var credit_union_agency "Agency and GSE-Backed"
gen credit_union_muni = oa47sto5_ffunds
label var credit_union_muni "Municipal Bonds"
gen credit_union_corp = oa47cof5_ffunds
label var credit_union_corp "Corporate and Foreign Bonds"
gen credit_union_loans = oa47abn0_ffunds
label var credit_union_loans "Loans"
gen credit_union_mut = oa47mfs5_ffunds
label var credit_union_mut "Mutual Fund Shares"
gen credit_union_other = oa47mft5_ffunds
label var credit_union_other "Other"

gen all_credit_union = oa47tao5_ffunds

*Finance Companies
gen finance_depos = oa61ccd3_ffunds
label var finance_depos "Checkable Deposits and Currency"
gen finance_time = oa61tid3_ffunds
label var finance_time "Time and Savings Deposits"
gen finance_corp = oa61cof3_ffunds
label var finance_corp "Corporate and Foreign Bonds"
gen finance_otherloan = oa61flb0_ffunds
label var finance_otherloan "Other Loans"
gen finance_mort = oa61mor0_ffunds
label var finance_mort "Mortgages"
gen finance_conscred = oa61cnc5_ffunds
label var finance_conscred "Consumer Credit"
gen finance_invest = oa61eqg3_ffunds
label var finance_invest "US Direct Investment Abroad"
gen finance_other = oa61mft5_ffunds
label var finance_other "Other"

gen all_finance = oa61tao5_ffunds

*Funding Corporations
gen funding_moneymark = oa50mms3_ffunds
label var funding_moneymark "Money Market Fund Shares"
gen funding_repo = oa50ahy3_ffunds
label var funding_repo "Security Repurchase Agreements"
gen funding_paper = oa50dpp5_ffunds
label var funding_paper "Open Market Paper"
gen funding_corp = oa50cof5_ffunds
label var funding_corp "Corporate and Foreign Bonds"
gen funding_otherloan = oa50cmb5_ffunds
label var funding_otherloan "Loans"
gen funding_equit = oa50cog5_ffunds
label var funding_equit "Corporate Equities"
gen funding_forbank= ol75jum3_ffunds
label var funding_forbank "Investment in Foreign Banks"
gen funding_brokers= ol66jum5_ffunds
label var funding_brokers "Investment in Brokers and Dealers"

gen all_funding = oa50tao5_ffunds

*Issuers of ABS
gen abs_treas= oa67tre3_ffunds
label var abs_treas "Treasury Securities"
gen abs_agency= oa67agi3_ffunds
label var abs_agency "Agency and GSE-Backed Securities"
gen abs_otherloans= oa67otl5_ffunds
label var abs_otherloans "Other Loans"
gen abs_mortgages= oa67mor5_ffunds
label var abs_mortgages "Mortgages"
gen abs_conscred= oa67cnc0_ffunds
label var abs_conscred "Consumer Credit"
gen abs_trade= oa67trd3_ffunds
label var abs_trade "Trade Credit"
gen abs_other= oa67aiw3_ffunds
label var abs_other "Other"

gen all_abs = oa67tao5_ffunds

*Percent breakdowns
foreach tag in insurance dealer reit credit_union finance funding abs{
	foreach var of varlist `tag'*{
		replace `var' = 100 * `var'/all_`tag'
	}
}

*Snapshot of last quarter in charted sample
keep if time == $charts_end

foreach v of varlist *{
	capture local l`v': variable label `v'
}
keep insurance* dealer* reit* credit_union* finance* funding* abs*
xpose, clear varname
levelsof _varname, clean local(vars)
gen component = ""
foreach v in `vars'{
	replace component = "`l`v''" if _varname == "`v'"
}

gsort -v1
gen sorter = _n
labmask sorter, values(component)

*Formatted the same way as sum stats in sum_stats_wholesample.do
foreach tag in insurance dealer reit credit_union finance funding abs{
	estpost tabstat v1 if strpos(_varname, "`tag'"), columns(statistics) by(sorter) nototal
	esttab using ../output/asset_breakdown_`tag'.tex, coeflabels(`e(labels)') ///
		main(mean 1) booktabs nonum replace title(" ") mtitle("This is a placeholder! (%)") ///
		noobs nonotes
}

*********************************************************************************
*
*					Pulling total financial assets from sectors
*					Calculating %in, %out for those sectors that were categorized						
*
*********************************************************************************

*Pulling total financial assets
clear
import haver OA51TAO5@FFUNDS OA54TAO5@FFUNDS OA65TAO0@FFUNDS OA70TAO5@FFUNDS OA76TAO5@FFUNDS OA73TAO5@FFUNDS OA66TAO5@FFUNDS ///
	OA64TAO5@FFUNDS OA79TAO5@FFUNDS OA63TAO5@FFUNDS OA71TAO5@FFUNDS OA59TAO5@FFUNDS OA40TAO5@FFUNDS OA41MOR5@FFUNDS FLBZB@FFUNDS FLAZB@FFUNDS ///
	OA56TAO5@FFUNDS	OA55TAO5@FFUNDS OA75TAO5@FFUNDS OA50TAO5@FFUNDS OA61TAO5@FFUNDS OA74TAO5@FFUNDS OA47TAO5@FFUNDS OA67TAO5@FFUNDS ///
	OL66TAO5@FFUNDS OL51TAO5@FFUNDS OL54TAO5@FFUNDS OL64TAO5@FFUNDS OL67AIJ5@FFUNDS OL47TAO5@FFUNDS OL50TAO5@FFUNDS OL61TAO5@FFUNDS FLIBC@FFUNDS

rename time qt_dt
rename oa51tao5_ffunds total_property
rename oa54tao5_ffunds total_life
rename oa65tao0_ffunds total_mutual
rename oa70tao5_ffunds total_depository
rename oa76tao5_ffunds total_us_depository
rename oa73tao5_ffunds total_holding_companies
rename oa66tao5_ffunds total_dealer
rename oa79tao5_ffunds total_domestic_financial
rename oa63tao5_ffunds total_money_market_mutual
rename oa71tao5_ffunds total_monetary_authority
rename oa59tao5_ffunds total_pension_funds
rename oa40tao5_ffunds total_gse
rename oa41mor5_ffunds total_agency_mortgage
rename oa56tao5_ffunds total_etf
rename oa55tao5_ffunds total_cef
rename oa75tao5_ffunds total_foreign_banks
rename oa50tao5_ffunds total_fundingcorp
rename oa61tao5_ffunds total_financecorp
rename oa74tao5_ffunds total_bank_affilitied
rename oa47tao5_ffunds total_credit_unions
rename oa67tao5_ffunds total_abs_issuers
rename ol66tao5_ffunds liab_total_dealer
rename ol51tao5_ffunds liab_total_property
rename ol54tao5_ffunds liab_total_life
gen liab_total_insurance = liab_total_property + liab_total_life
rename ol64tao5_ffunds liab_total_reit
gen liab_total_other = ol67aij5_ffunds + ol47tao5_ffunds + ol50tao5_ffunds + ol61tao5_ffunds
gen total_reit = (flbzb_ffunds + flazb_ffunds)/1000
gen total_financial_ours = total_domestic_financial + (total_reit - oa64tao5_ffunds) - total_gse - total_agency_mortgage - ///
	total_pension_funds - total_monetary_authority - total_mutual - total_money_market_mutual - total_etf - total_cef - total_foreign_banks - total_bank_affilitied

*Merge in categorizations
merge 1:1 qt_dt using ../temp/ffunds.dta, nogen
keep if qt_dt >= yq(2002, 1)
gen total_insurance = total_life + total_property // + total_liberty
gen total_other = total_financial_ours - total_reit - total_us_depository - total_holding_companies - total_dealer - total_insurance

*As a test - it all checks out
//gen total_other_test = total_abs_issuers+total_credit_unions+total_fundingcorp+total_financecorp

foreach tag in insurance dealer reit other{
	gen total_`tag'_orig = total_`tag'
	gen total_`tag'_full = total_`tag'
}

keep qt_dt insurance* mutual* dealer* total* insurance* reit* other* liab_total*

*Temporary output
save ../temp/ffunds.dta, replace
