*********************************************************************************
*
*	Summary Statistics for Y9C Sample
*	Collin Jones 
*	Purpose: Construct some summary statistics tables of BHC assets and liabilities
*
*********************************************************************************


clear
import excel ../input/categorizations.xlsx, sheet("Assets") firstrow
save ../temp/categorizations_assets, replace
clear
import excel ../input/categorizations.xlsx, sheet("Liabilities") firstrow
save ../temp/categorizations_liabilities, replace

use ../temp/categorizations_assets, clear
append using ../temp/categorizations_liabilities
save ../temp/categorizations, replace


global bankSample 1039502 1068025 1069778 1070345 1073757 1074156 1111435 1119794 ///
	1120754 1131787 1275216 1562859 1951350 2162966 2277860 2380443 3242838 3587146 3232316 1132449

program drop _all
program format_transposed
	foreach v of varlist _var*{
		local label_curr = `v'[1]
		local rename_curr = "`v'_" + `v'[2]
		label var `v' "`label_curr'"
		rename `v' `rename_curr'
		if "`label_curr'" == "Other"{
			order `rename_curr', last
			local helper = "other_" + `rename_curr'[2]
			rename `rename_curr' `helper'
		}
	}

end

program other_compile
	drop if _n == 1 | _n == 2
	destring other*, replace
	destring _var*, replace

	foreach bsheet in Assets Liabilities{
		foreach v of varlist *_`bsheet'{
			if `v'[_N] < 1{
				replace other_`bsheet' = other_`bsheet' + `v'
				drop `v'
				}
			}
		}
end

program gen_dates
	gen qt_dt = yq(2002 + int((_n-1)/4) * (_n>4), 1 + mod(_n-1, 4))
	format qt_dt %tq
	order qt_dt, first
	destring(_var*), replace
end 

clear

use ../temp/Contagion_preprocessing

gen keep_obs = 0

foreach bank of numlist $bankSample {
	replace keep_obs = 1 if entity == `bank'
}

replace BHCKG320 = assets_htm_mbs if qt_dt <yq(2009, 2)

drop BHCKG336 BHCKG340 BHCKG344
replace BHCKC026 = assets_in_htm_abs

replace BHCKG323 = assets_in_afs_mbs

drop BHCKG339 BHCKG343 BHCKG347
replace BHCKC027 = assets_in_afs_abs

drop BHCKG384 BHCKG385
replace BHCKG383 = assets_in_oth_debt

replace BHCKG319 = assets_in_afs_cmo_unc
replace BHCKG316 = assets_in_htm_cmo_unc

drop BHCKK146 BHCKK154
replace BHCKG308 = assets_in_htm_mbs_unc

drop BHCKK149 BHCKK157
replace BHCKG311 = assets_in_afs_mbs_unc

replace BHCM3541 = assets_in_trad_asst_unc
replace BHCM3543 = assets_in_deriv_unc

drop BHCKK198
replace BHCKG381 = assets_in_trad_mbs_unc

replace BHCKG209 = BHCK3546 if qt_dt < yq(2009, 1)

replace BHCKG300 = BHCKG300 + BHCKG304 + BHCKK142
drop BHCKG304 BHCKK142
replace BHCKG300 = BHCK1698 + BHCK1703 if qt_dt <yq(2009, 2)

replace BHCKG312 = BHCKG312 + BHCKK150
drop BHCKK150
replace BHCKG312 = BHCK1714+BHCK1718 - assets_in_htm_cmo_unc if qt_dt < yq(2009, 2)

replace BHCKG303 = BHCKG303 + BHCKG307 + BHCKK145
drop BHCKG307 BHCKK145
replace BHCKG303 = BHCK1702 + BHCK1707 if qt_dt < yq(2009, 2)

replace BHCKG315 = BHCKG315 + BHCKK153
drop BHCKK153
replace BHCKG315 = BHCK1717 + BHCK1732 - assets_in_afs_cmo_unc if qt_dt < yq(2009, 2)

replace BHCKG379 = BHCKG379 + BHCKG380 + BHCKK197
drop BHCKG380 BHCKK197
replace BHCKG379 = BHCK3534 + BHCK3535 if qt_dt < yq(2009, 2)

replace BHCKG386 = BHCK3537 - assets_in_oth_debt if qt_dt < yq(2009, 2)

replace BHCM3531 = BHCK3531 if qt_dt<yq(2008, 1)
replace BHCM3532 = BHCK3532 if qt_dt<yq(2008, 1)
replace BHCM3533 = BHCK3533 if qt_dt<yq(2008, 1)

keep if keep_obs == 1
drop keep_obs
//qui:tostring *, force replace

ds, has(type numeric)
keep `r(varlist)' qt_dt
order entity, first
order qt_dt, first
//collapse (sum) row_loaded_tstmp- liabs_in_frac, by(qt_dt)
xpose, clear varname
//sxpose, clear firstnames
rename _varname tag
foreach v of varlist v*{
	if `v'[2] == 1039502{
		local qt = `v'[1]
		rename `v' JPMC_`qt'
	}
	else if `v'[2] == 1068025{
		local qt = `v'[1]
		rename `v' KeyCorp_`qt'
	}
	else if `v'[2] == 1069778{
		local qt = `v'[1]
		rename `v' PNC_`qt'
	}
	else if `v'[2] == 1070345{
		local qt = `v'[1]
		rename `v' Fifth_Third_`qt'
	}
	else if `v'[2] == 1073757{
		local qt = `v'[1]
		rename `v' BofA_`qt'
	}
	else if `v'[2] == 1111435{
		local qt = `v'[1]
		rename `v' State_Street_`qt'
	}
	else if `v'[2] == 1119794{
		local qt = `v'[1]
		rename `v' US_Bancorp_`qt'
	}
	else if `v'[2] == 1120754{
		local qt = `v'[1]
		rename `v' Wells_Fargo_`qt'
	}
	else if `v'[2] == 1131787{
		local qt = `v'[1]
		rename `v' Suntrust_`qt'
	}
	else if `v'[2] == 1275216{
		local qt = `v'[1]
		rename `v' American_Express_`qt'
	}
	else if `v'[2] == 1562859{
		local qt = `v'[1]
		rename `v' Ally_`qt'
	}
	else if `v'[2] == 1951350{
		local qt = `v'[1]
		rename `v' Citigroup_`qt'
	}
	else if `v'[2] == 2162966{
		local qt = `v'[1]
		rename `v' Morgan_Stanley_`qt'
	}
	else if `v'[2] == 2277860{
		local qt = `v'[1]
		rename `v' Capital_One_`qt'
	}	
	else if `v'[2] == 2380443{
		local qt = `v'[1]
		rename `v' Goldman_Sachs_`qt'
	}
	else if `v'[2] == 3242838{
		local qt = `v'[1]
		rename `v' Regions_`qt'
	}
	else if `v'[2] == 3587146{
		local qt = `v'[1]
		rename `v' BONY_Mellon_`qt'
	}
	else if `v'[2] == 3232316{
		local qt = `v'[1]
		rename `v' HSBC_North_America_`qt'
	}
	else if `v'[2] == 1132449{
		local qt = `v'[1]
		rename `v' RBC_USA_`qt'
	}
	else if `v'[2] == 1074156{
		local qt = `v'[1]
		rename `v' BBandT_`qt'
	}
}

/*
JPMC KeyCorp PNC Fifth_Third BofA BBandT State_Street US_Bancorp Wells_Fargo Suntrust
American_Express Ally Citigroup Morgan_Stanley Capital_One Goldman_Sachs Regions BONY_Mellon
HSBC_North_America RBC_USA
*/

merge 1:1 tag using ../temp/categorizations
keep if _merge == 3 | tag == "BHCK2170" | tag == "BHCK2948"

replace BalanceSheet = "Assets" if tag == "BHCK2170"
replace Component = "ZZZ All Assets" if tag == "BHCK2170"

replace BalanceSheet = "Liabilities" if tag == "BHCK2948"
replace Component = "ZZZ All Liabilities" if tag == "BHCK2948"

global start_date = yq(2002, 1) 
global end_date = yq(2016, 3) 

preserve

foreach x of numlist $start_date(1)$end_date{
	local to_name = "y" + string(year(dofq(`x'))) + "q" + string(quarter(dofq(`x')))
	egen `to_name' = rowtotal(*_`x')
}
keep y*q* tag BalanceSheet Component

collapse (sum) y*q*, by (BalanceSheet Component)

sort BalanceSheet Component
foreach x of varlist y20*{
	bys BalanceSheet: replace `x' = 100 * `x'/`x'[_N]
}
bys BalanceSheet (Component): drop if _n == _N

egen avg = rowmean(y*q*)
gsort - avg
drop avg

tostring y*, replace force

sxpose, clear

format_transposed
other_compile
gen_dates

estpost sum *_Assets
esttab using ../output/sum_stat_all.tex, cell("mean(fmt(2)) sd(fmt(2)) min(fmt(2)) max(fmt(2))") replace title("Assets: Whole Sample") tex width("13cm") label
estpost sum *_Liabilities
esttab using ../output/sum_stat_all.tex, cell("mean(fmt(2)) sd(fmt(2)) min(fmt(2)) max(fmt(2))") append title("Liabilities: Whole Sample'") tex width("13cm") label

//foreach bsheet in Assets Liabilities{
//	sutex *_`bsheet', labels file("sum_stat_test_`bsheet'") replace minmax digits(2)
//}
restore

local banknames JPMC KeyCorp PNC Fifth_Third BofA BBandT State_Street US_Bancorp Wells_Fargo Suntrust ///
	American_Express Ally Citigroup Morgan_Stanley Capital_One Goldman_Sachs Regions BONY_Mellon ///
	HSBC_North_America RBC_USA
disp "`banknames'"

foreach bank in `banknames'{
	preserve
	keep `bank'* tag BalanceSheet Component
	collapse (sum) `bank'*, by (BalanceSheet Component)
	sort BalanceSheet Component
	foreach v of varlist `bank'*{
		bys BalanceSheet: replace `v' = 100 * `v'/`v'[_N]
	}
	bys BalanceSheet (Component): drop if _n == _N
	egen avg = rowmean(`bank'*)
	gsort - avg
	drop avg
	tostring `bank'*, replace force
	sxpose, clear
	format_transposed
	other_compile
	gen_dates
	local bname_tex = subinstr("`bank'", "_", " ", 100)
	estpost sum *_Assets
	esttab using ../output/`bank'.tex, cell("mean(fmt(2)) sd(fmt(2)) min(fmt(2)) max(fmt(2))") replace title("Assets: `bname_tex'") tex width(13cm) label 
	estpost sum *_Liabilities
	esttab using ../output/`bank'.tex, cell("mean(fmt(2)) sd(fmt(2)) min(fmt(2)) max(fmt(2))") append title("Liabilities: `bname_tex'") tex width(13cm) label
	restore
}

//	Some cross-sectional moment business

foreach val in sd mean{
	preserve
	drop _merge
	collapse (sum) *_*, by (BalanceSheet Component)
	sort BalanceSheet Component
	foreach x of varlist *_*{
		bys BalanceSheet: replace `x' = 100 * `x'/`x'[_N]
	}
	bys BalanceSheet (Component): drop if _n == _N
	egen avg = rowmean(*_*)
	gsort - avg
	drop avg

	foreach x of numlist $start_date(1)$end_date{
		egen sd_`x' = row`val'(*_`x')
	}

	tostring sd* , replace force
	keep sd* Component BalanceSheet
	sxpose, clear
	format_transposed
	drop if _n == 1 | _n == 2
	destring *, replace

	estpost sum *_Assets
	esttab using ../output/assets_cx_`val'.tex, cell("mean(fmt(2))") replace title(" ") width("13cm") label booktabs nonum noobs
	estpost sum *_Liabilities
	esttab using ../output/liabilities_cx_`val'.tex, cell("mean(fmt(2))") replace title(" ") width("13cm") label booktabs nonum noobs
	restore
	}

preserve
drop _merge
sort BalanceSheet Component
collapse (sum) *_*, by (BalanceSheet Component)

foreach x of numlist $start_date(1)$end_date{
	egen all_`x' = rowtotal(*_`x')
}

foreach x of varlist *_*{
	bys BalanceSheet: replace `x' = 100 * `x'/`x'[_N]
}
bys BalanceSheet: drop if _n == _N
keep all_* Component BalanceSheet
egen avg = rowmean(all_*)
gsort - avg
drop avg
tostring all_*, replace force
sxpose, clear
format_transposed
drop if _n == 1 | _n == 2
destring *, replace
estpost sum *_Assets
esttab using ../output/assets_ts.tex, ///
	cell("mean(fmt(2)) sd(fmt(2))") replace title(" ") width("13cm") label booktabs nonum noobs
estpost sum *_Liabilities
esttab using ../output/liabilities_ts.tex, ///
	cell("mean(fmt(2)) sd(fmt(2))") replace title(" ") width("13cm") label booktabs nonum noobs
restore

preserve
keep *_$end_date Component InOutUnclear BalanceSheet tag
keep if Component != "ZZZ All Assets" & Component != "ZZZ All Liabilities"
//expand 2 if InOutUnclear != "In" & InOutUnclear != "Out"
//expand 2 if InOutUnclear != "In" & InOutUnclear != "Out"
expand 2
egen val_all = rowtotal(*_$end_date)
bys tag: gen n = _n 
replace val_all = .5 * val_all if InOutUnclear != "In" & InOutUnclear != "Out"
replace val_all = 0 if InOutUnclear == "In" & n == 2
replace val_all = 0 if InOutUnclear == "Out" & n == 1
gen In_Out = InOutUnclear
replace In_Out = "In" if n == 1
replace In_Out = "Out" if n == 2
keep val_all Component BalanceSheet In_Out
collapse (sum) val_all, by(Component BalanceSheet In_Out)

bys BalanceSheet: egen tot_side = sum(val_all)
bys BalanceSheet In_Out: egen tot_inout = sum(val_all)

gen perc_side = 100 * val_all/tot_side
replace Component = "Other" if perc_side <= 1
gen perc_inout = 100 * val_all/tot_inout

bys BalanceSheet In_Out: gen n = _n


collapse (sum) val_all perc_side perc_inout (firstnm) tot_side (firstnm) tot_inout , by(Component BalanceSheet In_Out)

gsort -perc_inout
gen sorter = _n
replace sorter = 100 if Component == "Other"
labmask sorter, values(Component)

foreach side in Assets Liabilities{


	estpost tabstat perc_inout if BalanceSheet == "`side'" & In_Out == "Out", by(sorter) nototal 
	esttab using ../output/`side'_Out.tex, main(mean 2) coeflabels(`e(labels)') booktabs nonum label replace title(" ") mtitle("`side': Out (%)") noobs nonotes
	
	estpost tabstat perc_inout if BalanceSheet == "`side'" & In_Out == "In", by(sorter) nototal 
	esttab using ../output/`side'_In.tex, main(mean 2) coeflabels(`e(labels)') booktabs nonum label replace title(" ") mtitle("`side': In (%)") noobs nonotes
	}

collapse (sum) val_all, by (BalanceSheet In_Out)
bys BalanceSheet: egen tot_side = sum(val_all)
gen perc_inout = 100*val_all/tot_side
restore
