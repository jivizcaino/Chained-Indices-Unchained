version 16.0
capture log close
clear all
set more off
global DataPath `"/Users/akos/Dropbox/AkosBertRichard/Investment/Submissions/FinalSubmission/Supplementary"'
cd "$DataPath/data"
import excel "../rawdata/naics.xlsx", sheet("Sheet1") firstrow clear
rename desc  description
tempfile naics
save `naics'
*******
***
*** Consumption and investment value added requirement from the IO-tables
***               broken down by goods and services
***
******
use "Data-IO-FEValueAdded1947-1962.dta", clear
merge m:1 naics using `naics'
generate _sector="GOOD" if Industry<34 & _merge==3
replace  _sector="SERV" if Industry>33 & _merge==3
generate VAX=VAPIE+VAPIS+VAPIP+VAGIE+VAGIS+VAGIP
generate VAC=VAGDP-VAX
collapse (sum) VAC VAX, by(_sector year) 
tempfile d1
save `d1'
use "Data-IO-FEValueAdded1963-1996.dta", clear
merge m:1 naics using `naics'
generate _sector="GOOD" if Industry<34 & _merge==3
replace  _sector="SERV" if Industry>33 & _merge==3
generate VAX=VAPIE+VAPIS+VAPIP+VAGIE+VAGIS+VAGIP
generate VAC=VAGDP-VAX
collapse (sum) VAC VAX, by(_sector year) 
tempfile d2
save `d2'
use "Data-IO-FEValueAdded1997-2017.dta", clear
merge m:1 naics using `naics'
generate _sector="GOOD" if Industry<34 & _merge==3
replace  _sector="SERV" if Industry>33 & _merge==3
generate VAX=VAPIE+VAPIS+VAPIP+VAGIE+VAGIS+VAGIP
generate VAC=VAGDP-VAX
collapse (sum) VAC VAX, by(_sector year) 
tempfile d3
save `d3'
use `d1', clear
append using `d2'
append using `d3'
sort year _sector
foreach sec in GOOD SERV {
	foreach var in C X  {
		by year: egen _VA`var'_`sec'=sum(VA`var'*(_sector=="`sec'"))
	}
}
sort  _sector year
keep if _sector=="GOOD"
drop _sector  VAC VAX
tempfile vagdp
save `vagdp'
*******
***
*** WORLD KLEMS value added data and prices for goods and services
***              NIPA consumption and investment
***
******
merge 1:1 year using "Data-WORLDKLEMS.dta", generate(_merge_1)
merge 1:1 year using "Data-NIPA.dta", generate(_merge_2)

generate VAX_GOOD=round(X_TOT*_VAX_GOOD/(_VAX_GOOD+_VAX_SERV),1)
generate VAC_GOOD=round((VA_TOT-X_TOT)*(_VAC_GOOD)/(_VAC_GOOD+_VAC_SERV),1)
generate VAX_SERV=X_TOT-VAX_GOOD
generate VAC_SERV=VA_TOT-X_TOT-VAC_GOOD
rename C_TOT _C_TOT
generate  C_TOT=VA_TOT-X_TOT
replace C_GOOD=round(C_TOT*C_GOOD/_C_TOT,1)
replace C_SERV=C_TOT-C_GOOD


drop _* 

order POP GDP GDP_QI GDP_P PCE PIE GCE GIE NX VA_TOT VA_GOOD VA_SERV VA_TOT_P VA_GOOD_P VA_SERV_P VAX_GOOD VAC_GOOD VAX_SERV VAC_SERV LAB_TOT LAB_GOOD LAB_SERV LAB_S_TOT VA_TOT_QI VA_GOOD_QI VA_SERV_QI ///
      CAP_TOT_QI CAP_GOOD_QI CAP_SERV_QI LAB_TOT_QI LAB_GOOD_QI LAB_SERV_QI LAB_S_GOOD LAB_S_SERV  C_TOT ///
	  C_GOOD C_SERV C_TOT_QI C_TOT_P C_GOOD_QI C_GOOD_P C_SERV_QI C_SERV_P X_TOT X_TOT_QI X_TOT_P , after(year)

label variable C_TOT       "Consumption defined as GDP-X, current prices, millions of dollars"


label variable VAC_GOOD    "Consumption value added requirement, services, millions of dollars"
label variable VAC_SERV    "Consumption value added requirement, services, millions of dollars"

label variable VAX_GOOD    "Investment value added requirement, services, millions of dollars"
label variable VAX_SERV    "Investment value added requirement, services, millions of dollars"

tsset year, yearly

save "../paper/Data-MainDatabase.dta", replace


exit
