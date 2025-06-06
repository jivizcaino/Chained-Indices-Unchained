version 16.0
capture log close
clear all
set more off
global DataPath `"/Users/akos/Dropbox/AkosBertRichard/Investment/Submissions/FinalSubmission/Supplementary"'
cd "$DataPath/data"
*****
**
** Extending WORLD KLEMS data for the period 2015-2017 
**
*****
** Industry names and codes
*****
import excel "../rawdata/naics.xlsx", sheet("Sheet1") firstrow clear
drop if Industry==78
replace Industry=_n
replace naics=stritrim(naics)
tempfile naics
save `naics'
*****
** Industry value added data from BEA
*****
use "../rawdata/GDPbyInd_VA_1947-2017.dta", clear
keep if naics=="GF" |naics=="GFG" | naics=="GFE" | naics=="GSL" | naics=="GSLG" | naics=="GSLE"
keep if year>2013 & year<2018
sort year Industry
by year: egen VA_GF=sum(VA*(naics=="GF"))
by year: egen VA_GSL=sum(VA*(naics=="GSL"))
by year: egen VA_QI_GF=sum(VA_QI*(naics=="GF"))
by year: egen VA_QI_GSL=sum(VA_QI*(naics=="GSL"))
sort Industry year
by Industry: generate  _dlVA_QI_GF=ln(VA_QI)-ln(VA_QI[1])
by Industry: generate _dlVA_QI_GSL=ln(VA_QI)-ln(VA_QI[1])
generate _VA_GF_S=VA/VA_GF if naics=="GFE" | naics=="GFG"
generate _VA_GSL_S=VA/VA_GSL if naics=="GSLE" | naics=="GSLG"
drop if naics=="GF" | naics=="GSL"
keep naics year _dlVA_QI_GF _dlVA_QI_GSL _VA_GF_S _VA_GSL_S
tempfile gov
save `gov'

*****
** BEA Integrated Industry-Level Production Account (KLEMS)
**    expanding the data set from 63 to 65 industries 
*****

use "../rawdata/BEA-BLS-industry-level-production-account-1998-2017.dta", clear
drop GO GO_QI II II_QI LABCOL LABNOCOL CAPART CAPIT CAPO CAPRD CAPSW IIE IIM IIS LABCOL_QI LABNOCOL_QI  ///
     CAPART_QI CAPIT_QI CAPO_QI CAPRD_QI CAPSW_QI IIE_QI IIM_QI IIS_QI LPGO_I H_EMP_I TFPGO_I 
keep if year>2013 
local ngf=_N
qui sum year
local ngl=_N+(`r(max)'-2013)
expand 2 if Industry==62
expand 2 if Industry==63
replace Industry=64 if desc=="State and local government"
replace Industry=65 if desc=="State and local government" & _n>`ngl'
replace Industry=63 if desc=="Federal government" & _n>`ngf'
sort Industry year
replace naics="GFG"  if Industry==62
replace naics="GFE"  if Industry==63
replace naics="GSLG" if Industry==64
replace naics="GSLE" if Industry==65
replace desc="Federal general government"  if Industry==62
replace desc="Federal government enterprises"  if Industry==63
replace desc="State and local general government" if Industry==64
replace desc="State and local government enterprises" if Industry==65
sort Industry year
foreach var in VA VA_QI LAB_QI CAP_QI {
	by Industry: generate _dl`var'=ln(`var')-ln(`var'[1])
}
generate _LAB_S=LAB/VA
rename VA _VA
keep naics year  _dlVA _dlVA_QI _dlLAB_QI _dlCAP_QI _LAB_S _VA
tempfile beabls
save `beabls'

*****
** 
** Extending the WORLD KLEMS data for 2015-2017
** 
*****

use "../rawdata/usa_wk_mar_2017.dta", clear
drop GO II GO_Q II_Q

merge 1:1 year naics using `beabls', generate(_merge1)
merge 1:1 year naics using `gov', generate(_merge2)

sort naics year
replace Industry=Industry[_n-1] if Industry[_n]==. 
replace desc=desc[_n-1] if desc[_n]==""
sort Industry year
****
** Nominal variables
**** 
by Industry: egen _VA_2014=sum(VA*(year==2014))
sort year Industry 
by year: egen  _VA_GF_2014=sum(_VA_2014*(substr(naics,1,2)=="GF"))
by year: egen _VA_GSL_2014=sum(_VA_2014*(substr(naics,1,3)=="GSL"))
sort Industry year
replace VA=_VA_GF_2014*exp(_dlVA) *_VA_GF_S  if substr(naics,1,2)=="GF"  & year>2014
replace VA=_VA_GSL_2014*exp(_dlVA)*_VA_GSL_S if substr(naics,1,3)=="GSL" & year>2014
replace VA=_VA_2014*exp(_dlVA) if year>2014 & Industry<=61
replace LAB=_LAB_S*VA if year>2014
replace CAP=VA-LAB if year>2014
**** 
** Real variables
**** 
foreach var in VA LAB CAP {
	by Industry: egen _`var'_Q_2014=sum(`var'_Q*(year==2014))
	replace `var'_Q=_`var'_Q_2014*exp(_dl`var'_QI) if year>2014
}
foreach var in GF GSL {
	replace VA_Q=_VA_Q_2014*exp(_dlVA_QI_`var') if year>2014 & substr(naics,1,2)=="`var'"
}

**
** Quantity indexes 1947=1, and quantities at 1947 dollars
**

foreach var in VA LAB CAP {
	by Industry: generate _dl`var'_Q=ln(`var'_Q)-ln(L1.`var'_Q)
	by Industry: generate `var'_QI=exp(sum(_dl`var'_Q))
	by Industry: replace `var'_Q=`var'_QI*`var'[1]
}
drop _*
tempfile wk
save `wk'
******
***
*** Aggregating the economy into total economy, goods and service sectors
***
***
******
***
*** Appending space for the the three aggregates
***
******
drop _all
set obs 71
generate year=1946+_n
tempfile year
save `year'
use `wk', clear
sort Industry year
local i=65
foreach var in TOT GOOD SERV {
    local i=`i'+1
	append using `year'
	replace naics="`var'" if naics==""
	replace Industry=`i' if Industry==.
}
******
***
*** Defining the three aggregatestes 
***
******
generate Sector1="TOT"  if Industry<=65
generate Sector2="GOOD" if Industry<=65
replace  Sector2="SERV" if Industry>=27 & Industry<=65
replace desc="Total economy" if Industry==66
replace desc="Goods" if Industry==67
replace desc="Services" if Industry==68
order Sector1 Sector2, after(desc)
******
***
*** Total economy, goods and services at current prices
***
******
tsset Industry year, yearly
sort year Industry 
foreach var in VA LAB CAP  {
	by year:  egen _`var'_TOT=sum(`var'*(Sector1=="TOT"))
	replace  `var'=_`var'_TOT if naics=="TOT"
}
foreach sec in GOOD SERV  {
	foreach var in VA LAB CAP  {
		by year:  egen _`var'_`sec'=sum(`var'*(Sector2=="`sec'"))
		replace  `var'=_`var'_`sec' if naics=="`sec'"
	}
}
******
***
*** Industry shares in total economy averaged over t and t-1
***
******
sort Industry year 
foreach var in  VA LAB CAP {
	generate _`var'_S_TOT=0.5*(`var'/_`var'_TOT+L1.`var'/L1._`var'_TOT) if year>1947
}
******
***
*** Sector shares in total economy averaged over t and t-1
***
******
sort year Industry 
foreach sec in GOOD SERV  {
	foreach var in  VA LAB CAP {
		by year: egen _`var'_S_`sec'=sum(_`var'_S_TOT*(Sector2=="`sec'")) if year>1947
	}
}
******
**
** Industry quantity changes 
**
******
sort Industry year 
foreach var in  VA LAB CAP {
    generate _dl`var'_QI=ln(`var'_QI)-ln(L1.`var'_QI) 
}
******
***
*** Total economy aggregate of industry quantity changes
***
sort year Industry 
foreach var in  VA LAB CAP {
	by year: egen _dl`var'_QI_TOT=sum(_`var'_S_TOT*_dl`var'_QI*(Sector1=="TOT"))  if year>1947
}
******
***
***  Good and service aggregates of industry quantity changes
***
******
foreach sec in GOOD SERV  {
	foreach var in  VA LAB CAP {
		by year: egen _dl`var'_QI_`sec'=sum(_`var'_S_TOT*_dl`var'_QI/_`var'_S_`sec'*(Sector2=="`sec'"))  if year>1947
	}
}
******
***
***  Income shares 
***
******
sort Industry year 
generate LAB_S=0.5*(LAB/VA+L1.LAB/L1.VA)
******
***
***  Total economy, goods and service quantity indexes 1947=1
***
******
foreach sec in  TOT GOOD SERV  {
	foreach var in  VA LAB CAP  {
		by Industry: generate _d`var'_QI_`sec'=exp(sum(_dl`var'_QI_`sec'))
	}
}
foreach sec in  TOT GOOD SERV  {
	foreach var in  VA LAB CAP  {
		replace `var'_QI=_d`var'_QI_`sec' if naics=="`sec'"
	}
}
drop _*
******
***
***  Keeping total economy, goods and service nominal values and quantities in 1947 dollars
***
******
keep if naics=="TOT" | naics=="GOOD" | naics=="SERV"
drop Industry desc Sector1-Sector2 
tempfile wk
save `wk', replace
local wklist  "VA CAP LAB VA_QI CAP_QI LAB_QI LAB_S"
foreach var in `wklist' {
    use `wk', clear
	keep naics year `var' 
	reshape wide `var', i(year) j(naics, string)
	ds, has(format *g) 
	foreach subvar in `r(varlist)' {
		local sec = substr("`subvar'",length("`var'")+1,length("`subvar'"))
		rename `subvar' `var'_`sec'
	}
	tempfile `var'file
	save ``var'file'
}
use `VAfile', clear
foreach var in CAP LAB VA_QI CAP_QI LAB_QI LAB_S {
	merge 1:1 year using ``var'file', generate(_merge_`var')
}

foreach var in VA CAP LAB  {
	foreach sec in TOT GOOD SERV {
		rename  `var'_QI_`sec' `var'_`sec'_QI
	}
}
foreach sec in TOT GOOD SERV {
	generate VA_`sec'_P=VA_`sec'/VA_`sec'[1]/VA_`sec'_QI
}

drop _*

drop  CAP_GOOD CAP_SERV CAP_TOT


order VA_TOT VA_GOOD VA_GOOD VA_SERV VA_TOT_P VA_GOOD_P VA_SERV_P LAB_TOT LAB_GOOD LAB_SERV VA_TOT_QI VA_GOOD_QI VA_SERV_QI CAP_TOT_QI ///
      CAP_SERV_QI CAP_GOOD_QI LAB_TOT_QI LAB_GOOD_QI LAB_SERV_QI LAB_S_TOT LAB_S_GOOD LAB_S_SERV, after(year)


label variable VA_TOT       "GDP, current prices, millions of dollars"
label variable VA_GOOD      "Value added, current prices, goods, millions of dollars"
label variable VA_SERV      "Value added, current prices, services, millions of dollars"

label variable LAB_TOT      "Labour income, current prices, millions of dollars"
label variable LAB_GOOD     "Labour income, current prices, goods, millions of dollars"
label variable LAB_SERV     "Labour income, current prices, services, millions of dollars"

label variable VA_TOT_QI     "Value added, volume index, aggregate, 1947 = 1"
label variable VA_GOOD_QI    "Value added, volume index, goods, 1947 = 1"
label variable VA_SERV_QI    "Value added, volume index, services, 1947 = 1"

label variable VA_TOT_P    "Value added, price index, aggregate, 1947 = 1"
label variable VA_GOOD_P   "Value added, price index, goods, 1947 = 1"
label variable VA_SERV_P   "Value added, price index, services, 1947 = 1"

label variable LAB_TOT_QI  "Efficiency hours worked by persons engaged, volume index, aggregate, 1947 = 1"
label variable LAB_GOOD_QI "Efficiency hours worked by persons engaged, volume index, goods, 1947 = 1"
label variable LAB_SERV_QI "Efficiency hours worked by persons engaged, volume index, services, 1947 = 1"

label variable LAB_S_TOT   "Labour income share in value added, aggregate, t and t-1 average"
label variable LAB_S_GOOD  "Labour income share in value added, goods, t and t-1 average"
label variable LAB_S_SERV  "Labour income share in value added, services, t and t-1 average"

label variable CAP_TOT_QI  "Capital services, volume index, aggregate, 1947 = 1"
label variable CAP_GOOD_QI "Capital services, volume index, goods, 1947 = 1"
label variable CAP_SERV_QI "Capital services, volume index, services, 1947 = 1"

save "Data-WORLDKLEMS.dta", replace

exit




