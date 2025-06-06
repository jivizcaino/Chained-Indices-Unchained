version 16.0
capture log close
clear all
set more off
global DataPath `"/Users/akos/Dropbox/AkosBertRichard/Investment/Submissions/FinalSubmission/Supplementary"'
cd "$DataPath/data"
****
*
* Download the NIPA data from the BEA web site
*
****
import delimited "https://apps.bea.gov/national/Release/TXT/SeriesRegister.txt", varnames(1) case(preserve) encoding(ISO-8859-1) clear
tempfile SeriesRegister
save `SeriesRegister'

import delimited "https://apps.bea.gov/national/Release/TXT/NipaDataA.txt", varnames(1) case(preserve) encoding(ISO-8859-1) clear
destring Value, replace ignore(`","') float
replace SeriesCode=strupper(SeriesCode)
replace SeriesCode=strtrim(SeriesCode)
merge m:m SeriesCode using `SeriesRegister', generate(_merge1)

keep if _merge1==3

drop _*
rename Period year
sort SeriesCode year
****
*
* Select NIPA series used to for Data-MainDatabase.dta 
*
****

keep if SeriesCode=="A191RC" | SeriesCode=="A191RA" | SeriesCode=="DPCERC" | SeriesCode=="DPCERA" | SeriesCode=="A014RC" ///
      | SeriesCode=="A955RC" | SeriesCode=="A955RA" | SeriesCode=="A782RC" | SeriesCode=="B230RC" | SeriesCode=="DGDSRC"  ///
	  | SeriesCode=="DSERRC" | SeriesCode=="DSERRA" | SeriesCode=="DGDSRA" | SeriesCode=="A253RC" | SeriesCode=="A019RC" ///
	  | SeriesCode=="A646RC" | SeriesCode=="A255RC" | SeriesCode=="B656RC" | SeriesCode=="B253RA" | SeriesCode=="A006RC" ///
	  | SeriesCode=="B646RA" | SeriesCode=="B255RA" | SeriesCode=="B656RA" | SeriesCode=="W172RC" | SeriesCode=="W172RA"

keep year Value SeriesCode
reshape wide Value, i(year) j(SeriesCode, string)
rename ValueA191RA GDP_QI
rename ValueA191RC GDP
rename ValueA014RC INI
rename ValueA955RC GCE
rename ValueA955RA GCE_QI
rename ValueA782RC GIE
rename ValueDPCERC PCE
rename ValueDPCERA PCE_QI
rename ValueDGDSRC PCE_GOOD
rename ValueDGDSRA PCE_GOOD_QI
rename ValueDSERRC PCE_SERV
rename ValueDSERRA PCE_SERV_QI
rename ValueA006RC PIE
rename ValueA019RC NX
rename ValueA253RC EX_GOOD
rename ValueA646RC EX_SERV
rename ValueA255RC IM_GOOD
rename ValueB656RC IM_SERV
rename ValueB253RA EX_GOOD_QI
rename ValueB646RA EX_SERV_QI
rename ValueB255RA IM_GOOD_QI
rename ValueB656RA IM_SERV_QI
rename ValueW172RC X_TOT
rename ValueW172RA X_TOT_QI
rename ValueB230RC POP
keep if year>1946 & year<2018
foreach var in GDP X_TOT PCE GCE PCE_GOOD PCE_SERV EX_GOOD EX_SERV IM_GOOD IM_SERV {
	replace  `var'_QI=`var'_QI/100
	generate `var'_P=`var'/`var'[66]/`var'_QI
	generate `var'_Q=`var'_QI*`var'[66]
}
tsset year, yearly

******
*** Calculating quantities with Fisher-index
******
foreach var in GDP X_TOT PCE_GOOD PCE_SERV GCE EX_GOOD EX_SERV IM_GOOD IM_SERV {
	generate _`var'_LPCQ=L1.`var'_P*   `var'_Q
	generate _`var'_LPLQ=L1.`var'_P*L1.`var'_Q
	generate _`var'_CPCQ=   `var'_P*   `var'_Q
	generate _`var'_CPLQ=   `var'_P*L1.`var'_Q
}

generate _dC_TOT_QI =((_GDP_LPCQ-_X_TOT_LPCQ)/(_GDP_LPLQ-_X_TOT_LPLQ)*(_GDP_CPCQ-_X_TOT_CPCQ)/(_GDP_CPLQ-_X_TOT_CPLQ))^0.5
generate _dC_GOOD_QI=((_PCE_GOOD_LPCQ+_EX_GOOD_LPCQ-_IM_GOOD_LPCQ)/(_PCE_GOOD_LPLQ+_EX_GOOD_LPLQ-_IM_GOOD_LPLQ)*(_PCE_GOOD_CPCQ+_EX_GOOD_CPCQ-_IM_GOOD_CPCQ)/(_PCE_GOOD_CPLQ+_EX_GOOD_CPLQ-_IM_GOOD_CPLQ))^0.5
generate _dC_SERV_QI=((_PCE_SERV_LPCQ+_GCE_LPCQ+_EX_SERV_LPCQ-_IM_SERV_LPCQ)/(_PCE_SERV_LPLQ+_GCE_LPLQ+_EX_SERV_LPLQ-_IM_SERV_LPLQ)*(_PCE_SERV_CPCQ+_GCE_CPCQ+_EX_SERV_CPCQ-_IM_SERV_CPCQ)/(_PCE_SERV_CPLQ+_GCE_CPLQ+_EX_SERV_CPLQ-_IM_SERV_CPLQ))^0.5

******
*** Calculating the Fisher-index ends
******
generate C_TOT_QI=exp(sum(ln(_dC_TOT_QI)))
generate C_TOT_P=(GDP-X_TOT)/(GDP[1]-X_TOT[1])/C_TOT_QI
generate C_TOT=GDP-X_TOT

generate C_GOOD_QI=exp(sum(ln(_dC_GOOD_QI)))
generate C_GOOD_P=(PCE_GOOD+INI+EX_GOOD-IM_GOOD)/(PCE_GOOD[1]+INI[1]+EX_GOOD[1]-IM_GOOD[1])/C_GOOD_QI
generate C_GOOD=PCE_GOOD+INI+EX_GOOD-IM_GOOD

generate C_SERV_QI=exp(sum(ln(_dC_SERV_QI)))
generate C_SERV_P=(PCE_SERV+GCE+EX_SERV-IM_SERV)/(PCE_SERV[1]+GCE[1]+EX_SERV[1]-IM_SERV[1])/C_SERV_QI
generate C_SERV=PCE_SERV+GCE+EX_SERV-IM_SERV


keep year POP GDP GDP_QI GDP_P PCE PIE GCE GIE NX C_TOT C_TOT_QI C_TOT_P C_GOOD C_GOOD_QI C_GOOD_P C_SERV C_SERV_QI C_SERV_P X_TOT X_TOT_QI X_TOT_P

foreach var in GDP_QI GDP_P X_TOT_QI X_TOT_P {
    rename `var' _`var'
	generate `var'=_`var'/_`var'[1]
}

order POP GDP GDP_QI GDP_P PCE PIE GCE GIE NX X_TOT X_TOT_QI X_TOT_P C_TOT C_TOT_QI C_TOT_P C_GOOD C_GOOD_QI C_GOOD_P C_SERV C_SERV_QI C_SERV_P , after(year)


drop _*

label variable POP        "Population, midperiod, thousands"
label variable GDP        "GDP, millions of dollars"
label variable GDP_QI     "GDP, chain-type quantity indexes for GDP, 2009=1"
label variable GDP_P      "GDP, chain-type price indexes for GDP, 2009=1"

label variable PCE        "Personal consumption expenditures, millions of current dollars"
label variable PIE        "Private investment expenditures, millions of current dollars"
label variable GCE        "Government consumption expenditures, millions of current dollars"
label variable GIE        "Government investment expenditures, millions of current dollars"
label variable NX         "Net exports, millions of dollars"

label variable C_TOT      "Consumption expenditures, millions of current dollars"
label variable C_TOT_QI   "Chain-type quantity indexes for consumption, 1947=1"
label variable C_TOT_P    "Chain-type price indexes for consumption, 1947=1"
label variable C_GOOD     "Goods consumption expenditures, millions of current dollars"
label variable C_GOOD_QI  "Chain-type quantity indexes for goods consumption, 1947=1"
label variable C_GOOD_P   "Chain-type price indexes for goods consumption, 1947=1"
label variable C_SERV     "Service consumption expenditures, millions of current dollars"
label variable C_SERV_QI  "Chain-type quantity indexes for service consumption, 1947=1"
label variable C_SERV_P   "Chain-type price indexes for service consumption, 1947=1"

label variable X_TOT      "Fixed investment expenditures, millions of current dollars"
label variable X_TOT_QI   "Chain-type quantity indexes for fixed investment, 1947=1"
label variable X_TOT_P    "Chain-type price indexes for fixed investment, 1947=1"


save "Data-NIPA.dta", replace

exit

