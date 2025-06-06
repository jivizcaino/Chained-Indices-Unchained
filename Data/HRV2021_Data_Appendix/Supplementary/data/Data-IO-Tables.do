version 16.0
capture log close
clear all
set more off
global DataPath `"/Users/akos/Dropbox/AkosBertRichard/Investment/Submissions/FinalSubmission/Supplementary"'
cd "$DataPath/rawdata"
***
***     Constructing equipment and IPP investment for the government for periods 
***   1947-1962 and 1963-1996 by assigning government investment expenditures on
***      commodity "Professional, scientific, and technical services"
***
*** Use tables 1947-1962
***
use "IOUseBeforeRedefinitions1947-1962.dta", clear
generate _FinalUse8=FinalUse8+FinalUse11
generate _FinalUse9=FinalUse9+FinalUse12
generate _FinalUse10=FinalUse10+FinalUse13
drop FinalUse8 FinalUse9 FinalUse10 FinalUse11 FinalUse12 FinalUse13
rename _FinalUse8 FinalUse8
rename _FinalUse9 FinalUse9
rename _FinalUse10 FinalUse10
generate FinalUse11=0 if FinalUse10!=.
replace  FinalUse11=FinalUse10 if naics=="54"
replace  FinalUse10=0 if naics=="54"
label variable FinalUse8  "General government consumption expenditures"
label variable FinalUse9  "General government gross investment in structures"
label variable FinalUse10 "General government gross investment in equipment"
label variable FinalUse11 "General government gross investment in intellectual property products"
order FinalUse8 FinalUse9 FinalUse10 FinalUse11, after(FinalUse7)

sort year Industry
by year: egen _FinalUse10=sum(FinalUse10*(Industry<49))
by year: egen _FinalUse11=sum(FinalUse11*(Industry<49))
replace FinalUse10=_FinalUse10 if Industry==51 
replace FinalUse11=_FinalUse11 if Industry==51 
drop _*
save "../data/Data-IO-TablesUse1947-1962.dta", replace
***
*** Use tables 1963-1996
***
use "IOUseBeforeRedefinitions1963-1996.dta", clear
generate _FinalUse8=FinalUse8+FinalUse11+FinalUse14
generate _FinalUse9=FinalUse9+FinalUse12+FinalUse15
generate _FinalUse10=FinalUse10+FinalUse13+FinalUse16
drop FinalUse8-FinalUse16
rename _FinalUse8 FinalUse8
rename _FinalUse9 FinalUse9
generate FinalUse10=_FinalUse10 
replace  FinalUse10=0 if naics=="5415" | naics=="5412OP"
generate FinalUse11=0 
replace FinalUse11=_FinalUse10 if naics=="5415" | naics=="5412OP"
drop _*
sort year Industry
by year: egen _FinalUse10=sum(FinalUse10*(Industry<68))
by year: egen _FinalUse11=sum(FinalUse11*(Industry<68))
replace FinalUse10=_FinalUse10 if Industry==70 
replace FinalUse11=_FinalUse11 if Industry==70 

label variable FinalUse8  "General government consumption expenditures"
label variable FinalUse9  "General government gross investment in structures"
label variable FinalUse10 "General government gross investment in equipment"
label variable FinalUse11 "General government gross investment in intellectual property products"
order FinalUse8 FinalUse9 FinalUse10 FinalUse11, after(FinalUse7)
drop _*

save "../data/Data-IO-TablesUse1963-1996.dta", replace
***
*** Use tables 1997-2017
***
use "IOUseBeforeRedefinitions1997-2017.dta", clear
replace FinalUse9=FinalUse9+FinalUse13+FinalUse17
replace FinalUse10=FinalUse10+FinalUse14+FinalUse18
replace FinalUse11= FinalUse11+ FinalUse15+ FinalUse19
replace FinalUse12= FinalUse12+ FinalUse16+ FinalUse20
label variable FinalUse9  "General government consumption expenditures"
label variable FinalUse10  "General government gross investment in structures"
label variable FinalUse11 "General government gross investment in equipment"
label variable FinalUse12 "General government gross investment in intellectual property products"
drop FinalUse13 FinalUse14 FinalUse15 FinalUse16 FinalUse17 FinalUse18 FinalUse19 FinalUse20
save "../data/Data-IO-TablesUse1997-2017.dta", replace
***
*** Make tables 1947-2017
***
use  "IOMakeBeforeRedefinitions1947-1962.dta", clear
save "../data/Data-IO-TablesMake1947-1962.dta", replace
use  "IOMakeBeforeRedefinitions1963-1996.dta", clear
save "../data/Data-IO-TablesMake1963-1996.dta", replace
use  "IOMakeBeforeRedefinitions1997-2017.dta", clear
save "../data/Data-IO-TablesMake1997-2017.dta", replace

exit

