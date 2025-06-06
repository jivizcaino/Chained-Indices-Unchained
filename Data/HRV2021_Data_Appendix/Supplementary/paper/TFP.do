version 16.0
capture log close
clear all
set more off
global DataPath  `"/Users/Nacho/Dropbox/Omar&Nacho/Baqaee_Burstein/Note on Structural Change/HRV2021_Data_Appendix/Supplementary"'
cd "$DataPath/paper"
graph set window fontface "Times New Roman"
*set scheme lean2
use "_Data-MainDatabase1.dta", clear

generate VAX_S_GOOD=0.5*(VAX_GOOD/X_TOT+L1.VAX_GOOD/L1.X_TOT)
generate VAX_S_SERV=0.5*(VAX_SERV/X_TOT+L1.VAX_SERV/L1.X_TOT)

generate VAC_S_GOOD=0.5*(VAC_GOOD/C_TOT+L1.VAC_GOOD/L1.C_TOT)
generate VAC_S_SERV=0.5*(VAC_SERV/C_TOT+L1.VAC_SERV/L1.C_TOT)

generate VAX_GOOD_QI=VAX_GOOD/VAX_GOOD[1]/VA_GOOD_P
generate VAX_SERV_QI=VAX_SERV/VAX_SERV[1]/VA_SERV_P

generate VAC_GOOD_QI=VAC_GOOD/VAC_GOOD[1]/VA_GOOD_P
generate VAC_SERV_QI=VAC_SERV/VAC_SERV[1]/VA_SERV_P

generate VA_SERV_GOOD_P=VA_SERV_P/VA_GOOD_P
generate VA_SERV_X_P=VA_SERV_P/X_TOT_P
generate VA_GOOD_X_P=VA_GOOD_P/X_TOT_P
generate C_X_P=C_TOT_P/X_TOT_P

egen mLAB_S_TOT=mean(LAB_S_TOT)
egen mLAB_S_GOOD=mean(LAB_S_GOOD)
egen mLAB_S_SERV=mean(LAB_S_SERV)


*******
**
**
** Bottom up approach to calculate the effective investment specific TFP 
**
**
*******

*******
**
** Calculating TFP for the aggregate economy and goods and services sector
**  with growth accounting using constant aggregate factor shares 
**
******

foreach sec in TOT GOOD SERV {
	generate dlTFPVA_`sec'_I=ln(VA_`sec'_QI)-ln(L1.VA_`sec'_QI)-(1-mLAB_S_TOT)*(ln(CAP_`sec'_QI)-ln(L1.CAP_`sec'_QI))-mLAB_S_TOT*(ln(LAB_`sec'_QI)-ln(L1.LAB_`sec'_QI))
}

generate TFPVA_TOT_I =exp(sum(dlTFPVA_TOT_I ))
generate TFPVA_GOOD_I=exp(sum(dlTFPVA_GOOD_I))
generate TFPVA_SERV_I=exp(sum(dlTFPVA_SERV_I))

*******
**
** Calculating exogenous investment specific technical change with growth accounting
**
******

generate dlA_X_I=ln(X_TOT_QI)-ln(L1.X_TOT_QI)-VAX_S_GOOD*(ln(VAX_GOOD_QI)-ln(L1.VAX_GOOD_QI))-VAX_S_SERV*(ln(VAX_SERV_QI)-ln(L1.VAX_SERV_QI))

generate A_X_I=exp(sum(dlA_X_I))

*******
**
** Calculating the effective investment specific TFP 
**
******

generate calA_X_I_BU=A_X_I*1/(omega_x/TFPVA_GOOD_I+(1-omega_x)/TFPVA_SERV_I)

*******
**
** Effective consumption specific TFP if preferences are log-CES
**
******

generate calA_C_I=1/(omega_c/TFPVA_GOOD_I+(1-omega_c)/TFPVA_SERV_I)

*******
**
** Top down approach to calculate the effective investment specific TFP 
**
*******

*******
**
** Calculating effective investment specific TFP with growth accounting using constant aggregate factor shares and GDP is unit of the investment good
**
******

generate dlcalA_X_I_TD=(ln(VA_TOT)-ln(L1.VA_TOT)-(ln(X_TOT_P)-ln(L1.X_TOT_P)))-(1-mLAB_S_TOT)*(ln(CAP_TOT_QI)-ln(L1.CAP_TOT_QI))-mLAB_S_TOT*(ln(LAB_TOT_QI)-ln(L1.LAB_TOT_QI))

generate calA_X_I_TD=exp(sum(dlcalA_X_I_TD))


label variable VAX_S_GOOD                "Goods value added share in investment expenditure avergared over t and t-1"
label variable VAX_S_SERV                "Services value added share in investment expenditure avergared over t and t-1"
label variable VAC_S_GOOD                "Goods value added share in consumption expenditure avergared over t and t-1"
label variable VAC_S_SERV                "Services value added share in consumption expenditure avergared over t and t-1"


label variable VAX_GOOD_QI               "Goods quantity in investment expenditure, 1947 = 1"
label variable VAX_SERV_QI               "Services quantity in investment expenditure, 1947 = 1"
label variable VAC_GOOD_QI               "Goods quantity in consumption expenditure, 1947 = 1"
label variable VAC_SERV_QI               "Services quantity iin consumption expenditure, 1947 = 1"

label variable VA_SERV_GOOD_P            "Value added, services relative to goods price index, 1947 = 1"
label variable VA_GOOD_X_P               "Value added, goods relative to investment price index, 1947 = 1"
label variable VA_SERV_X_P               "Value added, services relative to investment price index, 1947 = 1"
label variable C_X_P                     "Consumption relative to investment price index, 1947 = 1"

label variable TFPVA_TOT_I               "TFP index, aggregate, 1947 = 1"
label variable TFPVA_GOOD_I              "TFP index, goods, 1947 = 1"
label variable TFPVA_SERV_I              "TFP index, services, 1947 = 1"

label variable A_X_I                     "Exogenous investment specific TFP index, 1947 = 1"
label variable calA_X_I_BU               "Effective investment specific TFP index, buttom up approach, 1947 = 1"
label variable calA_X_I_TD               "Effective investment specific TFP index, top down approach, 1947 = 1"

label variable calA_C_I                  "Effective consumption specific TFP index, log-CES, 1947 = 1"


label variable mLAB_S_TOT                "Average labour share, total economy"
label variable mLAB_S_GOOD               "Average labour share, goods sector"
label variable mLAB_S_SERV               "Average labour share, services sector"

drop  dlTFPVA_TOT_I dlTFPVA_GOOD_I dlTFPVA_SERV_I dlA_X_I dlcalA_X_I_TD

save "_Data-MainDatabase2.dta", replace


exit






