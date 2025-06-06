version 16.0
capture log close
clear all
set more off
global DataPath `"/Users/akos/Dropbox/AkosBertRichard/Investment/Submissions/FinalSubmission/Supplementary"'
cd "$DataPath/paper"
set more off
graph set window fontface "Times New Roman"
set scheme lean2
use "_Data-MainDatabase3.dta", clear

*******
**
** Assumption 2
** 
*******
local delta 0.08
local rho   0.04
generate dlcalA_X_I_BU=ln(calA_X_I_BU)-ln(L1.calA_X_I_BU)
generate dlTFPVA_SERV_I=ln(TFPVA_SERV_I)-ln(L1.TFPVA_SERV_I)

tssmooth ma madlcalA_X_I_BU = dlcalA_X_I_BU, window(4 1)
tssmooth ma madlTFPVA_SERV_I = dlTFPVA_SERV_I, window(4 1)


generate A2_1=chi*((1-mLAB_S_TOT)/mLAB_S_TOT*madlcalA_X_I_BU+madlTFPVA_SERV_I)
generate A2_2=`delta'+`rho'+chi*(madlcalA_X_I-madlTFPVA_SERV_I)+(1-chi)*madlcalA_X_I/mLAB_S_TOT


**
** Condition 1 in Assumptiom 2 holds if the maximum of A2_1<rho=0.04 
** Condition 2 in Assumptiom 2 holds if the minimum of A2_2>0 
**

sum A2_1 A2_2

*******
**
** Assumption 6
** 
*******

generate A6=(C_PC/X_TOT_P)^chi-(1-chi)/(1-gamma)*nu*(VA_SERV_P/X_TOT_P)^gamma*(VA_SERV_P/X_TOT_P)^(chi-gamma)

**
** Assumptiom 6 holds if the minimum of A6>0
**

sum A6

exit
