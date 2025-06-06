version 16.0
capture log close
clear all
set more off
global DataPath  `"/Users/Nacho/Dropbox/Omar&Nacho/Baqaee_Burstein/Note on Structural Change/HRV2021_Data_Appendix/Supplementary"'
cd "$DataPath/paper"
graph set window fontface "Times New Roman"
*set scheme lean2
use "Data-MainDatabase.dta", clear

generate VAX_GOOD_S=VAX_GOOD/X_TOT
generate VAX_SERV_S=VAX_SERV/X_TOT

generate VAC_GOOD_S=VAC_GOOD/C_TOT
generate VAC_SERV_S=VAC_SERV/C_TOT

generate C_PC=C_TOT/POP


********
**
**                            Estimate investment technology  
** we use the parmetrization epsilon_x=exp(epsiloncex) and estimate epsilonxex to ensure that epsilon_x>=0
**
********
nl (VAX_GOOD_S={omega_x}*VA_GOOD_P^(1-exp({epsilonxex}))/({omega_x}*VA_GOOD_P^(1-exp({epsilonxex}))+(1-{omega_x})*VA_SERV_P^(1-exp({epsilonxex})) ) ) , ///
   vce(robust) noconstant initial( epsilonxex -1  omega_x 0.5  )   
predict VAX_GOOD_S_hat, yhat
generate VAX_SERV_S_hat=1-VAX_GOOD_S_hat
generate omega_x=_b[/omega_x]
generate epsilon_x=exp(_b[/epsilonxex])

generate VAX_GOOD_QI_hat=(VAX_GOOD_S_hat*C_TOT)/(VAX_GOOD_S_hat[1]*C_TOT[1])/VA_GOOD_P
generate VAX_SERV_QI_hat=(VAX_SERV_S_hat*C_TOT)/(VAX_SERV_S_hat[1]*C_TOT[1])/VA_SERV_P
generate VAX_SERV_GOOD_hat=VAX_SERV_QI_hat/VAX_GOOD_QI_hat

***                                 Goodness of fit 
***  measured by ratio of the percentage increase in the expenditure share in the model relative to the data 1947-2017
egen dVAX_SERV_S_hat_2017=sum((VAX_SERV_S_hat-VAX_SERV_S_hat[1])/(VAX_SERV_S-VAX_SERV_S[1])*(year==2017))


********
**
**                               Estimate log CES preferences
** we use the parmetrization epsilon_c=exp(epsiloncex) and estimate epsiloncex to ensure that epsilon_c>=0
********
nl (VAC_GOOD_S={omega_c}*VA_GOOD_P^(1-exp({epsiloncex}))/({omega_c}*VA_GOOD_P^(1-exp({epsiloncex}))+(1-{omega_c})*VA_SERV_P^(1-exp({epsiloncex})) ) ), ///
   	 vce(robust) noconstant initial( epsiloncex -1  omega_c 0.5 )  
predict VAC_GOOD_S_CES_hat, yhat
generate VAC_SERV_S_CES_hat=1-VAC_GOOD_S_CES_hat
generate omega_c=_b[/omega_c]
generate epsilon_c=exp(_b[/epsiloncex])

generate VAC_GOOD_CES_QI_hat=(VAC_GOOD_S_CES_hat*C_TOT)/(VAC_GOOD_S_CES_hat[1]*C_TOT[1])/VA_GOOD_P
generate VAC_SERV_CES_QI_hat=(VAC_SERV_S_CES_hat*C_TOT)/(VAC_SERV_S_CES_hat[1]*C_TOT[1])/VA_SERV_P
generate VAC_SERV_GOOD_CES_hat=VAC_SERV_CES_QI_hat/VAC_GOOD_CES_QI_hat

***                                 Goodness of fit 
***  measured by ratio of the percentage increase in the expenditure share in the model relative to the data 1947-2017
egen dVAC_SERV_S_CES_hat_2017=sum((VAC_SERV_S_CES_hat-VAC_SERV_S_CES_hat[1])/(VAC_SERV_S-VAC_SERV_S[1])*(year==2017))

********
**
** Estimated PIGL preferences
**
** Goods and service prices relative to investment prices, nominal consumption per capita deflated with the investment price
**
********
nl (VAC_GOOD_S={nu}*(VA_GOOD_P/X_TOT_P/(VA_SERV_P/X_TOT_P))^{gamma}*(VA_SERV_P/X_TOT_P/(C_PC/X_TOT_P))^{chi} ), vce(robust) noconstant initial(nu 0.5  gamma 0.5 chi 0.5 )  
predict VAC_GOOD_S_PIGL_hat, yhat
generate VAC_SERV_S_PIGL_hat=1-VAC_GOOD_S_PIGL_hat
generate chi=_b[/chi]
generate gamma=_b[/gamma]
generate nu=_b[/nu]

generate VAC_GOOD_PIGL_QI_hat=(VAC_GOOD_S_PIGL_hat*C_TOT)/(VAC_GOOD_S_PIGL_hat[1]*C_TOT[1])/VA_GOOD_P
generate VAC_SERV_PIGL_QI_hat=(VAC_SERV_S_PIGL_hat*C_TOT)/(VAC_SERV_S_PIGL_hat[1]*C_TOT[1])/VA_SERV_P
generate VAC_SERV_GOOD_PIGL_hat=VAC_SERV_PIGL_QI_hat/VAC_GOOD_PIGL_QI_hat


***                                 Goodness of fit 
***  measured by ratio of the percentage increase in the expenditure share in the model relative to the data 1947-2017
gen dVAC_SERV_S_PIGL_hat_2017=sum((VAC_SERV_S_PIGL_hat-VAC_SERV_S_PIGL_hat[1])/(VAC_SERV_S-VAC_SERV_S[1])*(year==2017))




label variable VAX_GOOD_S                "Goods value added share in investment expenditure"
label variable VAX_SERV_S                "Services value added share in investment expenditure"
label variable VAC_GOOD_S                "Goods value added share in consumption expenditure"
label variable VAC_SERV_S                "Services value added share in consumption expenditure"

label variable C_PC                      "Nominal consumption expenditure per capita"

label variable VAX_GOOD_S_hat            "Predicted goods value added share in investment expenditure"
label variable VAX_SERV_S_hat            "Predicted services value added share in investment expenditure"

label variable omega_x                   "Weights on goods in the investment technology"
label variable epsilon_x                 "Elasticity of substitution between goods and services in the investment technology"

label variable VAX_GOOD_QI_hat           "Predicted goods quantity in investment expenditure"
label variable VAX_SERV_QI_hat           "Predicted services quantity in investment expenditure"
label variable VAX_SERV_GOOD_hat         "Predicted services over goods quantity in investment expenditure"

label variable VAC_GOOD_S_CES_hat        "Predicted goods value added share in consumption expenditure with log-CES preferences"
label variable VAC_SERV_S_CES_hat        "Predicted services value added share in consumption expenditure with log-CES preferences"

label variable omega_c                   "Weights on goods in the log-CES preferences"
label variable epsilon_c                 "Elasticity of substitution between goods and services in the log-CES preferences"

label variable VAC_GOOD_CES_QI_hat       "Predicted goods quantity in consumption expenditure with log-CES preferences"
label variable VAC_SERV_CES_QI_hat       "Predicted services quantity iin consumption expenditure with log-CES preferences"
label variable VAC_SERV_GOOD_CES_hat     "Predicted services over goods quantity in consumption expenditure with log-CES preferences"

label variable VAC_GOOD_S_PIGL_hat       "Predicted goods value added share in consumption expenditure with log-PIGL preferences"
label variable VAC_SERV_S_PIGL_hat       "Predicted services value added share in consumption expenditure with log-PIGL preferences"

label variable chi                       "Elasticity of goods expenditure share with respect to total expenditures in PIGL preferences"
label variable gamma                     "Elasticity of goods expenditure share with respect to relative prices in PIGL preferences"
label variable nu                        "Scale parameter for goods expenditure share in PIGL preferences"

label variable VAC_GOOD_PIGL_QI_hat      "Predicted goods quantity in consumption expenditure with PIGL preferences"
label variable VAC_SERV_PIGL_QI_hat      "Predicted services quantity iin consumption expenditure with PIGL preferences"
label variable VAC_SERV_GOOD_PIGL_hat    "Predicted services over goods quantity in consumption expenditure with PIGL preferences"


label variable dVAX_SERV_S_hat_2017      "Goodness of fitfor services value added share in investment expenditure"
label variable dVAC_SERV_S_CES_hat_2017  "Goodness of fit for services value added share in consumption expenditure with log-CES"
label variable dVAC_SERV_S_PIGL_hat_2017 "Goodness of fit  for services value added share in consumption expenditure with PIGL"


save "_Data-MainDatabase1.dta", replace
