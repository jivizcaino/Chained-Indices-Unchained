version 16.0
capture log close
clear all
set more off
global DataPath `"/Users/akos/Dropbox/AkosBertRichard/Investment/Submissions/FinalSubmission/Supplementary"'
cd "$DataPath/paper"
set more off
graph set window fontface "Times New Roman"
set scheme lean2
use "_Data-MainDatabase2.dta", clear
*******
**
** Graphs
**
*******

******
**
** Figure 1: Share of consumption in GDP
**
*******
generate C1_S=(PCE+GCE+NX)/(PCE+PIE+GCE+GIE+NX)
generate C2_S=(PCE+GCE)/(PCE+PIE+GCE+GIE)
generate C3_S=PCE/(PCE+PIE)
graph twoway ///
	(tsline C1_S C2_S C3_S, lwidth(medthick medthick medthick) lpattern(solid dash shortdash)) ,  ///
	graphregion(margin(tiny)) plotregion(margin(zero))                                ///
	ylabel(0.4(0.1)0.9, format(%9.2fc) labsize(large))  yscale(range(0.4 .9001))                         ///
	ytitle("")                                                                        ///
	xlabel(1950(10)2010, labsize(large)) tmtick(1947(1)2017)                                          ///
	xtitle("")                                                                        ///
	legend(order(1 "Net export included in consumption" 2 "Net export distributed pro rata"  3 "Net export and government consumption expenditures" "distributed pro rata") ///
	       rows(3) position(12) ring(0) bplacement(b) region(fcolor(none)) bmargin(zero) size(large))     ///
	title("Share of consumption expenditures in GDP", size(vlarge))                                       ///
	name(Graph01, replace)
******
**
** Figure 1: Final investment expenditure
**
*******
local var1 VAX_GOOD_S
local var2 VAX_SERV_S
local xcor1=year[12]
local ycor1=`var1'[12]
local xcor2=year[12]
local ycor2=`var2'[12]
graph twoway ///
	(tsline VAX_GOOD_S VAX_SERV_S, lwidth(medthick medthick) lpattern(solid dash) ) , ///
	graphregion(margin(tiny)) plotregion(margin(zero))                                ///
	ylabel(0(0.1)0.9, format(%9.1fc) labsize(large))  yscale(range(0 .9001))                         ///
	ytitle("")                                                                        ///
	xlabel(1950(10)2010, labsize(large)) tmtick(1947(1)2017)                                          ///
	xtitle("")                                                                        ///
	legend(off )                                                                      ///
	title("Final investment expenditure", size(vlarge))                                        ///
	text(`ycor1' `xcor1'  "Goods value added share",    place(ne) color(black) justification(left) margin(tiny) alignment(top) size(large))    ///
	text(`ycor2' `xcor2'  "Services value added share", place(se) color(black) justification(left) margin(tiny) size(large))  ///
	name(Graph02, replace)
******
**
** Figure 1: Final consumption expenditure
**
*******
local var1 VAC_GOOD_S
local var2 VAC_SERV_S
local xcor1=year[30]
local ycor1=`var1'[30]
local xcor2=year[30]
local ycor2=`var2'[30]
graph twoway ///
	(tsline `var1' `var2' , lwidth(medthick medthick) lpattern(solid dash)) ,  ///
	graphregion(margin(tiny)) plotregion(margin(zero))                                ///
	ylabel(0(0.1)0.9, format(%9.1fc) labsize(large))  yscale(range(0 .9001))                         ///
	ytitle("")                                                                        ///
	xlabel(1950(10)2010,labsize(large)) tmtick(1947(1)2017)                                          ///
	xtitle("")                                                                        ///
	legend(off)                                                                       ///
	title("Final consumption expenditure", size(vlarge))                                       ///
	text(`ycor1' `xcor1'  "Goods value added share",    place(sw) color(black) justification(right) margin(tiny) alignment(top) size(large))	///
	text(`ycor2' `xcor2'  "Services value added share", place(nw) color(black) justification(right) margin(tiny) alignment(top) size(large)) ///
	name(Graph03, replace)
*******
**
** Figure 2: Price of services relative to goods, price of services relative to investment, price of goods relative to investment
**
*******
local var1 VA_SERV_GOOD_P
local var2 VA_SERV_X_P
local var3 VA_GOOD_X_P
local xcor1=year[52]
local ycor1=`var1'[52]
local xcor2=year[44]
local ycor2=`var2'[44]
local xcor3=year[38]
local ycor3=`var3'[38]
graph twoway ///
    (tsline  `var1' `var2' `var3' , lwidth(medthick medthick  medthick)),                                 ///
	graphregion(margin(tiny)) plotregion(margin(zero))                                ///
	ylabel(0.80 1.00 1.25 1.56 1.95, format(%9.2fc) labsize(large))  yscale(range(0.8 1.95001) log)    ///
	ytitle("")                                                                        ///
	xlabel(1950(10)2010, labsize(large)) tmtick(1947(1)2017)                                          ///
	xtitle("")                                                                        ///
	legend(off )                                                                      ///
	title(" ", size(vlarge))                                                                        ///
	subtitle("1947=1, log scale", size(large))                                                   ///
	text(`ycor1' `xcor1'   "Price of services relative to goods",        place(nw) color(black) justification(right) size(large))       ///
	text(`ycor2' `xcor2'   "Price of services relative" "to investment", place(se) color(black) justification(center) margin(tiny) size(large)) ///
	text(`ycor3' `xcor3'   "Price of goods relative to investment",      place(n)  color(black) justification(center) margin(tiny) size(large)) ///
	name(Graph04, replace)
*******
**
** Figure 2: Price of consumption relatvie to investment
**
******* 
local var1 C_X_P
local xcor1=year[50]
local ycor1=`var1'[50]
graph twoway ///
    (tsline  `var1', lwidth(medthick)),                                 ///
	graphregion(margin(tiny)) plotregion(margin(zero))                                ///
	ylabel(0.80 1.00 1.25 1.56 1.95, format(%9.2fc) labsize(large))  yscale(range(0.8 1.95001) log)    ///
	ytitle("")                                                                        ///
	xlabel(1950(10)2010, labsize(large)) tmtick(1947(1)2017)                                          ///
	xtitle("")                                                                        ///
	title(" ", size(vlarge))                                                                        ///
	legend(off )                                                                      ///
	subtitle("1947=1, log scale", size(large))                                                   ///
	text(`ycor1' `xcor1'  "Price of consumption relative to investment", place(nw) color(black) justification(right)  margin(tiny) size(large)) ///
	name(Graph05, replace)
		
*******
**
** Figure 2: Real services value added over real goods value added
**
*******
generate VAC_SERV_GOOD=VAC_SERV_QI/VAC_GOOD_QI
generate VAX_SERV_GOOD=VAX_SERV_QI/VAX_GOOD_QI
local var1 VAC_SERV_GOOD
local var2 VAX_SERV_GOOD
local xcor1=year[47]
local ycor1=`var1'[47]
local xcor2=year[64]
local ycor2=`var2'[64]
graph twoway ///
	(tsline `var1' `var2' , lwidth(medthick medthick) lpattern(solid dash)), ///
	graphregion(margin(tiny)) plotregion(margin(zero))                                ///
	ylabel(0.75(0.25)2.25, format(%9.2fc) labsize(large))  yscale(range(0.75 2.25001))                    ///
	ytitle("")                                                                        ///
	xlabel(1950(10)2010, labsize(large)) tmtick(1947(1)2017)                                          ///
	xtitle("")                                                                        ///
	legend(off)                                                                       ///
	subtitle("1947=1", size(large))                                                           ///
	title("Real services value added over real goods value added", size(vlarge))               ///
	text(`ycor1' `xcor1'  "Final consumption expenditure", place(nw) color(black) justification(right) margin(tiny) size(large)) ///
	text(`ycor2' `xcor2'  "Final investment expenditure",  place(nw) color(black) justification(right) margin(tiny) size(large))  ///
	name(Graph06, replace)
*******
**
** Figure 3: Shares in investment: Nominal expenditure shares
**
*******	
local var1 VAX_GOOD_S
local var2 VAX_SERV_S
local xcor1=year[12]
local ycor1=`var1'[12]
local xcor2=year[12]
local ycor2=`var2'[12]
graph twoway ///
	(tsline `var1' VAX_GOOD_S_hat, lwidth(medthick medthick) lpattern(solid dash)  lcolor(black black))    ///
	(tsline `var2' VAX_SERV_S_hat, lwidth(medthick medthick) lpattern(solid dash)  lcolor(black black)),   ///
	graphregion(margin(tiny)) plotregion(margin(zero))                                ///
	ylabel(0(0.1)0.9, format(%9.1fc) labsize(large))  yscale(range(0 .9001))                         ///
	ytitle("")                                                                        ///
	xlabel(1950(10)2010, labsize(large)) tmtick(1947(1)2017)                                          ///
	xtitle("")                                                                        ///
	title("Nominal expenditure shares", size(vlarge))               ///
	legend(order(1 "Data" 2 "Model")   rows(1) position(6) ring(0) region(fcolor(none)) size(large))                  ///
	text(`ycor1' `xcor1' "Goods value added share",     place(ne) color(black) justification(left) margin(tiny) fcolor(none) size(large))     ///
	text(`ycor2' `xcor2' "Services value added share",  place(se) color(black) justification(left) margin(tiny) fcolor(none) size(large))  ///
	name(Graph07, replace)

*******
**
** Figure 3: Shares in investment: Real services value added over real goods value added
**
*******
local var1 VAX_SERV_GOOD
local var2 VAX_SERV_GOOD_hat
graph twoway ///
	(tsline `var1' `var2', lwidth(medthick medthick) lpattern(solid dash) lcolor(black black)), ///
	graphregion(margin(tiny)) plotregion(margin(zero))                                ///
	ylabel(0.6(0.1)1.5, format(%9.1fc) labsize(large))  yscale(range(0.6 1.50001))                    ///
	ytitle("")                                                                        ///
	xlabel(1950(10)2010, labsize(large)) tmtick(1947(1)2017)                                          ///
	xtitle("")                                                                        ///
	legend(order(1 "Data" 2 "Model")   rows(1) position(6) ring(0) region(fcolor(none)) size(large) )                  ///
	title("Real services value added over real goods value added, 1947=1", size(vlarge))               ///
	name(Graph07a, replace)
	
*******
**
** Figure 4: Shares of consumption: log-CES preferences, Nominal shares
**
*******	
local var1 VAC_GOOD_S
local var2 VAC_SERV_S
local xcor1=year[48]
local ycor1=`var1'[48]
local xcor2=year[48]
local ycor2=`var2'[48]
graph twoway ///
	(tsline  `var1' VAC_GOOD_S_CES_hat, lwidth(medthick medthick) lpattern(solid dash)  lcolor(black black))    ///
	(tsline  `var2' VAC_SERV_S_CES_hat, lwidth(medthick medthick) lpattern(solid dash)  lcolor(black black)),   ///
	graphregion(margin(tiny)) plotregion(margin(zero))                                ///
	ylabel(0(0.1)0.9, format(%9.1fc) labsize(large))  yscale(range(0 .9001))                         ///
	ytitle("")                                                                        ///
	xlabel(1950(10)2010, labsize(large)) tmtick(1947(1)2017)                                          ///
	xtitle("")                                                                        ///
	title("log-CES preferences" "Nominal shares", size(vlarge))                            ///
	legend(order(1 "Data" 2 "Model")   rows(1) position(9) ring(0) bplacement(c) region(fcolor(none)) size(large) bmargin(zero))     ///
	text(`ycor1' `xcor1'  "Goods value added share",    place(sw) color(black)  justification(right)  margin(tiny) fcolor(none) alignment(top) size(large))	          ///
	text(`ycor2' `xcor2'  "Services value added share", place(nw) color(black)  justification(right)  margin(tiny) fcolor(none) alignment(top) size(large))           ///
	name(Graph08, replace)
*******
**
** Figure 4: Shares of consumption: log-CES preferences, real services value added over real goods value added
**
*******
local var1 VAC_SERV_GOOD
local var2 VAC_SERV_GOOD_CES_hat
graph twoway ///
	(tsline `var1' `var2', lwidth(medthick medthick medthick medthick) lpattern(solid dash) lcolor(black black)), ///
	graphregion(margin(tiny)) plotregion(margin(zero))                                ///
	ylabel(0.6(0.2)2.4, format(%9.1fc)  labsize(large))  yscale(range(0.6 2.4001))                    ///
	ytitle("")                                                                        ///
	xlabel(1950(10)2010, labsize(large)) tmtick(1947(1)2017)                                          ///
	xtitle("")                                                                        ///
	legend(order(1 "Data" 2 "Model")   rows(1) position(6) ring(0) region(fcolor(none)) size(large) bmargin(zero))                  ///
	title("log-CES preferences" "Real services value added over real goods value added, 1947=1", size(vlarge))               ///
	name(Graph08a, replace)
*******
**
** Figure 4: Shares of consumption: PIGL preferences, nominal shares
**
*******	
local var1 VAC_GOOD_S
local var2 VAC_SERV_S
local xcor1=year[45]
local ycor1=`var1'[45]
local xcor2=year[45]
local ycor2=`var2'[45]
graph twoway ///
	(tsline VAC_GOOD_S VAC_GOOD_S_PIGL_hat, lwidth(medthick medthick) lpattern(solid dash)  lcolor(black black))    ///
	(tsline VAC_SERV_S VAC_SERV_S_PIGL_hat, lwidth(medthick medthick) lpattern(solid dash)  lcolor(black black)),   ///
	graphregion(margin(tiny)) plotregion(margin(zero))                               ///
	ylabel(0(0.1)0.9, format(%9.1fc) labsize(large))  yscale(range(0 .9001))                         ///
	ytitle("")                                                                        ///
	xlabel(1950(10)2010, labsize(large)) tmtick(1947(1)2017)                                          ///
	xtitle("")                                                                        ///
	title("PIGL preferences" "Nominal shares", size(vlarge) )                                        ///
	legend(order(1 "Data" 2 "Model")   rows(1) position(9) ring(0) bplacement(c) size(large) bmargin(zero))     ///
	text(`ycor1' `xcor1'  "Goods value added share",    place(sw) color(black) justification(right) margin(tiny) fcolor(none) alignment(top) size(large)) ///
	text(`ycor2' `xcor2'  "Services value added share", place(nw) color(black) justification(right) margin(tiny) fcolor(none) alignment(top) size(large)) ///
	name(Graph09, replace)
*******
**
** Figure 4: Shares of consumption: PIGL preferences, real services value added over real goods value added
**
*******
local var1 VAC_SERV_GOOD
local var2 VAC_SERV_GOOD_PIGL_hat
graph twoway ///
	(tsline `var1' `var1a' `var2' `var2a', lwidth(medthick medthick medthick medthick) lpattern(solid dash solid dash) lcolor(black black black black)), ///
	graphregion(margin(tiny)) plotregion(margin(zero))                                ///
	ylabel(0.6(0.2)2.4, format(%9.1fc) labsize(large))  yscale(range(0.6 2.4001))                    ///
	ytitle("")                                                                        ///
	xlabel(1950(10)2010, labsize(large)) tmtick(1947(1)2017)                                          ///
	xtitle("")                                                                        ///
	legend(order(1 "Data" 2 "Model")   rows(1) position(6) ring(0) region(fcolor(none)) size(large) bmargin(zero))                  ///
	title("PIGL preferences" "Real services value added over real goods value added, 1947=1" , size(vlarge))               ///
	name(Graph09a, replace)

*******
**
** Figure 5: Contribution of changes in income and relative prices under PIGL utility
**
*******	
generate VAC_GOOD_S_INC=nu*(VA_GOOD_P[32]/VA_SERV_P[32])^gamma*(VA_SERV_P/C_PC)^chi
generate VAC_SERV_S_INC=1-VAC_GOOD_S_INC
generate VAC_GOOD_S_PRICE=nu*(VA_GOOD_P/VA_SERV_P)^gamma*(VA_SERV_P[32]/C_PC[32])^chi
generate VAC_SERV_S_PRICE=1-VAC_GOOD_S_PRICE
local var1 VAC_GOOD_S
local var2 VAC_SERV_S
local xcor1=year[50]
local ycor1=`var1'[50]
local xcor2=year[50]
local ycor2=`var2'[50]
graph twoway ///
	(tsline VAC_GOOD_S VAC_GOOD_S_INC VAC_GOOD_S_PRICE, lwidth(medthick medthick medthick) lpattern(solid dash shortdash)  lcolor(black black black))   ///
	(tsline VAC_SERV_S VAC_SERV_S_INC VAC_SERV_S_PRICE, lwidth(medthick medthick medthick) lpattern(solid dash shortdash)  lcolor(black black black)),  ///
	graphregion(margin(small)) plotregion(margin(zero))                                ///
	ylabel(0(0.1)0.9, format(%9.1fc) labsize(large))  yscale(range(0 0.9001))                         ///
	ytitle("")                                                                         ///
	xlabel(1950(10)2010, labsize(large)) tmtick(1947(1)2017)                                           ///
	xtitle("")                                                                         ///
	legend(order(1 "Data" 2 "Income Effect" 3 "Price Effect")   rows(1) position(9) ring(0) region(fcolor(none)) bplacement(c) size(large))                        ///
	text(`ycor1' `xcor1'  "Goods value added share",    place(sw) color(black) justification(right) margin(tiny) alignment(top) size(large))  ///
	text(`ycor2' `xcor2'  "Services value added share", place(nw) color(black) justification(right) margin(tiny) alignment(top) size(large))  ///
	name(Graph10, replace)	
*******
**
** Figure 6: Sectoral TFPs: Exogenous sector TFPs
**
*******	
local var1 TFPVA_GOOD_I
local var2 TFPVA_SERV_I
local var3 A_X_I
local xcor1=year[70]
local ycor1=`var1'[70]
local xcor2=year[70]
local ycor2=`var2'[70]
local xcor3=year[30]
local ycor3=`var3'[30]
graph twoway ///
    (tsline `var1' `var2' `var3', lwidth(medthick medthick medthick) lpattern(solid dash shortdash) lcolor(black black black)), ///
	graphregion(margin(tiny)) plotregion(margin(zero))                                                    ///
	ylabel(0.85 1 1.176 1.384	1.628 1.916 2.254, format(%9.1fc) labsize(large))  yscale(range(0.85 2.254001) log)      ///
	ytitle("")                                                                        ///
	xlabel(1950(10)2010, labsize(large)) tmtick(1947(1)2017)                                          ///
	xtitle("")                                                                        ///
	legend(off )                                                                      ///
	subtitle("1947=1, log scale", size(large))                                                     ///
	text(`ycor1' `xcor1'  "Goods TFP", place(nw) color(black)  justification(right) margin(small) size(large))           ///
	text(`ycor2' `xcor2'  "Services TFP ", place(nw) color(black) justification(right) margin(medsmall) size(large))     ///
	text(`ycor3' `xcor3'  "Exogenous investment-specific TFP ", place(se) color(black) justification(left) size(large))  ///
	name(Graph11, replace)	
*******
**
** Figure 6: Sectoral TFPs: Effective investment-specific TFP
**
*******	

generate lcalA_X_I_BU=ln(calA_X_I_BU)
reg lcalA_X_I_BU year
predict lcalA_X_I_BU_hat
generate calA_X_I_BU_hat=exp(lcalA_X_I_BU_hat)

local var1 calA_X_I_BU
local xcor1=year[34]
local ycor1=`var1'[34]
graph twoway ///
    (tsline calA_X_I_BU calA_X_I_BU_hat, lwidth(medthick medthick medthick) lpattern(solid dash shortdash) lcolor(black black black)), ///
	graphregion(margin(tiny)) plotregion(margin(zero))                                                  ///
	ylabel(0.85 1 1.176 1.384	1.628 1.916 2.254, format(%9.1fc) labsize(large))  yscale(range(0.85 2.254001) log)    ///
	ytitle("")                                                                        ///
	xlabel(1950(10)2010, labsize(large)) tmtick(1947(1)2017)                                          ///
	xtitle("")                                                                        ///
	legend(off )                                                                      ///
	subtitle("1947=1, log scale", size(large) )                                                     ///
	text(`ycor1' `xcor1'   "Effective Investment-specific TFP", place(se) color(black)  justification(center) margin(small) size(large))          ///
	name(Graph12, replace)	

*******
**
** Figure 6: Sectoral TFPs:  Effective investment-specific TFP, top down and bottom up measure
**
*******	
local var1 calA_X_I_BU
local var2 calA_X_I_TD
local xcor1=year[35]
local ycor1=`var1'[35]
local xcor2=year[41]
local ycor2=`var2'[41]
graph twoway ///
    (tsline `var1' `var2', lwidth(medthick medthick medthick) lpattern(solid dash shortdash) lcolor(black black black)), ///
	graphregion(margin(tiny)) plotregion(margin(zero))                                        ///
	ylabel(0.85 1 1.176 1.384	1.628 1.916 2.254, format(%9.1fc) labsize(large))  yscale(range(0.85 2.254001) log)       ///
	ytitle("")                                                                        ///
	xlabel(1950(10)2010, labsize(large)) tmtick(1947(1)2017)                                          ///
	xtitle("")                                                                        ///
	legend(off )                                                                      ///
	subtitle("1947=1, log scale", size(large))                                                   ///
	text(`ycor1' `xcor1'   "Bottom up (paper)", place(s) color(black)  justification(right) margin(medsmall) size(large))          ///
	text(`ycor2' `xcor2'   "Top down", place(nw) color(black)  justification(right) margin(small) size(large))          ///
	name(Graph12a, replace)	
*******
**
** Figure 7: GDP per capita growth: Numeraire investment (model) versus Fisher index (NIPA)
**
*******	
generate GDP_QI_X=GDP/GDP[1]/X_TOT_P
local var1 GDP_QI_X
local xcor1=year[70]
local ycor1=`var1'[70]
local var2 GDP_QI
local xcor2=year[55]
local ycor2=`var2'[55]
graph twoway ///
    (tsline GDP_QI GDP_QI_X, lwidth(medthick medthick medthick) lpattern(solid dash shortdash)), ///
	graphregion(margin(tiny)) plotregion(margin(zero))                                        ///
	ylabel(1.00	1.70	2.87	4.87	8.26	14, format(%9.1fc) labsize(large))  yscale(range(0.99 14.0001) log)       ///
	ytitle("")                                                                        ///
	xlabel(1950(10)2010, labsize(large)) tmtick(1947(1)2017)                                          ///
	xtitle("")                                                                        ///
	legend(off )                                                                      ///
	subtitle("1947=1, log scale", size(large))                                                   ///
	text(`ycor1' `xcor1'   "Numeraire investment", place(nw) color(black)  justification(right) margin(small) size(large))          ///
	text(`ycor2' `xcor2'   "Fisher-index",         place(se) color(black)  justification(right) margin(zero) size(large))          ///
	name(Graph15, replace)
*******
**
** Figure 8: The components of investment-specific technical progres: Decomposition relative to services TFP
**
*******	
generate VA_SERV_X_P_MOD=calA_X_I_BU/TFPVA_SERV_I
generate A_X_TFPVA_SERV_I=A_X_I/TFPVA_SERV_I
local var1 VA_SERV_X_P_MOD
local var2 A_X_TFPVA_SERV_I	
local xcor1=year[70]
local ycor1=`var1'[70]
local xcor2=year[55]
local ycor2=`var2'[55]
graph twoway ///
    (tsline `var1' `var2', lwidth(medthick medthick) lpattern(solid dash)), ///
	graphregion(margin(tiny)) plotregion(margin(zero))                                        ///
	ylabel(0.711 0.843 1.000 1.186 1.406 1.667 , format(%9.1fc) labsize(large))  yscale(range(0.711 1.667001) log)       ///
	ytitle("")                                                                        ///
	xlabel(1950(10)2010, labsize(large)) tmtick(1947(1)2017)                                          ///
	xtitle("")                                                                        ///
	legend(off )                                                                      ///
	subtitle("1947=1, log scale", size(large))                                                   ///
	name(Graph13, replace)
*******
**
** Figure 8: The components of investment-specific technical progres: Decomposition relative to consumption TFP
**
*******	
generate VA_C_X_P_MOD=calA_X_I_BU/calA_C_I
generate A_X_calA_C_I=A_X_I/calA_C_I
local var1 VA_C_X_P_MOD
local var2 A_X_calA_C_I
local xcor1=year[70]
local ycor1=`var1'[70]
local xcor2=year[55]
local ycor2=`var2'[55]
graph twoway ///
    (tsline `var1' `var2', lwidth(medthick medthick) lpattern(solid dash)), ///
	graphregion(margin(tiny)) plotregion(margin(zero))                                        ///
	ylabel(0.711 0.843 1.000 1.186 1.406 1.667, format(%9.1fc) labsize(large))  yscale(range(0.711 1.667001) log)       ///
	ytitle("")                                                                        ///
	xlabel(1950(10)2010, labsize(large)) tmtick(1947(1)2017)                                          ///
	xtitle("")                                                                        ///
	legend(off )                                                                      ///
	subtitle("1947=1, log scale", size(large))                                                   ///
	name(Graph14, replace)	


save "_Data-MainDatabase3.dta", replace


cd "../"


forvalues i=1(1)9 {
	graph export  "InvestmentGraph0`i'.pdf", name(Graph0`i') as(pdf) replace
	graph export  "InvestmentGraph0`i'.ps",  name(Graph0`i') fontface("Times New Roman") as(eps) replace
}

forvalues i=10(1)12 {
	graph export  "InvestmentGraph`i'.pdf", name(Graph`i') as(pdf) replace
	graph export  "InvestmentGraph`i'.ps",  name(Graph`i') fontface("Times New Roman") as(eps) replace
}

graph export  "InvestmentGraph07a.pdf", name(Graph07a) as(pdf) replace 
graph export  "InvestmentGraph08a.pdf", name(Graph08a) as(pdf) replace 
graph export  "InvestmentGraph09a.pdf", name(Graph09a) as(pdf) replace 
graph export  "InvestmentGraph12a.pdf", name(Graph12a) as(pdf) replace 

graph export  "InvestmentGraph07a.ps", name(Graph07a) fontface("Times New Roman") as(eps) replace 
graph export  "InvestmentGraph08a.ps", name(Graph08a) fontface("Times New Roman") as(eps) replace 
graph export  "InvestmentGraph09a.ps", name(Graph09a) fontface("Times New Roman") as(eps) replace 
graph export  "InvestmentGraph12a.ps", name(Graph12a) fontface("Times New Roman") as(eps) replace 

graph export  "InvestmentGraph13.eps", name(Graph13)  fontface("Times New Roman") as(eps) replace 
graph export  "InvestmentGraph14.eps", name(Graph14)  fontface("Times New Roman") as(eps) replace 

graph export  "InvestmentGraph15.ps", name(Graph15)   fontface("Times New Roman") as(eps) replace 


graph export  "InvestmentGraph15.pdf", name(Graph15) as(pdf) replace 



exit





