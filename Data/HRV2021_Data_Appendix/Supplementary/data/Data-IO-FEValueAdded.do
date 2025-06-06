version 16.0
capture log close
clear all
set more off
global DataPath `"/Users/akos/Dropbox/AkosBertRichard/Investment/Submissions/FinalSubmission/Supplementary"'
cd "$DataPath/data"
set type float
***
***  Constructing consumption and investment value added first by calcuating the total requirement matrix
***
forvalues year=1947/1962 {
	use "Data-IO-TablesMake1947-1962", clear
	keep if year==`year'
	mkmat Commodity* if _n<_N,  matrix(M) rownames(naics)
	matrix UnitInd=J(1,rowsof(M),1)
	matrix UnitCom=J(colsof(M),1,1)
	matrix q=(UnitInd*M)'
	matrix g=M*UnitCom
	ds, has(varl *Scrap*)
	generate scrap=`r(varlist)'/TotalIndustryOutput
	mkmat scrap if _n<_N, matrix(p) rownames(naics)
	mkmat Commodity* if _n<=rowsof(M),  matrix(M) rownames(naics)
	use "Data-IO-TablesUse1947-1962", clear
	drop Industry
	keep if year==`year'
	mkmat Industry* if _n<=colsof(M), matrix(U) rownames(naics)
	mkmat Final*    if _n<=colsof(M), matrix(E) rownames(naics)
	matrix UnitFin=J(colsof(E),1,1)
	matrix dc=q-U*UnitInd'-E*UnitFin
	svmat dc, names(dc)
	ds, has(varl *Personal*)
	replace `r(varlist)'=`r(varlist)'+dc1
	matrix va=g-(UnitCom'*U)'
	matrix vaM=diag(inv(diag(g))*va)
	matrix B=U*inv(diag(g))
	matrix D=M*inv(diag(q))
	matrix W=D
	matrix TRM=W*inv(I(colsof(B*W))-B*W)
	generate PCE=FinalUse1
	generate PIS=FinalUse2
	generate PIE=FinalUse3
	generate PIP=FinalUse4
	generate INI=FinalUse5
	generate GCE=FinalUse8
	generate GIS=FinalUse9
	generate GIE=FinalUse10
	generate GIP=FinalUse11
	generate EX=FinalUse6
	generate IM=FinalUse7
	foreach var in PCE PIE PIS PIP INI GCE GIE GIS GIP EX IM {
		mkmat `var' if _n<=rowsof(E), matrix(`var') rownames(naics)
	}
	foreach var in PCE PIE PIS PIP INI GCE GIE GIS GIP EX IM {
		matrix `var'_va=vaM*TRM*`var'
	}	
	matrix define V`year'=(PCE_va, PIE_va, PIS_va, PIP_va, INI_va, GCE_va, GIE_va, GIS_va, GIP_va, EX_va, IM_va)
	matrix drop M U q g p dc B D W PCE PIE PIS PIP INI GCE GIE GIS EX IM PCE_va GCE_va GIE_va GIS_va GIP_va PIE_va PIS_va PIP_va INI_va EX_va IM_va
}
***
***
***
forvalues year=1963/1996 {
	use "Data-IO-TablesMake1963-1996", clear
	keep if year==`year'
	mkmat Commodity* if _n<_N,  matrix(M) rownames(naics)
	matrix UnitInd=J(1,rowsof(M),1)
	matrix UnitCom=J(colsof(M),1,1)
	matrix q=(UnitInd*M)'
	matrix g=M*UnitCom
	ds, has(varl *Scrap*)
	generate scrap=`r(varlist)'/TotalIndustryOutput
	mkmat scrap if _n<_N, matrix(p) rownames(naics)
	mkmat Commodity* if _n<=rowsof(M),  matrix(M) rownames(naics)
	use "Data-IO-TablesUse1963-1996", clear
	drop Industry
	keep if year==`year'
	mkmat Industry* if _n<=colsof(M), matrix(U) rownames(naics)
	mkmat Final*    if _n<=colsof(M), matrix(E) rownames(naics)
	matrix UnitFin=J(colsof(E),1,1)
	matrix dc=q-U*UnitInd'-E*UnitFin
	svmat dc, names(dc)
	ds, has(varl *Personal*)
	replace `r(varlist)'=`r(varlist)'+dc1
	matrix va=g-(UnitCom'*U)'
	matrix vaM=diag(inv(diag(g))*va)
	matrix B=U*inv(diag(g))
	matrix D=M*inv(diag(q))
	matrix W=D
	matrix TRM=W*inv(I(colsof(B*W))-B*W)
	generate PCE=FinalUse1
	generate PIS=FinalUse2
	generate PIE=FinalUse3
	generate PIP=FinalUse4
	generate INI=FinalUse5
	generate GCE=FinalUse8
	generate GIS=FinalUse9
	generate GIE=FinalUse10
	generate GIP=FinalUse11
	generate EX=FinalUse6
	generate IM=FinalUse7
	foreach var in PCE PIE PIS PIP INI GCE GIE GIS GIP EX IM {
		mkmat `var' if _n<=rowsof(E), matrix(`var') rownames(naics)
	}
	foreach var in PCE PIE PIS PIP INI GCE GIE GIS GIP EX IM {
		matrix `var'_va=vaM*TRM*`var'
	}	
	matrix define V`year'=(PCE_va, PIE_va, PIS_va, PIP_va, INI_va, GCE_va, GIE_va, GIS_va, GIP_va, EX_va, IM_va)
	matrix drop M U q g p dc B D W PCE PIE PIS PIP INI GCE GIE GIS EX IM PCE_va GCE_va GIE_va GIS_va GIP_va PIE_va PIS_va PIP_va INI_va EX_va IM_va
}
***
***
***
forvalues year=1997/2017 {
	use "Data-IO-TablesMake1997-2017", clear
	keep if year==`year'
	mkmat Commodity* if _n<_N,  matrix(M) rownames(naics)
	matrix UnitInd=J(1,rowsof(M),1)
	matrix UnitCom=J(colsof(M),1,1)
	matrix q=(UnitInd*M)'
	matrix g=M*UnitCom
	ds, has(varl *Scrap*)
	generate scrap=`r(varlist)'/TotalIndustryOutput
	mkmat scrap if _n<_N, matrix(p) rownames(naics)
	mkmat Commodity* if _n<=rowsof(M),  matrix(M) rownames(naics)
	use "Data-IO-TablesUse1997-2017", clear
	drop Industry
	keep if year==`year'
	mkmat Industry* if _n<=colsof(M), matrix(U) rownames(naics)
	mkmat Final*    if _n<=colsof(M), matrix(E) rownames(naics)
	matrix UnitFin=J(colsof(E),1,1)
	matrix dc=q-U*UnitInd'-E*UnitFin
	svmat dc, names(dc)
	ds, has(varl *Personal*)
	replace `r(varlist)'=`r(varlist)'+dc1
	matrix va=g-(UnitCom'*U)'
	matrix vaM=diag(inv(diag(g))*va)
	matrix B=U*inv(diag(g))
	matrix D=M*inv(diag(q))
	matrix W=D
	matrix TRM=W*inv(I(colsof(B*W))-B*W)
	generate PCE=FinalUse1
	generate PIS=FinalUse2+FinalUse5
	generate PIE=FinalUse3
	generate PIP=FinalUse4
	generate INI=FinalUse6
	generate GCE=FinalUse9
	generate GIS=FinalUse10
	generate GIE=FinalUse11
	generate GIP=FinalUse12
	generate EX=FinalUse7
	generate IM=FinalUse8
	foreach var in PCE PIE PIS PIP INI GCE GIE GIS GIP EX IM {
		mkmat `var' if _n<=rowsof(E), matrix(`var') rownames(naics)
	}
	foreach var in PCE PIE PIS PIP INI GCE GIE GIS GIP EX IM {
		matrix `var'_va=vaM*TRM*`var'
	}	
	matrix define V`year'=(PCE_va, PIE_va, PIS_va, PIP_va, INI_va, GCE_va, GIE_va, GIS_va, GIP_va, EX_va, IM_va)
	matrix drop M U q g p dc B D W PCE PIE PIS PIP INI GCE GIE GIS EX IM PCE_va GCE_va GIE_va GIS_va GIP_va PIE_va PIS_va PIP_va INI_va EX_va IM_va
}
**
**
**
forvalues year=1947/2017 {
	drop _all
	svmat2 V`year', names(col) rnames(naics)
	generate year=`year'
	order naics year, before(PCE)
	label variable PCE  "Personal consumption expenditures, gross value added requirement"
	label variable PIE  "Private investment in equipment, gross value added requirement"
	label variable PIS  "Private investment in structures, gross value added requirement"
	label variable PIP  "Private investment in intellectual property products, gross value added requirement"
	label variable INI  "Inventory investment, gross value added requirement"
	label variable GCE  "Government consumption expenditures, gross value added requirement"
	label variable GIE  "Government investment in equipment, gross value added requirement"
	label variable GIS  "Government investment in structures, gross value added requirement"
	label variable GIP  "Government investment in intellectual property products, gross value added requirement"
	label variable EX   "Exports gross value added requirement"
	label variable IM   "Imports gross value added requirement"	
	foreach var in PCE PIE PIS PIP INI GCE GIE GIS GIP EX IM {
		rename `var' VA`var'
	}
	generate id=_n
	generate VAGDP=VAPCE+VAPIE+VAPIS+VAPIP+VAINI+VAGCE+VAGIE+VAGIS+VAGIP+VAEX+VAIM
	label variable VAGDP   "GDP, gross value added requirement"	
	tempfile V`year'
	save `V`year'', replace
}

*******
***
*** Final Expenditure value added  
***
*******
** 1947-1962
*******
use `V1947', clear
forvalues year=1948/1962 {
	append using `V`year''
}
drop id
save "Data-IO-FEValueAdded1947-1962.dta", replace
*******
** 1963-1997
*******
use `V1963', clear
forvalues year=1964/1996 {
	append using `V`year''
}
drop id
save "Data-IO-FEValueAdded1963-1996.dta", replace
*******
** 1997-2017
*******
use `V1997', clear
forvalues year=1998/2017 {
	append using `V`year''
}
drop id
save "Data-IO-FEValueAdded1997-2017.dta", replace






exit
