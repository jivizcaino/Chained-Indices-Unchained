cls 
clear all
* Set format to show at least 5 decimal places
cd "\Users\Nacho\Dropbox\Omar&Nacho\Baqaee_Burstein\Note on Structural Change\"

use "C:\Users\Nacho\Dropbox\Omar&Nacho\Baqaee_Burstein\Note on Structural Change\HRV2021_Data_Appendix\Supplementary\paper\Data-MainDatabase.dta" 
 
foreach var of varlist _all {
    capture confirm numeric variable `var'
    if !_rc {
        recast float `var'   // reduce precision
        format `var' %9.0g   // typical float format
    }
}
export excel using "HRV_Data_float.xlsx", replace firstrow(variables)

* Store variable names
ds, has(varlabel)
local vars `r(varlist)'
local nvars: word count `vars'

* Create a macro list of variable names and their labels
local i = 1
foreach v of local vars {
    local name_`i' "`v'"
    local label_`i' "`: variable label `v''"
    local ++i
}

* Now clear and create a dataset with these
clear
set obs `nvars'

gen str32 varname = ""
gen str120 varlabel = ""

* Fill in the new dataset
forvalues j = 1/`nvars' {
    replace varname = "`name_`j''" in `j'
    replace varlabel = "`label_`j''" in `j'
}

* Export to a new sheet
export excel using "HRV_Data_NEW.xlsx", sheet("var_labels") sheetreplace


