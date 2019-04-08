//Neil Davies 05/12/18
//This combines the data from R plm and Gib's function for each of the outcome to generate a plot:

import delimited results/r_within.txt, delimiter(space) stripquote(yes) clear 
rename outcome exposurex
rename exposure outcome
rename exposure exposure
gen method="plm_within"
save "workingdata/plm_within",replace

import delimited results/r_non_within.txt, delimiter(space) stripquote(yes) clear 
gen method="standard_mr"
rename outcome exposurex
rename exposure outcome
rename exposure exposure
save "workingdata/standard_mr",replace

import delimited results/results_do_mr_wf.txt, delimiter(space) stripquote(yes) clear 
gen method="do_mr_wf"
save "workingdata/do_mr_wf",replace

rename b beta
rename pval pvalue
append using "workingdata/standard_mr",
append using "workingdata/plm_within",

foreach i in eduyears2 out_highbloodpressure out_diabetes out_bmi out_height{
	cap:replace outcome="`i'" if outcome=="`i'_res"
	cap:replace exposure="`i'" if exposure=="`i'_res"
	}

save "workingdata/all_results_r",replace

//Append all the results from Stata

use "results/ber within                                                                         mi_eduyears2",clear
foreach j in  eduyears2  out_diabetes  out_highbloodpressure {
	foreach i in bmi height{
		append using "results/`i'_`j'",
		}
	}
	
duplicates drop
order  depvar var cmd

gen order=round((_n+2)/5)
bys order:gen within=_n
order order within

rename depvar outcome
rename var exposure

rename coef beta
rename stderr se 
rename pval pvalue

keep order-N F cdf

gen method="OLS" if within==5
replace method="Family FE" if within==4
replace method="Standard MR full sample" if within==3
replace method="Standard MR sibling sample" if within==2
replace method="Within family MR FE" if within==1

append using "workingdata/all_results_r"

replace within =6-within



sort exposure outcome within method


drop order 
gen order=round((_n+3.99)/8)
replace within =6 if method=="standard_mr"
replace within =7 if method=="plm_within"
replace within =8 if method=="do_mr_wf"


order order within
sort ord
