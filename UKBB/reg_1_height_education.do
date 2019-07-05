//Neil Davies 22/11/18
//This runs the ivreg analysis for the height->education

//Define programs
//Run the four regressions (standard OLS, standard MR with allele score on full sample, standard MR on one of the siblings
// and within family MR.


cap prog drop regressions
prog def regressions
args outcome exposure instrument
preserve
if "`outcome'"!="`exposure'"{

	reg `outcome' `exposure' pc1-pc20 cov_age cov_male if n_sibs>1 & n_sibs<8 & n_sibs!=. ,ro cluster(famid)
	regsave `exposure' using "results/`outcome'_`exposure'",detail(all) pval replace

	xtreg `outcome' `exposure' pc1-pc20 cov_age cov_male if n_sibs>1 & n_sibs<8 & n_sibs!=.,ro cluster(famid) fe i(famid)
	regsave `exposure' using "results/`outcome'_`exposure'",detail(all) pval append
	
	ivreg2 `outcome' (`exposure'=`instrument') pc1-pc20  cov_age cov_male if n_sibs==1| n_sibs>8,ro cluster(famid)
	regsave `exposure' using "results/`outcome'_`exposure'",detail(all) pval append
	
	keep if `outcome'!=.
	bys famid:gen n_sibs2=_N
	drop if  n_sibs2==1	
			
	ivreg2 `outcome' (`exposure'=`instrument') pc1-pc20  cov_age cov_male if n_sibs>1 & n_sibs<8 ,ro  cluster(famid)
	regsave `exposure' using "results/`outcome'_`exposure'",detail(all) pval append

	xtivreg `outcome' (`exposure'=`instrument') pc1-pc20  cov_age cov_male if n_sibs>1 & n_sibs<8 ,i(famid) fe  vce(cluster famid )
	regsave `exposure' using "results/`outcome'_`exposure'",detail(all) pval append
	}
else{
}
end


//1. Open data and run the IPD analysis in UKBB

use "workingdata/analysisdata", clear

//Set famid = to unique value for the individuals not in families
replace famid =_n+22658 if famid ==.

//Run regressions across all the exposures and outcomes
ds eduyears3 out_highbloodpressure2 out_diabetes
foreach i in `r(varlist)'{
	regressions `i' out_bmi locke_69_wprs
	regressions `i' out_height wood_386_wprs
	}

//2. Clean results and meta-analyse across UKBB and HUNT
//Clean the HUNT non-family results
//The OLS results not allowing for a family FE
/*
The following files were sent by Ben from the analysis on HUNT data
#OLS
lm()
r_non_within_ols_v7

#OLS family FE
plm()
r_within_ols_v7.txt

#Standard MR (non-sibs)
ivreg()
r_non_within_v8.txt

#Standard MR (sibs)
ivreg()
r_non_within_v7.txt

# MR within families
plm()
r_within_v7.txt
*/
//Clean each of the files and create a single results file.
foreach i in r_non_within_ols_v7 r_within_ols_v7 r_non_within_v8 r_non_within_sibs_v7  r_within_v7{
	import delimited results/`i'.txt, delimiter(space) clear 
	drop in 1
	ds *
	cap: rename v1 exposure 
	cap: rename v2 outcome
	cap: rename v3 N 
	cap: rename v4 coef
	cap: rename v5 stderr
	cap: drop v6
	
	cap: rename n N 
	cap: rename beta coef
	cap: rename se stderr
	cap: drop pvalue
	
	local type=`type'+1
	gen type=`type'
	gen study=2
	save "workingdata/hunt_`i'",replace
	count
	}
	
foreach i in r_non_within_ols_v7 r_within_ols_v7 r_non_within_v8 r_non_within_sibs_v7 {
	append using "workingdata/hunt_`i'",
	}
compress
save "workingdata/hunt_ipd",replace
	
//Clean UK Biobank results so that they're in the same format as the HUNT results.
#delimit ;
use "results/out_bmi_eduyears3.dta",clear;
foreach i in 
eduyears3 out_bmi out_highbloodpressure2 out_height out_diabetes{;
	cap:append using "results/`i'_eduyears3";
	cap:	append using "results/`i'_out_bmi";
	cap:	append using "results/`i'_out_height";
	};
#delimit cr
rename var exposure
rename depvar outcome
order outcome exposure

duplicates drop
drop if exposure==outcome

gen type=1 if cmd=="regress"
replace type=2 if cmd=="xtreg"
replace type=3 if cmd=="ivreg2" & N>40000
replace type=4 if cmd=="ivreg2" & N<40000
replace type=5 if cmd=="xtivreg"

sort  outcome exposure type

//Add index for outcome
levels outcome
gen outcome_index=.
foreach i in `r(levels)'{
	local x=`x'+1
	replace outcome_index=`x' if outcome=="`i'"
	}
order outcome exposure type outcome_index
//Input all the OLS estimates into a matrix for one exposure
drop if outcome=="out_sys_bp"| outcome=="out_dia_bp"| outcome=="cov_male"| outcome=="cov_age"

gen cont=0
foreach i in eduyears3 out_sys_bp out_dia_bp out_height out_bmi out_dia_bp cov_num_sisters cov_num_brothers cov_comp_bodysize8 cov_comp_height8 cov_birthweight{
	replace cont=1 if outcome=="`i'"
	}
drop if (outcome=="out_height" & exposure=="out_bmi")|(exposure=="out_height" & outcome=="out_bmi")|(exposure=="out_height" & outcome=="out_highbloodpressure2")|(exposure=="out_height" & outcome=="out_diabetes")
drop if exposure=="eduyears3"
//Generate index per analysis:
gen index=.
replace index=round((_n+1.9)/5)
order index

keep index outcome exposure type coef stderr pval N N_g cmd cmdline cont

save "workingdata/ukbb_ipd_results",replace

append using "workingdata/hunt_ipd"
replace study=1 if study==.

drop if outcome=="out_diabetes" & exposure=="out_height"
drop if outcome=="out_highbloodpressure" & exposure=="out_height"

//Meta-analyse the 5 IPD analyses:
//Setup results variables
gen m_index=.
gen m_outcome=""
gen m_exposure=""
gen m_study_outcome=""
gen m_study_exposure=""
gen m_type=. 
gen m_outcome_index=.
gen m_coef=.
gen m_stderr=.
gen m_pval=.
gen m_N=.
gen m_het_chi=.
gen m_het_pvalue=.

cap prog drop save_results
prog def save_results
args index exposure outcome study_outcome study_exposure type outcome_index coef stderr pval row het_chi het_pvalue
replace m_index=`index' in `row'
replace m_outcome="`outcome'" in `row'
replace m_exposure="`exposure'" in `row'
replace m_study_outcome="`study_outcome'" in `row'
replace m_study_exposure="`study_exposure'" in `row'
replace m_type=`type' in `row'
replace m_outcome_index=`outcome_index' in `row'
replace m_coef=`coef' in `row'
replace m_stderr=`stderr' in `row'
replace m_pval=`pval' in `row'
replace m_het_chi=`het_chi' in `row'
replace m_het_pvalue=`het_pvalue' in `row'
end

//Generate the counts for each analysis
replace outcome="eduyears3" if outcome=="eduyears2"
replace outcome="out_highbloodpressure2" if outcome=="out_highbloodpressure"
bys outcome exposure type: egen total_n=sum(N)

//BMI on education
forvalues i=1(1)5{
	local j=`j'+1
	di "`j'"
	metan coef stderr if substr(outcome,1,8)=="eduyears" & exposure=="out_bmi" & type==`i'
	save_results 8 out_bmi eduyears same same `i' 1 r(ES) r(seES) r(p_z) `j' r(het) r(p_het)
	tabstat total_n if substr(outcome,1,8)=="eduyears" & exposure=="out_bmi" & type==`i',save
	replace m_N=el(r(StatTotal),1,1) in `j'
	}

//Height on education
forvalues i=1(1)5{
	local j=`j'+1
	di "`j'"
	metan coef stderr if substr(outcome,1,8)=="eduyears" & exposure=="out_height" & type==`i'
	save_results 8 out_height eduyears same same `i' 1 r(ES) r(seES) r(p_z) `j' r(het) r(p_het)
	tabstat total_n if substr(outcome,1,8)=="eduyears" & exposure=="out_height" & type==`i',save
	replace m_N=el(r(StatTotal),1,1) in `j'
	}
	
//BMI on high blood pressure
forvalues i=1(1)5{
	local j=`j'+1
	di "`j'"
	metan coef stderr if substr(outcome,1,8)=="out_diab" & exposure=="out_bmi" & type==`i'
	save_results 8 out_bmi out_diabetes  same same `i' 1 r(ES) r(seES) r(p_z) `j' r(het) r(p_het)
	tabstat total_n if substr(outcome,1,8)=="out_diab" & exposure=="out_bmi" & type==`i',save
	replace m_N=el(r(StatTotal),1,1) in `j'
	}

//BMI on diabetes
forvalues i=1(1)5{
	local j=`j'+1
	di "`j'"
	metan coef stderr if substr(outcome,1,8)=="out_high" & exposure=="out_bmi" & type==`i'
	save_results 8 out_bmi out_highbloodpressure2  same same `i' 1 r(ES) r(seES) r(p_z) `j' r(het) r(p_het)
	tabstat total_n if substr(outcome,1,8)=="out_high" & exposure=="out_bmi" & type==`i',save
	replace m_N=el(r(StatTotal),1,1) in `j'
	}
keep m_*
drop if m_type==.
compress

rename m_* *
save "workingdata/meta_analysed_ipd",replace

//Clean results for inclusion in text

//gen confidence intervals
gen beta=coef
format %9.2f beta 
replace beta=beta*100 if outcome=="out_diabetes"|outcome=="out_highbloodpressure2"
gen lci=beta-1.96*stderr if !(outcome=="out_diabetes"|outcome=="out_highbloodpressure2")
gen uci=beta+1.96*stderr if !(outcome=="out_diabetes"|outcome=="out_highbloodpressure2")
replace lci=beta-1.96*stderr*100 if outcome=="out_diabetes"|outcome=="out_highbloodpressure2"
replace uci=beta+1.96*stderr*100 if outcome=="out_diabetes"|outcome=="out_highbloodpressure2"
gen double pval2=pval
replace  pval2=2*normal(-abs(coef/stderr)) 
format 
