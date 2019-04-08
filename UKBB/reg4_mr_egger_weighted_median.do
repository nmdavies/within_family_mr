//Neil Davies 01/04/19
//This runs the MR-Egger, weighted median and mode estimators for the cross sample results:

use "workingdata/ukbb_out_highbloodpressure2",clear
foreach i in eduyears2 out_bmi out_diabetes out_height{
	joinby SNP using  "workingdata/ukbb_`i'"
	}
keep SNP UKBB_EAF UKBB_effect_allele UKBB_other_allele sample_size ukbb_*
compress
save "workingdata/ukbb_results",replace	

joinby SNP using "workingdata/hunt_results",unmatched(master) _merge(SSS)

//Harmonize to the UKBB effect alleles
forvalues i=1(1)464{
	local ukbb_effect_allele=UKBB_effect_allele[`i']
	local ukbb_other_allele=UKBB_other_allele[`i']
	local hunt_effect_allele=HUNT_effect_allele[`i']
	local hunt_other_allele=HUNT_other_allele[`i']
	
	di "UKBB effect=`ukbb_effect_allele' other=`ukbb_other_allele' HUNT effect=`hunt_effect_allele' other=`hunt_other_allele'"
	
	if "`ukbb_effect_allele'"=="`hunt_other_allele'" & "`ukbb_other_allele'"=="`hunt_effect_allele'"{
		foreach j in HUNT_eduyears2_beta HUNT_out_bmi_beta HUNT_out_diabetes_beta HUNT_out_height_beta{
			replace `j'=`j'*-1 in `i'
			}
		di "Flipped"
		}
	else{
		if "`ukbb_effect_allele'"=="`hunt_effect_allele'" & "`ukbb_other_allele'"=="`hunt_other_allele'"{
			di "All ok"
			}
		else{
			di "ERROR UNMATCHED"
			}
		}
	}

gen outcome=""
gen exposure=""
gen study_outcome=""
gen study_exposure=""
gen beta=.
gen se=.
gen pval=.
gen n_snps=.	
	
//Run MR regression
mregger HUNT_eduyears2_beta HUNT_out_height_beta [aw=1/(HUNT_eduyears2_se^2)] if wood==1,
regsave using "results/summary_data_mr_egger_weighted_median_mode",replace detail(all) pval ci

rename HUNT_out_highbloodpressure_beta HUNT_out_hbp2_beta
rename HUNT_out_highbloodpressure_se HUNT_out_hbp2_se

rename ukbb_out_highbloodpressure2_beta ukbb_out_hbp2_beta
rename ukbb_out_highbloodpressure2_se ukbb_out_hbp2_se

foreach i in HUNT ukbb{
	foreach j in HUNT ukbb{ 
		foreach k in eduyears2 out_hbp2 out_diabetes{
			mregger `i'_`k'_beta `j'_out_height_beta [aw=1/(`i'_`k'_se^2)] if wood==1,
			regsave using "results/summary_data_mr_egger_weighted_median_mode",append detail(all) pval ci
			mregger `i'_`k'_beta `j'_out_bmi_beta [aw=1/(`i'_`k'_se^2)] if locke==1,
			regsave using "results/summary_data_mr_egger_weighted_median_mode",append detail(all) pval ci
			}
		}
	}
	
foreach i in HUNT ukbb{
	foreach j in HUNT ukbb{ 
		foreach k in eduyears2 out_hbp2 out_diabetes{
			mrmodal `i'_`k'_beta `i'_`k'_se `j'_out_height_beta `j'_out_height_se if wood==1,
			regsave using "results/summary_data_mr_egger_weighted_median_mode",append detail(all) pval ci
			mrmodal `i'_`k'_beta `i'_`k'_se `j'_out_bmi_beta `j'_out_bmi_se if locke==1,
			regsave using "results/summary_data_mr_egger_weighted_median_mode",append detail(all) pval ci
			}
		}
	}
	
foreach i in HUNT ukbb{
	foreach j in HUNT ukbb{ 
		foreach k in eduyears2 out_hbp2 out_diabetes{
			mrmedian `i'_`k'_beta `i'_`k'_se `j'_out_height_beta `j'_out_height_se if wood==1,
			regsave using "results/summary_data_mr_egger_weighted_median_mode",append detail(all) pval ci
			mrmedian `i'_`k'_beta `i'_`k'_se `j'_out_bmi_beta `j'_out_bmi_se if locke==1,
			regsave using "results/summary_data_mr_egger_weighted_median_mode",append detail(all) pval ci
			}
		}
	}	

	
//Clean the results	
use results/summary_data_mr_egger_weighted_median_mode.dta,clear
duplicates drop
drop if strpos(cmdline,"HUNT" )==0| strpos(cmdline,"ukbb" )==0
drop if strpos(cmdline,"height")!=0 & strpos(cmdline,"diabetes")!=0
drop if strpos(cmdline,"height")!=0 & strpos(cmdline,"hbp2")!=0

//Meta-analyse the two cross estimates

rename coef Coef
rename stderr Stderr
rename pval Pval

gen index=.
gen outcome=""
gen exposure=""
gen study_outcome=""
gen study_exposure=""
gen type=. 
gen outcome_index=.
gen coef=.
gen stderr=.
gen pval=.
gen N=.
gen het_chi=.
gen het_pvalue=.
gen method=""

cap prog drop save_results
prog def save_results
args index exposure outcome study_outcome study_exposure type outcome_index coef stderr pval row het_chi het_pvalue method
replace index=`row' in `row'
replace outcome="`outcome'" in `row'
replace exposure="`exposure'" in `row'
replace study_outcome="`study_outcome'" in `row'
replace study_exposure="`study_exposure'" in `row'
replace type=`type' in `row'
replace outcome_index=`row' in `row'
replace coef=`coef' in `row'
replace stderr=`stderr' in `row'
replace pval=`pval' in `row'
replace het_chi=`het_chi' in `row'
replace het_pvalue=`het_pvalue' in `row'
replace method="`method'" in `row'
end

//Height on education
metan Coef Stderr if strpos(cmdline,"height")!=0 & strpos(cmdline,"eduyears")!=0 & strpos(cmdline,"mrmedian")!=0
save_results 1 out_height eduyears3 cross cross 5 1 r(ES) r(seES) r(p_z) 2 r(het) r(p_het) mrmedian

metan Coef Stderr if strpos(cmdline,"height")!=0 & strpos(cmdline,"eduyears")!=0 & strpos(cmdline,"mrmodal")!=0
save_results 2 out_height eduyears3 cross cross 5 2 r(ES) r(seES) r(p_z) 3 r(het) r(p_het) mrmodal

metan Coef Stderr if strpos(cmdline,"height")!=0 & strpos(cmdline,"eduyears")!=0 & strpos(cmdline,"mregger")!=0 & strpos(var,"slope")!=0
save_results 3 out_height eduyears3 cross cross 5 3 r(ES) r(seES) r(p_z) 4 r(het) r(p_het) mregger_slope

metan Coef Stderr if strpos(cmdline,"height")!=0 & strpos(cmdline,"eduyears")!=0 & strpos(cmdline,"mregger")!=0 & strpos(var,"_cons")!=0
save_results 4 out_height eduyears3 cross cross 5 4 r(ES) r(seES) r(p_z) 5 r(het) r(p_het) mregger_cons

//BMI on education 
metan Coef Stderr if strpos(cmdline,"bmi")!=0 & strpos(cmdline,"eduyears")!=0 & strpos(cmdline,"mrmedian")!=0
save_results 1 out_bmi eduyears3 cross cross 5 5 r(ES) r(seES) r(p_z) 7 r(het) r(p_het) mrmedian

metan Coef Stderr if strpos(cmdline,"bmi")!=0 & strpos(cmdline,"eduyears")!=0 & strpos(cmdline,"mrmodal")!=0
save_results 2 out_bmi eduyears3 cross cross 5 6 r(ES) r(seES) r(p_z) 8 r(het) r(p_het) mrmodal

metan Coef Stderr if strpos(cmdline,"bmi")!=0 & strpos(cmdline,"eduyears")!=0 & strpos(cmdline,"mregger")!=0 & strpos(var,"slope")!=0
save_results 3 out_bmi eduyears3 cross cross 5 7 r(ES) r(seES) r(p_z) 9 r(het) r(p_het) mregger_slope

metan Coef Stderr if strpos(cmdline,"bmi")!=0 & strpos(cmdline,"eduyears")!=0 & strpos(cmdline,"mregger")!=0 & strpos(var,"_cons")!=0
save_results 4 out_bmi eduyears3 cross cross 5 8 r(ES) r(seES) r(p_z) 10 r(het) r(p_het) mregger_cons

//BMI on diabetes 
metan Coef Stderr if strpos(cmdline,"bmi")!=0 & strpos(cmdline,"diabetes")!=0 & strpos(cmdline,"mrmedian")!=0
save_results 1 out_bmi diabetes cross cross 5 9 r(ES) r(seES) r(p_z) 12 r(het) r(p_het) mrmedian

metan Coef Stderr if strpos(cmdline,"bmi")!=0 & strpos(cmdline,"diabetes")!=0 & strpos(cmdline,"mrmodal")!=0
save_results 2 out_bmi diabetes cross cross 5 10 r(ES) r(seES) r(p_z) 13 r(het) r(p_het) mrmodal

metan Coef Stderr if strpos(cmdline,"bmi")!=0 & strpos(cmdline,"diabetes")!=0 & strpos(cmdline,"mregger")!=0 & strpos(var,"slope")!=0
save_results 3 out_bmi diabetes cross cross 5 11 r(ES) r(seES) r(p_z) 14 r(het) r(p_het) mregger_slope

metan Coef Stderr if strpos(cmdline,"bmi")!=0 & strpos(cmdline,"diabetes")!=0 & strpos(cmdline,"mregger")!=0 & strpos(var,"_cons")!=0
save_results 4 out_bmi diabetes cross cross 5 12 r(ES) r(seES) r(p_z) 15 r(het) r(p_het) mregger_cons

//BMI on high blood pressure 
metan Coef Stderr if strpos(cmdline,"bmi")!=0 & strpos(cmdline,"hbp2")!=0 & strpos(cmdline,"mrmedian")!=0
save_results 1 out_bmi hbp2 cross cross 5 13 r(ES) r(seES) r(p_z) 17 r(het) r(p_het) mrmedian

metan Coef Stderr if strpos(cmdline,"bmi")!=0 & strpos(cmdline,"hbp2")!=0 & strpos(cmdline,"mrmodal")!=0
save_results 2 out_bmi hbp2 cross cross 5 14 r(ES) r(seES) r(p_z) 18 r(het) r(p_het) mrmodal

metan Coef Stderr if strpos(cmdline,"bmi")!=0 & strpos(cmdline,"hbp2")!=0 & strpos(cmdline,"mregger")!=0 & strpos(var,"slope")!=0
save_results 3 out_bmi hbp2 cross cross 5 15 r(ES) r(seES) r(p_z) 19 r(het) r(p_het) mregger_slope

metan Coef Stderr if strpos(cmdline,"bmi")!=0 & strpos(cmdline,"hbp2")!=0 & strpos(cmdline,"mregger")!=0 & strpos(var,"_cons")!=0
save_results 4 out_bmi hbp2 cross cross 5 16 r(ES) r(seES) r(p_z) 20 r(het) r(p_het) mregger_cons

keep outcome-method
compress
drop if outcome==""
save workingdata/mregger_results,replace

use "workingdata/metan",clear

keep if study_exposure  =="cross"

append using  workingdata/mregger_results,
replace outcome_index=6 in 2
replace outcome_index=16 in 3
replace outcome_index=11 in 4

sort outcome_index
replace method="ivw" if method==""

gen uci=coef+1.96*stderr
gen lci=coef-1.96*stderr

replace outcome="diabetes" if outcome=="out_diabetes"
replace outcome="hbp2" if outcome=="out_highbloodpressure2"

//Rescale binary outcomes
foreach i in coef uci lci{
	replace `i'=`i'*100 if outcome=="diabetes"|outcome=="hbp2"
	}
save workingdata/mregger_results_for_r,replace
drop if method=="mregger_cons"
save workingdata/mregger_results_for_r_no_cons,replace
