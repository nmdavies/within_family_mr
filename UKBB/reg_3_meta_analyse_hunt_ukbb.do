//Neil Davies 21/12/18
//This imports the results from HUNT and meta-analyses with the UK BB data


cd "$path1/19_within_families_paper/"

import delimited results/individual_snps_v8_withinfo.txt, delimiter(space) clear
keep snp effect_alleledose other_alleledose
duplicates drop
save "workingdata/hunt_effect_allele_other",replace

import delimited results/individual_snps_v8_withinfo.txt, delimiter(space) clear
drop if trait=="0"
joinby snp using "workingdata/hunt_effect_allele_other",unmatched(master)
drop _m

rename snp SNP
rename eaf EAF

//Convert to wide format to use mrrobust package
levels trait 
foreach trait in `r(levels)'{
	preserve
	keep if trait =="`trait'"
	rename beta `trait'_beta
	rename se `trait'_se
	save "workingdata/hunt_`trait'",replace
	restore
	}
use "workingdata/hunt_out_highbloodpressure",clear
foreach i in eduyears2 out_bmi out_diabetes out_height{
	joinby SNP using  "workingdata/hunt_`i'"
	}
keep SNP EAF effect_allele other_allele sample_size-out_height_se out_highbloodpressure_beta out_highbloodpressure_se
compress
save "workingdata/hunt_results",replace	

//Identify the BMI and height variants
//69 Locke SNPs
use rsid using "workingdata/locke_69_coefficients.dta",clear
rename rsid SNP
save  "workingdata/locke_69_rsid.dta",replace

//386 Wood SNPs
use v2 using "workingdata/wood_386_coefficients.dta",clear
rename v2 SNP
save  "workingdata/wood_386_rsid.dta",replace

//Create indicator for Locke and Wood SNPs
use "workingdata/hunt_results", clear
joinby SNP using "workingdata/locke_69_rsid.dta",unmatched(master)

gen locke=(_m==3)
drop _m
joinby SNP using "workingdata/wood_386_rsid.dta",unmatched(master)
gen wood=(_m==3)
drop _m

rename out_highbloodpressure_beta-out_height_se HUNT_=

save "workingdata/hunt_results",replace

//Clean UKBB data to the same format:
use "workingdata/indi_snp_analysis_clean.dta",clear
duplicates drop
replace trait="eduyears2" if trait=="eduyears3"
keep if method=="Within FE"

levels trait 
foreach trait in `r(levels)'{
	preserve
	keep if trait =="`trait'"
	rename beta ukbb_`trait'_beta
	rename se ukbb_`trait'_se
	rename effect_allele UKBB_effect_allele
	rename other_allele UKBB_other_allele
	rename EAF UKBB_EAF
	save "workingdata/ukbb_`trait'",replace
	restore
	}
	
use "workingdata/ukbb_out_highbloodpressure2",clear
foreach i in eduyears2 out_bmi out_diabetes out_height{
	joinby SNP using  "workingdata/ukbb_`i'"
	}
keep SNP UKBB_EAF UKBB_effect_allele UKBB_other_allele sample_size ukbb_*
compress
save "workingdata/ukbb_results",replace	

joinby SNP using "workingdata/hunt_results",unmatched(master) _merge(SSS)

//Harmonize to the UKBB effect alleles
forvalues i=1(1)455{
	local ukbb_effect_allele=UKBB_effect_allele[`i']
	local ukbb_other_allele=UKBB_other_allele[`i']
	local hunt_effect_allele=HUNT_effect_allele[`i']
	local hunt_other_allele=HUNT_other_allele[`i']
	
	di "UKBB effect=`ukbb_effect_allele' other=`ukbb_other_allele' HUNT effect=`hunt_effect_allele' other=`hunt_other_allele'
	
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
mregger HUNT_eduyears2_beta HUNT_out_height_beta [aw=1/(HUNT_eduyears2_se^2)] if wood==1, ivw fe
regsave using "results/summary_data",replace detail(all) pval ci

rename HUNT_out_highbloodpressure_beta HUNT_out_highbloodpressure2_beta
rename HUNT_out_highbloodpressure_se  HUNT_out_highbloodpressure2_se

foreach i in HUNT ukbb{
	foreach j in HUNT ukbb{ 
		foreach k in eduyears2 out_highbloodpressure2 out_diabetes{
			mregger `i'_`k'_beta `j'_out_height_beta [aw=1/(`i'_`k'_se^2)] if wood==1, ivw fe
			regsave using "results/summary_data",append detail(all) pval ci
			mregger `i'_`k'_beta `j'_out_bmi_beta [aw=1/(`i'_`k'_se^2)] if locke==1, ivw fe
			regsave using "results/summary_data",append detail(all) pval ci
			}
		}
	}
	
use results/summary_data.dta,clear

duplicates drop

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

cap prog drop save_results
prog def save_results
args index exposure outcome study_outcome study_exposure type outcome_index coef stderr pval row het_chi het_pvalue
replace index=`index' in `row'
replace outcome="`outcome'" in `row'
replace exposure="`exposure'" in `row'
replace study_outcome="`study_outcome'" in `row'
replace study_exposure="`study_exposure'" in `row'
replace type=`type' in `row'
replace outcome_index=`outcome_index' in `row'
replace coef=`coef' in `row'
replace stderr=`stderr' in `row'
replace pval=`pval' in `row'
replace het_chi=`het_chi' in `row'
replace het_pvalue=`het_pvalue' in `row'
end

//Height on education
metan Coef Stderr if inlist(var,"HUNT_eduyears2_beta:HUNT_out_height_beta","ukbb_eduyears2_beta:ukbb_out_height_beta")
save_results 8 out_height eduyears3 same same 5 1 r(ES) r(seES) r(p_z) 1 r(het) r(p_het)
metan Coef Stderr if inlist(var,"HUNT_eduyears2_beta:ukbb_out_height_beta","ukbb_eduyears2_beta:HUNT_out_height_beta")
save_results 8 out_height eduyears3 cross cross 6 1 r(ES) r(seES) r(p_z) 2 r(het) r(p_het)

//BMI on education
metan Coef Stderr if inlist(var,"HUNT_eduyears2_beta:HUNT_out_bmi_beta","ukbb_eduyears2_beta:ukbb_out_bmi_beta")
save_results 5 out_bmi eduyears3 same same 5 1 r(ES) r(seES) r(p_z) 3 r(het) r(p_het)
metan Coef Stderr if inlist(var,"HUNT_eduyears2_beta:ukbb_out_bmi_beta","ukbb_eduyears2_beta:HUNT_out_bmi_beta")
save_results 5 out_bmi eduyears3 cross cross 6 1 r(ES) r(seES) r(p_z) 4  r(het) r(p_het)

//BMI on high blood pressure
metan Coef Stderr if inlist(var,"HUNT_out_highbloodpressure2_beta:HUNT_out_bmi_beta","ukbb_out_highbloodpressure2_beta:ukbb_out_bmi_beta")
save_results 7 out_bmi out_highbloodpressure2 same same 5 5 r(ES) r(seES) r(p_z) 5 r(het) r(p_het)
metan Coef Stderr if inlist(var,"HUNT_out_highbloodpressure2_beta:ukbb_out_bmi_beta","ukbb_out_highbloodpressure2_beta:HUNT_out_bmi_beta")
save_results 7 out_bmi  out_highbloodpressure2 cross cross 6 5 r(ES) r(seES) r(p_z) 6 r(het) r(p_het)

//BMI on diabetes
metan Coef Stderr if inlist(var,"HUNT_out_diabetes_beta:HUNT_out_bmi_beta","ukbb_out_diabetes_beta:ukbb_out_bmi_beta")
save_results 6 out_bmi out_diabetes same same 5 3 r(ES) r(seES) r(p_z) 7 r(het) r(p_het)
metan Coef Stderr if inlist(var,"HUNT_out_diabetes_beta:ukbb_out_bmi_beta","ukbb_out_diabetes_beta:HUNT_out_bmi_beta")
save_results 6 out_bmi out_diabetes cross cross 6 3 r(ES) r(seES) r(p_z) 8  r(het) r(p_het)

keep index-het_pvalue
drop if index==.

save "workingdata/metan",replace

//Clean and output the summary data
//Clean UKBB data to the same format:
use "workingdata/indi_snp_analysis_clean.dta",clear
replace trait="eduyears2" if trait=="eduyears3"

export delimited using "results/snp_summary_data_nmd_190319.csv", replace

//Meta-analyse the IPD results

import delimited results/individual_snps_v7.txt, delimiter(space) clear
drop if trait=="0"
joinby snp using "workingdata/hunt_effect_allele_other",unmatched(master)
drop _m



