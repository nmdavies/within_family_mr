//Neil Davies 06/12/18
//This runs the individual SNP analysis for the within families paper

//Load data for siblings
use "workingdata/analysisdata" if n_sibs !=1 & n_sibs <8, clear

//Create results file
xtreg out_height rs1000940_G_A pc1-pc20 cov_male cov_age ,ro i(famid ) fe  cluster(famid)
regsave rs1000940_G_A using "results/indi_snp_analysis",replace detail(all) pval ci
	
//Generate indicators for samples for each outcome
foreach outcome in out_bmi out_height out_diabetes out_highbloodpressure2 eduyears3{
	preserve
	
	//Drop families with only one person after accounting for missing outcomes
	keep if `outcome'!=.
	bys famid:gen n_sibs2=_N
	drop if  n_sibs2==1

	ds rs*
	foreach j in `r(varlist)'{
		xtreg `outcome' `j' pc1-pc20 cov_male cov_age ,ro i(famid ) fe cluster(famid)
		regsave `j' using "results/indi_snp_analysis",append detail(all) pval ci
		reg `outcome' `j' pc1-pc20 cov_male cov_age ,ro  cluster(famid)
		regsave `j' using "results/indi_snp_analysis",append detail(all) pval ci
		}
	restore
	}
	

//Get the EAFs
gen EAF=.
gen SNP=""
ds rs*
foreach j in `r(varlist)'{
	tabstat `j',save
	local k=`k'+1
	replace EAF=el(r(StatTotal),1,1)/2 in `k'
	replace SNP=word(subinstr("`j'","_"," ",2),1) in `k'
	}
keep SNP EAF
compress
drop if SNP==""
save "workingdata/EAF",replace

use "results/indi_snp_analysis",clear
duplicates drop
gen method="Within FE" if cmd=="xtreg"
replace method="Standard MR" if cmd=="regress"

rename depvar trait
gen SNP=word(subinstr(var,"_"," ",2),1)
gen effect_allele=word(subinstr(var,"_"," ",2),2)
gen other_allele=word(subinstr(var,"_"," ",2),3)

rename coef beta
rename std se
rename pval pvalue 

rename N sample_size

order trait SNP effect_allele other_allele beta se pvalue sample_size method
keep trait SNP effect_allele other_allele beta se pvalue sample_size method
compress

joinby SNP using "workingdata/EAF",unmatched(master)
drop _m
save "workingdata/indi_snp_analysis_clean",replace
