//Neil Davies 19/03/19
//This estiamtes the correlation between the within family and overall coefficients 


use "workingdata/indi_snp_analysis_clean",clear

keep if sample_size>20000

rename beta beta_within
rename se se_within
keep SNP trait beta_within se_within
save "workingdata/indi_snp_analysis_clean_within",replace

use "workingdata/indi_snp_analysis_clean",clear
keep if sample_size<20000

joinby SNP trait using "workingdata/indi_snp_analysis_clean_within.dta", _merge(XX) 


//Select the height SNPs
joinby SNP using "workingdata/locke_79_rsid.dta",unmatched(master)
gen locke_SNP=(_merge==3)

//386 Wood SNPs
joinby SNP using "workingdata/wood_387_rsid.dta", unmatched(master) _merge(XX2)
gen wood_SNP=(XX2==3)

reg  beta_within beta if trait=="out_bmi" & locke_SNP==1
test _b[beta]==1
di r(p)
reg  beta_within beta if trait=="out_height" & wood_SNP==1
test _b[beta]==1
di r(p)
reg beta_within beta if trait=="eduyears3"  & wood_SNP==1
test _b[beta]==1
di r(p)
reg beta_within beta if trait=="eduyears3"  & locke_SNP==1
test _b[beta]==1
di r(p)
reg beta_within beta if trait=="out_diabetes"& locke_SNP==1
test _b[beta]==1
di r(p)
reg beta_within beta if trait=="out_highbloodpressure2"& locke_SNP==1
test _b[beta]==1
di r(p)

tabstat  beta_within beta if trait=="out_bmi" & locke_SNP==1
tabstat  beta_within beta if trait=="out_height" & wood_SNP==1
