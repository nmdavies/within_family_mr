//Neil Davies 24/07/17
//This imports the individual SNP data into Stata and cleans it into additive format

cap prog drop import_snps
prog def import_snps
args file

import delimited "rawdata/extracted_snps/`file'", delimiter(tab, collapse) encoding(ISO-8859-1) clear
rename fid n_eid
drop iid pat mat phenotype
compress

save "rawdata/`file'.dta",replace 
end

fs "rawdata/extracted_snps/locke_bmi_79_*.raw"
foreach i in `r(files)' {
	di "`i'"
	import_snps `i'
	}


	
use "rawdata/locke_bmi_79_01.raw.dta",clear
forvalues k=2(1)22{
	if `k'<10{\
		local k="0`k'"
		}
	cap:joinby n_eid using "rawdata/locke_bmi_79_`k'.raw.dta",unmatched(both)
	cap:drop _m
	}
gen n=_n
compress
save "workingdata/locke_bmi_79_clean.dta",replace

//79 SNPs
import delimited "$path1/19_within_families_paper/snplists/locke_bmi_79.txt",clear
rename v1 rsid
save "workingdata/locke_bmi_79",replace

//Add in effect and other allele
import delimited rawdata/locke_bmi_79.snp-stats2.txt, delimiter(tab) clear varnames(1)

joinby rsid using "workingdata/locke_bmi_79",
gen n=_n
save "workingdata/locke_bmi_79.snp-stats",replace


use "workingdata/locke_bmi_79_clean.dta",clear

joinby n using "workingdata/locke_bmi_79.snp-stats",unmatched(both) _merge(X)

//Rename each SNP

local rsid="X"
while "`rsid'"!=""{
	cap{
		local i=`i'+1
		local rsid=rsid[`i']
		local allelea=allelea[`i']
		local alleleb=alleleb[`i']
		local alleleb_low=lower("`alleleb'")
		di "`allelea' `alleleb' `rsid'"
		rename `rsid'_`alleleb_low' `rsid'_`alleleb'_`allelea'
		}
	}
compress

drop n-total

save "workingdata/locke_79_clean_final.dta",replace

//Create the height PRS
//Clean coefficients from MR-Base
import delimited $path1/rawdata/locke_bmi_79_coefficients.txt, delimiter(space) clear 
gen nx=_n
compress
save "workingdata/locke_79_coefficients.dta",replace

use "workingdata/locke_79_clean_final.dta",clear
gen nx=_n
joinby nx using "workingdata/locke_79_coefficients.dta",unmatched(master)

gen locke_79_prs=0
gen locke_79_wprs=0

forvalues i=1(1)79{

	local rsid=snp[`i']
	local effect_allele=effect_alleleexposure[`i']
	local other_allele=other_alleleexposure[`i']
	local coef=betaexposure[`i']
	
	di "`rsid' `effect_allele' `other_allele' `coef'"

	cap:replace locke_79_prs=`rsid'_`effect_allele'_`other_allele'+locke_79_prs if `coef'>0
	cap:replace locke_79_prs=`rsid'_`other_allele'_`effect_allele'*-1+locke_79_prs  if `coef'>0
	
	cap:replace locke_79_prs=`rsid'_`other_allele'_`effect_allele'+locke_79_prs  if `coef'<0
	cap:replace locke_79_prs=`rsid'_`effect_allele'_`other_allele'*-1+locke_79_prs  if `coef'<0
	
	cap:replace locke_79_wprs=`coef'*`rsid'_`effect_allele'_`other_allele'+locke_79_wprs 
	cap:replace locke_79_wprs=`coef'*`rsid'_`other_allele'_`effect_allele'*-1+locke_79_wprs 
	}
compress

/*
forvalues i=1(1)79{
	local rsid=snp[`i']
	local effect_allele=effect_alleleexposure[`i']
	local other_allele=other_alleleexposure[`i']
	local coef=betaexposure[`i']

	//Harmonize to the same effect allele
	cap{
		ds `rsid'_`other_allele'_`effect_allele'
		}
	if _rc==0{	
		replace `rsid'_`other_allele'_`effect_allele'=2-`rsid'_`other_allele'_`effect_allele'
		rename `rsid'_`other_allele'_`effect_allele' `rsid'_`effect_allele'_`other_allele'
		}
	//Flip the negative coefficients
	if `coef'<0{
		di "AAA"
		replace `rsid'_`effect_allele'_`other_allele'=2-`rsid'_`effect_allele'_`other_allele'
		rename `rsid'_`effect_allele'_`other_allele' `rsid'_`other_allele'_`effect_allele'
		}
	}
*/
drop nx-v21
save "workingdata/locke_79_clean_final2.dta",replace
	

