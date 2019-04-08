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

fs "rawdata/extracted_snps/wood_height_387*.raw"
foreach i in `r(files)' {
	di "`i'"
	import_snps `i'
	}

use "rawdata/wood_height_387_01.raw.dta",clear
forvalues k=2(1)22{
	if `k'<10{\
		local k="0`k'"
		}
	cap:joinby n_eid using "rawdata/wood_height_387_`k'.raw.dta",unmatched(both)
	cap:drop _m
	}
gen n=_n
compress
save "workingdata/wood_height_387_clean.dta",replace

//387 SNPs
import delimited "$path1/snplists/wood_height_387.txt",clear
rename v1 rsid
save "workingdata/wood_height_387",replace

//Add in effect and other allele
import delimited $path1/rawdata/extracted_snps/wood_height_387.snp-stats2.txt, delimiter(tab) clear varnames(1)

joinby rsid using "workingdata/wood_height_387",
gen n=_n
save "workingdata/wood_height_387.snp-stats",replace


use "workingdata/wood_height_387_clean.dta",clear

joinby n using "workingdata/wood_height_387.snp-stats",unmatched(both) _merge(X)

//Rename each SNP
cap{
local rsid="X"
while "`rsid'"!=""{
	local i=`i'+1
	local rsid=rsid[`i']
	local allelea=allelea[`i']
	local alleleb=alleleb[`i']
	local alleleb_low=lower("`alleleb'")
	di "`allelea' `alleleb' `rsid'"
	rename `rsid'_`alleleb_low' `rsid'_`alleleb'_`allelea'
	}
compress
}

//Keep the common variant for the four triallelic SNPs 
drop rs11177669_C_G  rs12986437_G_A rs1659127_C_G rs16865791_T_C 

drop n-total

save "workingdata/wood_387_clean_final.dta",replace

//Create the height PRS
//Clean coefficients from MR-Base
import delimited "$path1/snplists/IV_height.txt", clear 
gen nx=_n
compress
save "workingdata/wood_387_coefficients.dta",replace

use "workingdata/wood_387_clean_final.dta",clear
gen nx=_n
joinby nx using "workingdata/wood_387_coefficients.dta",unmatched(master)

gen wood_387_prs=0
gen wood_387_wprs=0

forvalues i=1(1)386{

	local rsid=v2[`i']
	local effect_allele=v3[`i']
	local other_allele=v4[`i']
	local coef=v6[`i']
	
	di "`rsid' `effect_allele' `other_allele' `coef'"

	cap:replace wood_387_prs=`rsid'_`effect_allele'_`other_allele'+wood_387_prs if `coef'>0
	cap:replace wood_387_prs=`rsid'_`other_allele'_`effect_allele'*-1+wood_387_prs  if `coef'>0
	
	cap:replace wood_387_prs=`rsid'_`other_allele'_`effect_allele'+wood_387_prs  if `coef'<0
	cap:replace wood_387_prs=`rsid'_`effect_allele'_`other_allele'*-1+wood_387_prs  if `coef'<0
	
	cap:replace wood_387_wprs=`coef'*`rsid'_`effect_allele'_`other_allele'+wood_387_wprs 
	cap:replace wood_387_wprs=`coef'*`rsid'_`other_allele'_`effect_allele'*-1+wood_387_wprs 
	}
compress

//Recode all variants to be height increasing
/*
forvalues i=1(1)386{

	local rsid=v2[`i']
	local effect_allele=v3[`i']
	local other_allele=v4[`i']
	local coef=v6[`i']
	
	//Harmonize to the same effect allele
	cap{
		ds `rsid'_`other_allele'_`effect_allele'
		}
	if _rc==0{	
		replace `rsid'_`other_allele'_`effect_allele'=2-`rsid'_`other_allele'_`effect_allele'
		rename `rsid'_`other_allele'_`effect_allele' `rsid'_`effect_allele'_`other_allele'
		}
	//Flip the negative coefficients
	if `coef'<0 & "`rsid'"!="rs9825951"{
		di "AAA"
		replace `rsid'_`effect_allele'_`other_allele'=2-`rsid'_`effect_allele'_`other_allele'
		rename `rsid'_`effect_allele'_`other_allele' `rsid'_`other_allele'_`effect_allele'
		}
	}
*/	
drop nx-v16

save "workingdata/wood_387_clean_final2.dta",replace
	

