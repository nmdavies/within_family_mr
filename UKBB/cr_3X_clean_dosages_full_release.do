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

fs "rawdata/extracted_snps/all_snps*.raw"
foreach i in `r(files)' {
	di "`i'"
	import_snps `i'
	}

use "rawdata/all_snps_01.raw.dta",clear
forvalues k=2(1)22{
	if `k'<10{\
		local k="0`k'"
		}
	cap:joinby n_eid using "rawdata/all_snps_`k'.raw.dta",unmatched(both)
	cap:drop _m
	}
gen n=_n
compress
save "workingdata/all_snps_clean.dta",replace

//566 SNPs from Wood and Locke
import delimited "snplists/all_snps.txt",clear
rename v1 rsid
save "workingdata/all_snps",replace

//Add in effect and other allele
import delimited "rawdata/extracted_snps/all_snps.snp-stats2.txt", delimiter(tab) clear varnames(1)

joinby rsid using "workingdata/all_snps",
save "workingdata/all_snps.snp-stats",replace

//Identify 390 Wood SNPs (including 4 triallelic SNPs)
import delimited "snplists/wood_eur_snps.txt",clear
rename v1 rsid
joinby rsid using  "workingdata/all_snps.snp-stats", unmatched(both)
keep if _m==3
duplicates drop rsid-total,force
gen n=_n
save "workingdata/wood_eur_snps",replace

//Identify 69 Locke SNPs
import delimited "snplists/locke_eur_snps.txt",clear
rename v1 rsid
joinby rsid using  "workingdata/all_snps.snp-stats", unmatched(both)
keep if _m==3
duplicates drop rsid-total,force
gen n=_n
save "workingdata/locke_eur_snps",replace


use "workingdata/all_snps_clean.dta",clear
drop n
gen n=_n
joinby n using "workingdata/wood_eur_snps",unmatched(both) _merge(X)

//Rename each SNP
//cap{
local rsid="X"
while "`rsid'"!=""{
	local i=`i'+1
	local rsid=rsid[`i']
	local allelea=allelea[`i']
	local alleleb=alleleb[`i']
	local alleleb_low=lower("`alleleb'")
	di "`allelea' `alleleb' `rsid'"
	rename `rsid'_`alleleb_low' X`rsid'_`alleleb'_`allelea'
	order X`rsid'_`alleleb'_`allelea'
	}
	
compress
//}

//Keep the common minor allele for the four triallelic SNPs 
//Tag duplicates
duplicates tag rsid, gen(duplicate)

drop Xrs11177669_C_G  Xrs12986437_G_A Xrs1659127_C_G Xrs16865791_T_C 
drop X__ X

keep X* n_eid
compress
rename X* *
order n_eid
save "workingdata/wood_386_clean_final.dta",replace

//Create the height PRS
//Clean coefficients from MR-Base
import delimited snplists/wood_eur_effect.txt, delimiter(space) bindquote(nobind) varnames(nonames) stripquote(yes) rowrange(2) clear 

keep v7 v3 v10 v2
order v7 v2 v3 v10

rename v3 v3
rename v10 v4 
rename v2 v6
rename v7 v2

gen nx=_n
//compress
save "workingdata/wood_386_coefficients.dta",replace

use "workingdata/wood_386_clean_final.dta",clear
gen nx=_n
joinby nx using "workingdata/wood_386_coefficients.dta",unmatched(master)

gen wood_386_prs=0
gen wood_386_wprs=0

forvalues i=1(1)386{

	local rsid=v2[`i']
	local effect_allele=v3[`i']
	local other_allele=v4[`i']
	local coef=v6[`i']
	
	di "`rsid' `effect_allele' `other_allele' `coef'"

	cap:replace wood_386_prs=`rsid'_`effect_allele'_`other_allele'+wood_386_prs if `coef'>0
	cap:replace wood_386_prs=`rsid'_`other_allele'_`effect_allele'*-1+wood_387_prs  if `coef'>0
	
	cap:replace wood_386_prs=`rsid'_`other_allele'_`effect_allele'+wood_386_prs  if `coef'<0
	cap:replace wood_386_prs=`rsid'_`effect_allele'_`other_allele'*-1+wood_386_prs  if `coef'<0
	
	cap:replace wood_386_wprs=`coef'*`rsid'_`effect_allele'_`other_allele'+wood_386_wprs 
	cap:replace wood_386_wprs=`coef'*`rsid'_`other_allele'_`effect_allele'*-1+wood_386_wprs 
	}
compress

//Recode all variants to be height increasing

drop nx-v4
compress
save "workingdata/wood_386_clean_final2.dta",replace
	

//Create the BMI PRS

use "workingdata/all_snps_clean.dta",clear
drop n
gen n=_n
joinby n using "workingdata/locke_eur_snps",unmatched(both) _merge(X)

//Rename each SNP
local rsid="X"
while "`rsid'"!=""{
	local i=`i'+1
	local rsid=rsid[`i']
	local allelea=allelea[`i']
	local alleleb=alleleb[`i']
	local alleleb_low=lower("`alleleb'")
	di "`allelea' `alleleb' `rsid'"
	rename `rsid'_`alleleb_low' X`rsid'_`alleleb'_`allelea'
	order X`rsid'_`alleleb'_`allelea'
	}
	
compress

keep X* n_eid
drop X_ X
compress
save "workingdata/locke_69_clean_final.dta",replace

//Clean coefficients from MR-Base
import delimited snplists/locke_eur_effect.txt, delimiter(space) bindquote(nobind) varnames(nonames) stripquote(yes) rowrange(1) clear 

keep v3 v4 v5 v7

rename v3 rsid
rename v5 other_allele
rename v4 effect_allele
rename v7 beta

drop in 1

gen nx=_n
destring beta, replace
//compress
save "workingdata/locke_69_coefficients.dta",replace

use "workingdata/locke_69_clean_final.dta",clear
gen nx=_n
joinby nx using "workingdata/locke_69_coefficients.dta",unmatched(master)

gen locke_69_prs=0
gen locke_69_wprs=0

rename X* *

//Recode all variants to be BMI increasing
forvalues i=1(1)69{

	local rsid=rsid[`i']
	local effect_allele=effect_allele[`i']
	local other_allele=other_allele[`i']
	local coef=beta[`i']
	
	di "`rsid' `effect_allele' `other_allele' `coef'"

	cap:replace locke_69_prs=`rsid'_`effect_allele'_`other_allele'+locke_69_prs if `coef'>0
	cap:replace locke_69_prs=`rsid'_`other_allele'_`effect_allele'*-1+locke_69_prs  if `coef'>0
	
	cap:replace locke_69_prs=`rsid'_`other_allele'_`effect_allele'+locke_69_prs  if `coef'<0
	cap:replace locke_69_prs=`rsid'_`effect_allele'_`other_allele'*-1+locke_69_prs  if `coef'<0
	
	cap:replace locke_69_wprs=`coef'*`rsid'_`effect_allele'_`other_allele'+locke_69_wprs 
	cap:replace locke_69_wprs=`coef'*`rsid'_`other_allele'_`effect_allele'*-1+locke_69_wprs 
	}
compress

drop nx-beta
compress
save "workingdata/locke_69_clean_final2.dta",replace
	
