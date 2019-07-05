//Neil Davies 21/11/18
//This merges the phenotypic and genotypic data

use "workingdata/cleaned_phenotype_biobank",clear

//Merge in family indicator
joinby n_eid using "workingdata/family_ids",unmatched(master)
tab _m
//drop if _m!=3
//drop _m

//Merge in genetic data
joinby n_eid using "workingdata/wood_386_clean_final2",unmatched(master) _merge(_merge2)
compress
joinby n_eid using "workingdata/locke_69_clean_final2",unmatched(master) _merge(_merge3)
compress

//Match in the prinipal components
preserve
import delimited "rawdata/data.pca1-40.field_22009.txt", delimiter(space) clear 
drop v2
rename v1 n_eid
forvalues i=3(1)42{
	local j =`i'-2
	rename v`i' pc`j'
	}
save "workingdata/pcs",replace
restore
drop _merge*
joinby n_eid using "workingdata/pcs",unmatched(master)

//Final exclusion list
drop _merge
joinby n_eid using "rawdata/exclusions_181022.dta" , unmatched(master)
drop if _m==3

//Drop people missing height
drop if out_height==.

//Drop people missing education
drop if eduyear==.

//Drop if missing genetics
drop if wood_386_prs==.

//Drop anyone missing education, BMI or height and recalculate number of sibs
drop if eduyears3 ==. |out_bmi ==.|out_height ==.|out_diabetes ==.|out_highbloodpressure2==. 

//Gen indicator for number of siblings
bys famid :gen n_sibs=_N
compress
drop _m

save "workingdata/analysisdata",replace

export delimited using "workingdata/analysisdata.csv", delimiter(tab) replace

drop rs*
export delimited using "workingdata/analysisdata_no_snps.csv", delimiter(tab) replace
