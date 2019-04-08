//Neil Davies 06/12/18
//This creates the Okbay score using 75 SNPs chosen using MR-base

cd "$path1/19_within_families_paper"

use "$path3/workingdata/okbay_hill_snps_clean.dta",clear

gen allele_score_ea=0

joinby n using "workingdata/okbay_coef",unmatched(master)

forvalues i=1(1)75{
	local rsid=okbay_rsid[`i']
	local effect=okbay_effect[`i']
	local other=okbay_other[`i']
	local beta=okbay_beta[`i']
	cap:ds `rsid'_`effect'_`other'
	if _rc==111{
		cap:ds `rsid'_`other'_`effect'
		if _rc!=111{
			di "Not reversed `rsid' effect=`effect' other=`other' beta=`beta'"
			replace allele_score_ea=`rsid'_`other'_`effect'*`beta'+allele_score_ea
			}
		else{
			di "Error `rsid' effect=`effect' other=`other' beta=`beta'"
			}
		}
	else{
		di "Reversed `rsid' effect=`effect' other=`other' beta=`beta'"
		replace allele_score_ea=`rsid'_`effect'_`other'*`beta'*-1+allele_score_ea
		}
	}
	
keep n_eid allele_score_ea
save "workingdata/allele_score_ea",replace
