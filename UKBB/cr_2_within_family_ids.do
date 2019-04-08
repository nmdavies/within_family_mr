//Neil Davies 19/11/18
//This cleans the sibling data from the KINSHIP data provided by UKBB

import delimited $path2/rawdata/ukbb_sibs.csv,clear


local f0 = "red"
local f1 = "green"
twoway (scatter ibs0 kinship if sibpair =="TRUE", msymbol(Oh) mcolor(`f0')) ///
       (scatter ibs0 kinship  if sibpair =="FALSE", msymbol(Oh) mcolor(`f1')) ///
       , legend(off)
*/
//The sibling pairs are identified below
//Generate a long file that list UKBB IDs + famids

keep if sibpair=="TRUE"
gen famid=.

local n=_N
forvalues i=1(1)`n'{
	local id1=id1[`i']
	local id2=id2[`i']
	replace famid=`i' if (id1==`id1'| id1==`id2'| id2==`id1'| id2==`id2') & famid==.
	di "`i'"
	}
//There are a few multi-sibling families. Fix the IDs:

bys id1: egen min=min(famid)
replace famid =min if min<famid

drop min 
bys id2: egen min=min(famid)
replace famid =min if min<famid


preserve 
keep id1 famid 
rename id1 n_eid

save "workingdata/sib1",replace

restore
keep id2 famid 
rename id2 n_eid
append using "workingdata/sib1"

//Drop duplicate observations
duplicates drop

//This gives a file that is one row per ID with a single family ID across all siblings:
//Generate within family ID, to restrict to 2 siblings per analysis restrict on within_fam_id<3
bys famid: gen within_fam_id=_n
compress
save "workingdata/family_ids",replace
