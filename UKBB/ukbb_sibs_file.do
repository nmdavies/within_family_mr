//Neil Davies 20/05/19
//This creates a ukbb_sibs.csv file with the IEU ids on it.


import delimited "rawdata/ukbb_sibs.csv", encoding(ISO-8859-1)
save "workingdata/ukbb_sibs.dta"

import delimited "rawdata/linker.fam", encoding(ISO-8859-1)clear
rename v1 ieu_id
rename v2 n_eid
save "workingdata/linker.fam.dta",replace


use  "workingdata/ukbb_sibs.dta",clear
rename id1 n_eid
joinby n_eid using "workingdata/linker.fam.dta",unmatched(both)
gen ieu_id1=ieu_id
drop ieu_id 
drop if _m==2
drop _m
rename n_eid id1
rename id2 n_eid
joinby n_eid using "workingdata/linker.fam.dta",unmatched(both)
gen ieu_id2=ieu_id

drop if _m==2
drop _merge 
drop ieu_id 
drop id1 n_eid 
drop v1
order ieu_id1 ieu_id2 
save "workingdata/ukbb_sibs_ieu_ids.dta"
