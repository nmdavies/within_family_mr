//Neil Davies 03/07/15
//This creates an indicator for birth month in the UK Biobank data



use "$path2/raw_data/biobank_phenotypes_nmd_150417.dta", clear

//Gen indicator for months and year of birth
gen year_month_birth=100*n_34_0_0+n_52_0_0
tab year_month_birth
replace year_month_birth=100*n_34_0_0+n_52_0_0

//Gen age at assessment centre visit
gen dob=mdy(n_52_0_0,15,n_34_0_0)
gen cov_age= ts_53_0_0-dob 

//Program to combine the fields 0 and 1 for each variable
cap prog drop merge_var
prog def merge_var
replace `1'_0_0=`1'_1_0 if (`1'_0_0==.|`1'_0_0<0 )& (`1'_1_0>0 & `1'_1_0!=.)
end

cap prog drop merge_var2
prog def merge_var2
replace `1'_0_0=`1'_1_`2' if (`1'_0_`2'==.|`1'_0_`2'<0 )& (`1'_1_`2'>0 & `1'_1_`2'!=.)
end

cap prog drop merge_svar
prog def merge_svar
replace `1'_0_0=`1'_1_0 if (`1'_0_0=="")& (`1'_1_0!="")
end

//Clean years of full time education
//Use Okbay method for defining years of education
forval i = 0/1 {
	forval j = 0/5 {
		g EA_`i'_`j' = 17 if n_6138_0_0 == 1
		replace EA_`i'_`j' = 14 if n_6138_`i'_`j' == 2
		replace EA_`i'_`j' = 12 if n_6138_`i'_`j' == 3
		replace EA_`i'_`j' = 12 if n_6138_`i'_`j' == 4
		replace EA_`i'_`j' = 13 if n_6138_`i'_`j' == 5
		replace EA_`i'_`j' = 15 if n_6138_`i'_`j' == 6
		replace EA_`i'_`j' = 11 if n_6138_`i'_`j' == -7
		replace EA_`i'_`j' = . if n_6138_`i'_`j' == -3
		}
	}
// take max 
egen eduyears3 = rmax(EA_*_*)
drop EA_*


//Identify individuals who were not born in England:
merge_var n_1647
gen born_english=(n_1647_0_0==1 & n_20115_0_0==.)


//Generate dob variable
gen yob=n_34_0_0
gen mob=n_52_0_0

//Clean height and BMI data
//Anthropometry
merge_var n_21001
merge_var n_50
gen out_bmi=n_21001_0_0
gen out_height=n_50_0_0

//Gender
gen male=(n_31_0_0==1)

//Income
//See http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=738

merge_var n_738
gen out_income_under_18k=(n_738_0_0>1) if n_738_0_0>0 &n_738_0_0!=.
gen out_income_over_31k=(n_738_0_0>2) if n_738_0_0>0 &n_738_0_0!=.
gen out_income_over_52k=(n_738_0_0>3) if n_738_0_0>0 &n_738_0_0!=.
gen out_income_over_100k=(n_738_0_0>4) if n_738_0_0>0 &n_738_0_0!=.

//Blood pressure
merge_var n_4080
merge_var n_4079
egen out_sys_bp=rowmean(n_4080_0_1 n_4080_0_0)
egen out_dia_bp=rowmean(n_4079_0_1 n_4079_0_0)


//Diagnosed with hypertension
merge_var n_2966
gen out_highbloodpressure=(n_2966_0_0>0 &  n_2966_0_0!=.) if n_2966_0_0>0

//Blood pressure updated for consistency with HUNT:
gen out_highbloodpressure2=(out_highbloodpressure==1|out_sys_bp>140|out_dia_bp>90) if (out_highbloodpressure!=.|out_sys_bp!=.|out_dia_bp!=.)

//Diagnosed with diabetes
merge_var n_2443
merge_var n_2976
gen out_diabetes=(n_2443_0_0==1) if n_2443_0_0>=0 
replace out_diabete=. if n_2976_0_0<=21 & n_2976_0_0!=.

keep out_* male year_month_birth  eduyears3      born_english  yob mob dob  n_eid  cov_age
compress

//Limit the sample just to English born individuals
drop if born_english!=1

joinby n_eid using "working data/covariates",unmatched(master) _merge(XXX)
drop if XXX!=3

//Create month of birth variable

gen mobi=100*yob+mob

compress
drop male XXX


save "$path1/workingdata/cleaned_phenotype_biobank",replace


