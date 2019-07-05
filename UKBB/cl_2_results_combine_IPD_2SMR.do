//Neil Davies 04/07/19
//This merges the IPD and 2SMR versions of the analysis togeather and creates the plots:


//Create the 2SMR results using the summary level SNP results.
//Merge in the meta-analysed UKBB+HUNT results from reg_3_meta_analyse_hunt_ukbb.do and reg_1_height_education.do

use  "workingdata/meta_analysed_ipd",clear

append using "workingdata/metan"

replace type=7 if type==6 
replace type=6 if type==5 & _n>20

replace outcome="eduyears" if outcome=="eduyears3"

sort  outcome exposure type
replace N=61008 if N==.

gen double pval2=pval
replace  pval2=2*normal(-abs(coef/stderr)) 

//Data for text + exact p-values
gen lci=coef -stderr *1.96
gen uci=coef +stderr *1.96

order outcome exposure coef lci uci pval2
compress
save "results/final_table",replace

//Save data for export to R
keep outcome exposure N coef lci uci pval2  
replace coef=coef*10 if exposure=="out_height"
replace lci=lci*10 if exposure=="out_height"
replace uci=uci*10 if exposure=="out_height"

save "results/final_table_r",replace

use "results/final_table",clear
keep outcome exposure N coef stderr lci uci pval2  

//For text
replace coef =coef*100 if outcome =="out_diabetes"|outcome =="out_highbloodpressure"
replace lci =lci*100 if outcome =="out_diabetes"|outcome =="out_highbloodpressure"
replace uci =uci*100 if outcome =="out_diabetes"|outcome =="out_highbloodpressure"

//Rescale per 10cm increase in height
replace coef =coef*10 if exposure =="out_height"
replace lci =lci*10 if exposure =="out_height"
replace uci =uci*10 if exposure =="out_height"


format %9.2f coef lci uci
format %9.1e pval2
