//Neil Davies 22/11/18
//This runs the ivreg analysis for the height->education

//Define programs
//Run the four regressions (standard OLS, standard MR with allele score on full sample, standard MR on one of the siblings
// and within family MR.

use "workingdata/analysisdata" , clear

replace famid =_n+22658 if famid ==.

cap prog drop regressions
prog def regressions
args outcome exposure instrument
preserve
if "`outcome'"!="`exposure'"{
	reg `outcome' `exposure' pc1-pc20 cov_age cov_male if n_sibs>1 & n_sibs<8 & n_sibs!=. & (within==1|within==.),ro 
	regsave `exposure' using "results/`outcome'_`exposure'",detail(all) pval replace

	ivreg2 `outcome' (`exposure'=`instrument') pc1-pc20  cov_age cov_male if (within==1|within==.),ro cluster(famid)
	regsave `exposure' using "results/`outcome'_`exposure'",detail(all) pval append
	
	keep if `outcome'!=.
	bys famid:gen n_sibs2=_N
	drop if  n_sibs2==1	
			
	ivreg2 `outcome' (`exposure'=`instrument') pc1-pc20  cov_age cov_male if n_sibs>1 & n_sibs<8 &within==1,ro  cluster(famid)
	regsave `exposure' using "results/`outcome'_`exposure'",detail(all) pval append

	xtivreg `outcome' (`exposure'=`instrument') pc1-pc20  cov_age cov_male if n_sibs>1 & n_sibs<8 ,i(famid) fe  vce(cluster famid )
	regsave `exposure' using "results/`outcome'_`exposure'",detail(all) pval append
	}
else{
}
end

//Basic plot of the results
cap prog drop regressions_plot
prog def  regressions_plot
args outcome exposure instrument	
	
use "results/`outcome'_`exposure'",clear	

gen n=5-_n
sort n

//Plot results using coefplot
//Genetic scores
mkmat coef stderr ,matrix(results) rownames(n)
matrix results_iv=results'

#delimit ;

coefplot (matrix(results_iv) , se(2) ms(C) msize(vsmall) mc(edkblue) ciopts(lc(edkblue) lwidth(vthin))), 
		graphregion(color(white))  plotregion(lc(white)) grid(none) xline(0) ylabel(,labsize(vsmall)) 
		xlabel(,labsize(vsmall)) xtitle("`caption' on `j' (95%CI)",size(vsmall))
		xlabel(,noticks)  xtick(none) ytick(none) title("Outcome=`outcome', exposure=`exposure'")
		coeflabels(1="Ordinary least squares (OLS)"
	2="Mendelian randomization no FE - full sample"
	3="Mendelian randomization 1 sib from siblings sample no FE"
	4="Mendelian randomization siblings FE"
	, wrap(40));	
	#delimit cr	
graph export "results/`outcome'_`exposure'.eps", as(pdf) replace fontface("Calibri")
end

//Open data

use "workingdata/analysisdata", clear

//Set famid = to unique value for the individuals not in families
replace famid =_n+22658 if famid ==.

//Run regressions across all the exposures and outcomes
ds /*eduyears3*/ out_highbloodpressure2 /*out_diabetes*/
foreach i in `r(varlist)'{
	//regressions `i' eduyears3 allele_score_ea
	regressions `i' out_bmi locke_79_wprs
	//regressions `i' out_height wood_387_wprs
	}

//Plot the results using a basic plot
//#delimit ;
//foreach i in 
//out_bmi            out_highbloodpressure out_height     out_diabetes eduyears2{;
//	regressions_plot `i' eduyears2;
//	regressions_plot `i' out_bmi;
//	regressions_plot `i' out_height;
//	};

//Plot the results in a single plot

#delimit ;
use "results/out_bmi_eduyears3.dta",clear;
foreach i in 
eduyears3 out_bmi out_highbloodpressure2 out_height out_diabetes{;
	cap:append using "results/`i'_eduyears3";
	cap:	append using "results/`i'_out_bmi";
	cap:	append using "results/`i'_out_height";
	};
#delimit cr
rename var exposure
rename depvar outcome
order outcome exposure

duplicates drop
drop if exposure==outcome

gen type=1 if cmd=="regress"
replace type=2 if cmd=="ivreg2" & N>40000
replace type=3 if cmd=="ivreg2" & N<40000
replace type=4 if cmd=="xtivreg"


sort exposure outcome type

//Add index for outcome
levels outcome
gen outcome_index=.
foreach i in `r(levels)'{
	local x=`x'+1
	replace outcome_index=`x' if outcome=="`i'"
	}
order outcome exposure type outcome_index
//Input all the OLS estimates into a matrix for one exposure
drop if outcome=="out_sys_bp"| outcome=="out_dia_bp"| outcome=="cov_male"| outcome=="cov_age"

gen cont=0
foreach i in eduyears3 out_sys_bp out_dia_bp out_height out_bmi out_dia_bp cov_num_sisters cov_num_brothers cov_comp_bodysize8 cov_comp_height8 cov_birthweight{
	replace cont=1 if outcome=="`i'"
	}
drop if (outcome=="out_height" & exposure=="out_bmi")|(exposure=="out_height" & outcome=="out_bmi")|(exposure=="out_height" & outcome=="out_highbloodpressure2")|(exposure=="out_height" & outcome=="out_diabetes")

//Generate index per analysis:
gen index=.
replace index=round((_n+1.9)/4)
order index
drop if exposure=="eduyears3"
append using "workingdata/metan"
replace cont=1 if cont==. & outcome=="eduyears3"
replace cont=0 if cont==. 

replace index=6 if index==7 & N==.
replace index=7 if index==8 & N==.
sort index outcome exposure type

/*
mkmat coef stderr if type==1 & cont==1 ,matrix(results) rownames(index)
matrix results_ols=results'
mkmat coef stderr if type==2 & cont==1  ,matrix(results) rownames(index)
matrix results_mr_full_sample=results'
mkmat coef stderr if type==3 & cont==1 ,matrix(results) rownames(index)
matrix results_mr_sibsample=results'
mkmat coef stderr if type==4 & cont==1 ,matrix(results) rownames(index)
matrix results_within_mr=results'
mkmat coef stderr if type==5 & cont==1 ,matrix(results) rownames(index)
matrix results_within_sum_same_mr=results'
mkmat coef stderr if type==6 & cont==1 ,matrix(results) rownames(index)
matrix results_within_cross_mr=results'
*/
gen double pval2=pval
replace  pval2=2*normal(-abs(coef/stderr)) 
/*
#delimit ;

coefplot 	(matrix(results_ols) , se(2) ms(C) msize(vtiny) mc(edkblue) ciopts(lc(edkblue) lwidth(vthin)) offset(0.15))
			(matrix(results_mr_full_sample) , se(2) ms(S) msize(vtiny) mc(eltblue) ciopts(lc(eltblue) lwidth(vthin)) offset(0.05))
			(matrix(results_mr_sibsample) , se(2) ms(T) msize(vtiny) mc(rose) ciopts(lc(rose) lwidth(vthin)) offset(-0.05))
			(matrix(results_within_mr) , se(2) ms(D) msize(vtiny) mc(emerald) ciopts(lc(emerald) lwidth(vthin)) offset(-0.15))
			(matrix(results_within_sum_same_mr) , se(2) ms(D) msize(vtiny) mc(grey) ciopts(lc(grey) lwidth(vthin)) offset(-0.25))
			(matrix(results_within_cross_mr) , se(2) ms(D) msize(vtiny) mc(red) ciopts(lc(red) lwidth(vthin)) offset(-0.35))
		, 
		graphregion(color(white))  plotregion(lc(white)) grid(none) xline(0) ylabel(,labsize(vsmall)) 
		/* yscale(alt axis(2))*/		
		xlabel(,labsize(vsmall)) xtitle("Mean difference in outcome (95%CI)",size(vsmall))
		xlabel(,noticks)  xtick(none) ytick(none)
		coeflabels(5="BMI on education"
	4="Height on education"
	
	, wrap(40))
	/*
	groups(	4 = "N = 18"                                           
            5 = "N = 11"                                           
			, nogap angle(horizontal) labsize(vsmall))*/
	legend(off)
	;	
#delimit cr	
graph save  "results/all_cont_outcome",replace
graph export "results/all_cont_outcome.eps", as(pdf) replace fontface("Calibri")

mkmat coef stderr if type==1 & cont==0 ,matrix(results) rownames(index)
matrix results_ols=results'
mkmat coef stderr if type==2 & cont==0  ,matrix(results) rownames(index)
matrix results_mr_full_sample=results'
mkmat coef stderr if type==3 & cont==0 ,matrix(results) rownames(index)
matrix results_mr_sibsample=results'
mkmat coef stderr if type==4 & cont==0 ,matrix(results) rownames(index)
matrix results_within_mr=results'
mkmat coef stderr if type==5 & cont==0 ,matrix(results) rownames(index)
matrix results_within_sum_same_mr=results'
mkmat coef stderr if type==6 & cont==0 ,matrix(results) rownames(index)
matrix results_within_cross_mr=results'
#delimit ;

coefplot 	(matrix(results_ols) , se(2) ms(C) msize(vtiny) mc(edkblue) ciopts(lc(edkblue) lwidth(vthin)) offset(0.15))
			(matrix(results_mr_full_sample) , se(2) ms(S) msize(vtiny) mc(eltblue) ciopts(lc(eltblue) lwidth(vthin)) offset(0.05))
			(matrix(results_mr_sibsample) , se(2) ms(T) msize(vtiny) mc(rose) ciopts(lc(rose) lwidth(vthin)) offset(-0.05))
			(matrix(results_within_mr) , se(2) ms(D) msize(vtiny) mc(emerald) ciopts(lc(emerald) lwidth(vthin)) offset(-0.15))
			(matrix(results_within_sum_same_mr) , se(2) ms(D) msize(vtiny) mc(grey) ciopts(lc(grey) lwidth(vthin)) offset(-0.25))
			(matrix(results_within_cross_mr) , se(2) ms(D) msize(vtiny) mc(red) ciopts(lc(red) lwidth(vthin)) offset(-0.35))
, 
		graphregion(color(white))  plotregion(lc(white)) grid(none) xline(0) ylabel(,labsize(vsmall)) 
		xlabel(,labsize(vsmall)) xtitle("Risk difference in outcome (95%CI)",size(vsmall))
		xlabel(,noticks)  xtick(none) ytick(none)
		coeflabels(6="BMI on diabetes"
	7="BMI on high blood pressure"
	, wrap(40))
	legend(size(tiny) position(5) ring(0) label(2 "Ordinary least squares UKBB") label(4 "Full sample MR UKBB") label(6 "Sib sample MR UKBB") label(8 "Sib sample sib FE UKBB") label(10 "SMR - UKBB + HUNT same") label(12 "SMR - cross UKBB + HUNT cross"))
	;	
#delimit cr	
graph save  "results/all_bin_outcome",replace
graph export "results/all_bin_outcome.eps", as(pdf) replace fontface("Calibri")
graph combine "results/all_cont_outcome" "results/all_bin_outcome", col(1) imargin(0 0 0 0) ysize(8) xsize(6.26) graphregion(color(white))
graph save  "results/all_bin_combine",replace
graph export "results/all_bin_combine.eps", as(pdf) replace fontface("Calibri")
*/
//Data for text + exact p-values
gen lci=coef -stderr *1.96
gen uci=coef +stderr *1.96

order outcome exposure coef lci uci pval2
compress
save "results/final_table",replace

//Save data for export to R
keep outcome exposure N coef lci uci pval2  

save "results/final_table_r",replace

use "results/final_table",clear
keep outcome exposure N coef stderr lci uci pval2  

metan coef lci uci if outcome=="eduyears3" & exposure=="out_bmi", nooverall 

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
