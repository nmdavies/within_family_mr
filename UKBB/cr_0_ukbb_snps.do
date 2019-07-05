//Neil Davies 25/04/19
//This creates a list of all the SNPs in UKBB


forvalues i=1(1)22{
	if `i'<10{
		local j="0"+string(`i')
		}
	else{
		local j="`i'"
		}
	import delimited rawdata/data.chr`j'.snp-stats, delimiter(tab) varnames(10) colrange(2:2) clear
	compress
	save workingdata/snps_`j',replace
	}
	
use workingdata/snps_01,clear
forvalues i=2(1)22{
	if `i'<10{
		local j="0"+string(`i')
		}
	else{
		local j="`i'"
		}
	append using workingdata/snps_`j'
	}
save workingdata/snps,replace
