*Code for: Systematic Review and Meta-Analysis of the Associations Between Body Mass Index, Prostate Cancer, Advanced Prostate Cancer and Prostate Specific Antigen
*Code written by Sean Harrison between 11/2016 and 07/2019

*All file paths removed
*IPD analysis restricted to imputation and beyond (i.e. no data cleaning)

*IPD data is restricted, so not put on Github

*File list therefore includes:
	*Files from systematic review:
		*PCa categories
		*PCa MD
		*PSA categories
		*PSA MD
	*Files from AD analysis:
		*IPD PCa - complete
		*IPD PCa - complete (categorical)
		*IPD PSA - complete
		*IPD PSA - complete (categorical)
*These files should be enough to recreate all data processing and analyses (path names will need editing)

*Sean Harrison - 2019/07/08

********************************************************************************

global cd_folder ""
global cd_data "$cd_folder\Data"
global cd_graphs "$cd_folder\Graphs"
global cd_tables "$cd_folder\Tables"

cd "$cd_data"

*ssc install bspline
*ssc install metan
*ssc install metareg

*Merging and analysis of IPD
{
*Replace biopsy to 0 if the PSA is below the study threshold
use "All.dta", clear

replace biopsy = 1 if study == "ProtecT" & psa >= 3
replace biopsy = 0 if study == "ProtecT" & psa < 3

replace pca = 1 if gleason >2 & gleason <.
replace biopsy = 1 if pca == 1

replace pca = . if pca == 0 & biopsy == 0
replace pca = 0 if biopsy == 1 & pca != 1

rename study study2
gen study = 2 if study2 == "Krimpen"
replace study = 3 if study2 == "PCPT"
replace study = 4 if study2 == "PLCO"
replace study = 5 if study2 == "ProtecT"
label define study 2 "Krimpen" 3 "PCPT" 4 "PLCO" 5 "ProtecT"
label values study study
drop study2

replace psa = 0.02 if psa == 0
gen logpsa = ln(psa)
replace bmi = . if bmi >100 | bmi <12

*Generate T-score binary variable
gen pcat_binary = 1 if pcat > 6 & pcat < .
replace pcat_binary = 0 if pcat <=6

*Create advanced PCa measure (1 = advanced)
*Advanced = Gleason 7+, T3+, N= or M=1
gen advanced_pca = 1 if pca == 1 & ((pcat>=7 & pcat !=. )| (pcan == 1 & pcan !=.) | (pcam == 1 & pcam !=.))
replace advanced_pca = 0 if pca == 1 & advanced_pca == .
replace advanced_pca = 0 if pca == 0 

*Ethnicity - Krimpen are almost all white
replace ethnicity = 1 if study == 2
keep if ethnicity == 1

*Screening arm of PLCO only
drop if arm == 2

*Keep only control arm of PCPT
drop if treatment == 1

save "IPD All clear.dta", replace

**********

*IPD Meta-analysis without imputation (complete case analysis)

use "IPD All clear.dta", clear
set more off
drop if bmi == .
drop if psa == .
replace ethnicity = 1 if study == 2
drop if ethnicity != 1
replace pca = 0 if pca == . //Assume that all non-biopsied are also controls
replace advanced_pca = 0 if advanced_pca == . //Assume that all non-biopsied are also controls for advanced PCa too
gen psa5 = .
gen psa5_se = .
gen psa5_nopca = .
gen psa5_se_nopca = .
gen pca5 = .
gen pca5_se = .
gen advanced5 = .
gen advanced5_se = .
gen pcapsa = .
gen pcapsa_se = .
gen apcapsa = .
gen apcapsa_se = .
gen study2 = "Krimpen" in 2
replace study2 = "PCPT" in 3
replace study2 = "PLCO" in 4
replace study2 = "ProtecT" in 5

foreach i in 2 3 4 5 {
	regress logpsa bmi age pca if study == `i'
	replace psa5 = _b[bmi]*5 in `i'
	replace psa5_se = _se[bmi]*5 in `i'
	regress logpsa bmi age if study == `i'
	replace psa5_nopca = _b[bmi]*5 in `i'
	replace psa5_se_nopca = _se[bmi]*5 in `i'	
	logit pca bmi age familyhistory if study == `i'
	replace pca5 = _b[bmi]*5 in `i'
	replace pca5_se = _se[bmi]*5 in `i'
	logit advanced_pca bmi age familyhistory if study == `i'
	replace advanced5 = _b[bmi]*5 in `i'
	replace advanced5_se = _se[bmi]*5 in `i'
	regress psa age bmi pca advanced_pca if study == `i'
	replace pcapsa = _b[pca] in `i'
	replace pcapsa_se = _se[pca] in `i'
	replace apcapsa = _b[advanced_pca] in `i'
	replace apcapsa_se = _se[advanced_pca] in `i'
	
}

keep psa5-advanced5_se study2
rename study2 study
drop if psa5 == .

gen type = "Not Imputed"

preserve
keep pca5 pca5_se study type
rename pca5 md5
rename pca5_se md5_se
gen outcome = "All"
save "Complete case - pca.dta", replace
restore

keep advanced5 advanced5_se study type
rename advanced5 md5
rename advanced5_se md5_se
gen outcome = "Advanced"
save "Complete case - advanced pca.dta", replace

********************************************************************************

*IPD imputation and analyses

*Imputations
{

*Imputations - linear
*All
use "IPD All clear.dta", clear	

mi set flong
replace _mi_miss = 1 if pca == . | familyhistory == . | age == . | bmi == . | logpsa == .
mi register imputed pca familyhistory age bmi logpsa 
mi impute chained (logit) pca (logit) familyhistory (reg) age (reg) bmi (reg) logpsa = i.study, add(100) rseed(100) dots
save "Imputed data (linear).dta", replace

*Advanced
use "IPD All clear.dta", clear	
mi set flong
replace _mi_miss = 1 if advanced_pca == . | familyhistory == . | age == . | bmi == . | logpsa == .
mi register imputed familyhistory advanced_pca age bmi logpsa 
mi impute chained (logit) advanced_pca (logit) familyhistory (reg) age (reg) bmi (reg) logpsa = i.study, add(100) rseed(100) dots	
save "Imputed data (linear, advanced).dta", replace

*Imputations - categorical
*All
use "IPD All clear.dta", clear
gen bmi_cat = 0 if bmi <25
replace bmi_cat = 1 if bmi >= 25 & bmi<30
replace bmi_cat = 2 if bmi	>= 30 & bmi <.
mi set flong
replace _mi_miss = 1 if pca == . | familyhistory == . | age == . | bmi == . | logpsa == .
mi register imputed pca familyhistory age bmi_cat logpsa 
mi impute chained (logit) pca (logit) familyhistory (reg) age (ologit) bmi_cat (reg) logpsa = i.study, add(100) rseed(100) dots
save "Imputed data (categorical).dta", replace

*Advanced
use "IPD All clear.dta", clear	
gen bmi_cat = 0 if bmi <25
replace bmi_cat = 1 if bmi >= 25 & bmi<30
replace bmi_cat = 2 if bmi	>= 30 & bmi <.
mi set flong
replace _mi_miss = 1 if advanced_pca == . | familyhistory == . | age == . | bmi == . | logpsa == .
mi register imputed familyhistory advanced_pca age bmi_cat logpsa 
mi impute chained (logit) advanced_pca (logit) familyhistory (reg) age (ologit) bmi_cat (reg) logpsa = i.study, add(100) rseed(100) dots	
save "Imputed data (categorical, advanced).dta", replace
}
********************************************************************************
*Analyses
{
*LINEAR Analysis
use "Imputed data (linear).dta", clear

*Individual results - All (linear)

foreach k in 2 4 5 {
	dis "Start of analysis `k'"
	preserve
	keep if study == `k'

	mi estimate, dots: regress logpsa bmi pca age if study == `k'
	matrix b = e(b_mi)
	matrix V = e(V_mi)
	local psa5_`k' = b[1,1]*5 
	local psa5_se_`k' = sqrt(V[1,1])*5 
	local psa5_p_`k' = 2*normal(-abs(b[1,1]/sqrt(V[1,1]))) 
	
	mi estimate, dots: logit pca bmi age familyhistory if study == `k'
	matrix b = e(b_mi)
	matrix V = e(V_mi)
	local pca5_`k' = b[1,1]*5
	local pca5_se_`k' = sqrt(V[1,1])*5 
	local pca5_p_`k' = 2*normal(-abs(b[1,1]/sqrt(V[1,1])))
	
	qui sum pca if study == `k' & _mi_m != 0
	local ncases_`k' = r(N)*r(mean)/100
	local ncontrols_`k' = (r(N)-r(N)*r(mean))/100
	local ntotal_`k' = r(N)/100
	
	qui sum age if study == `k' & _mi_m != 0
	local age_`k' = r(mean)
	
	qui sum bmi if study == `k' & _mi_m != 0
	local bmi_all_`k' = r(mean)
	local bmi_all_sd_`k' = r(sd)
	
	qui sum bmi if study == `k' & pca == 1 & _mi_m != 0
	local bmi_cases_`k' = r(mean)
	local bmi_cases_sd_`k' = r(sd)
	
	qui sum bmi if study == `k' & pca == 0 & _mi_m != 0
	local bmi_controls_`k' = r(mean)
	local bmi_controls_sd_`k' = r(sd)
	
	restore
	dis "End of analysis `k'"
}

keep in 1/5

foreach var in psa5 pca5 {
	gen `var' = .
	gen `var'_se = .
	gen `var'_p = .
}

gen ncases = .
gen ncontrols = .
gen ntotal = .
gen bmi_cases = .
gen bmi_cases_sd = .
gen bmi_controls = .
gen bmi_controls_sd = .
gen bmi_all = .
gen bmi_all_sd = .

gen studyname = ""
replace studyname = "Krimpen" in 2
replace studyname = "PCPT" in 3
replace studyname = "PLCO" in 4
replace studyname = "ProtecT" in 5

gen studylocation = ""
replace studylocation = "Netherlands" in 2
replace studylocation = "USA" in 4
replace studylocation = "UK" in 5

gen type = "IPD"
drop ethnicity
gen ethnicity = "White" in 1/5
gen overallriskofbias = "Medium"

foreach k in 2 4 5 {
	foreach var in psa5 pca5 {
		qui replace `var' = ``var'_`k'' in `k'
		qui replace `var'_se = ``var'_se_`k'' in `k'
		qui replace `var'_p = ``var'_p_`k'' in `k'
	}
}

foreach k in 2 4 5 {
	foreach var of varlist ncases-bmi_all_sd {
		qui replace `var' = ``var'_`k'' in `k'
	}
}

drop in 3
drop in 1

keep psa5* pca5* ncases ncontrols ntotal bmi_cases* bmi_con* bmi_all* studyname type studylocation ethnicity overall

preserve
drop pca5* 
rename psa5 beta5
rename psa5_se beta5_se
rename psa5_p pvalue
gen whichvariablesareadjustedfor = "Age, ethnicity, prostate cancer"
save "IPD PSA - complete", replace
restore

drop psa5*
rename pca5 md5
rename pca5_se md5_se
rename pca5_p pvalue
gen outcome = "All"
gen whichvariablesareadjustedfor = "Age, ethnicity"
save "IPD PCa - complete", replace

*Individual results - advanced (linear)

use "Imputed data (linear, advanced).dta", clear

foreach k in 2 4 5 {
	dis "Start of analysis `k'"
	preserve
	keep if study == `k'
	
	mi estimate, dots: logit advanced_pca bmi age familyhistory if study == `k'
	matrix b = e(b_mi)
	matrix V = e(V_mi)
	local apca5_`k' = b[1,1]*5
	local apca5_se_`k' = sqrt(V[1,1])*5 
	local apca5_p_`k' = 2*normal(-abs(b[1,1]/sqrt(V[1,1])))
	
	qui sum advanced_pca if study == `k' & _mi_m != 0
	local ncases_`k' = r(N)*r(mean)/100
	local ncontrols_`k' = (r(N)-r(N)*r(mean))/100
	local ntotal_`k' = r(N)/100
	
	qui sum age if study == `k' & _mi_m != 0
	local age_`k' = r(mean)
	
	qui sum bmi if study == `k' & _mi_m != 0
	local bmi_all_`k' = r(mean)
	local bmi_all_sd_`k' = r(sd)
	
	qui sum bmi if study == `k' & pca == 1 & _mi_m != 0
	local bmi_cases_`k' = r(mean)
	local bmi_cases_sd_`k' = r(sd)
	
	qui sum bmi if study == `k' & pca == 0 & _mi_m != 0
	local bmi_controls_`k' = r(mean)
	local bmi_controls_sd_`k' = r(sd)
	
	restore
	dis "End of analysis `k'"
}

keep in 1/5

foreach var in psa5 apca5 {
	gen `var' = .
	gen `var'_se = .
	gen `var'_p = .
}

gen ncases = .
gen ncontrols = .
gen ntotal = .
gen bmi_cases = .
gen bmi_cases_sd = .
gen bmi_controls = .
gen bmi_controls_sd = .
gen bmi_all = .
gen bmi_all_sd = .

gen studyname = ""
replace studyname = "Krimpen" in 2
replace studyname = "PCPT" in 3
replace studyname = "PLCO" in 4
replace studyname = "ProtecT" in 5

gen studylocation = ""
replace studylocation = "Netherlands" in 2
replace studylocation = "USA" in 4
replace studylocation = "UK" in 5

gen type = "IPD"
drop ethnicity
gen ethnicity = "White" in 1/5
gen overallriskofbias = "Medium"

foreach k in 2 4 5 {
	foreach var in apca5 {
		qui replace `var' = ``var'_`k'' in `k'
		qui replace `var'_se = ``var'_se_`k'' in `k'
		qui replace `var'_p = ``var'_p_`k'' in `k'
	}
}

foreach k in 2 4 5 {
	foreach var of varlist ncases-bmi_all_sd {
		qui replace `var' = ``var'_`k'' in `k'
	}
}

drop in 3
drop in 1

keep apca5* ncases ncontrols ntotal bmi_cases* bmi_con* bmi_all* studyname type studylocation ethnicity overall
rename apca5 md5
rename apca5_se md5_se
rename apca5_p pvalue
gen outcome = "Advanced"
gen whichvariablesareadjustedfor = "Age, ethnicity, log(PSA)"
append using "IPD PCa - complete.dta"
save "IPD PCa - complete", replace

********************************************************************************

*CATEGORICAL analysis

set more off
use "Imputed data (categorical).dta", clear

foreach var in obese overweight {
	foreach x in pca apca psa {
		gen `var'_`x' = .
		gen `var'_`x'_se = .
		gen `var'_`x'_p = .
	}
}

*Ns and mean BMIs
foreach var in normal overweight obese {
	gen `var'_ncases = .
	gen `var'_ncontrols = .
	gen `var'_nacases = .
	gen `var'_nacontrols = .	
	gen `var'_ntotal = .
}
gen or_hr = "OR"

gen normal_bmi = .
gen overweight_bmi = .
gen obese_bmi = .
gen study_bmi_mean = .
gen study_bmi_sd = .
gen normal_logpsa = .
gen overweight_logpsa = .
gen obese_logpsa = .
gen normal_logpsa_sd = .
gen overweight_logpsa_sd = .
gen obese_logpsa_sd = .

gen ncases_total=.
gen nacases_total=.
gen ncontrols_total = .
gen nacontrols_total = .
gen ntotal = .

gen normal = 1 if bmi_cat == 0
replace normal = 0 if bmi_cat == 1 | bmi_cat == 2
gen overweight = 1 if bmi_cat == 1
replace overweight = 0 if bmi_cat == 0 | bmi_cat == 2
gen obese = 1 if bmi_cat == 2
replace obese = 0 if bmi_cat == 0 | bmi_cat == 1

save "Imputed data (categorical) - working.dta", replace

use "Imputed data (categorical) - working.dta", clear

*Individual results
foreach k in 2 4 5 {
	dis "Start of analysis `k'"
	preserve
	keep if study == `k'
	
	*Categorical
	
	mi estimate, dots: regress logpsa overweight obese age pca if study == `k'
	matrix b = e(b_mi)
	matrix V = e(V_mi)
	local overweight_psa_`k' = b[1,1]
	local overweight_psa_se_`k' = sqrt(V[1,1]) 
	local overweight_psa_p_`k' = 2*normal(-abs(b[1,1]/sqrt(V[1,1]))) 
	local obese_psa_`k' = b[1,2]
	local obese_psa_se_`k' = sqrt(V[2,2]) 
	local obese_psa_p_`k' = 2*normal(-abs(b[1,2]/sqrt(V[2,2])))
	
	mi estimate, dots: logit pca overweight obese age familyhistory if study == `k'
	matrix b = e(b_mi)
	matrix V = e(V_mi)
	local overweight_pca_`k' = b[1,1]
	local overweight_pca_se_`k' = sqrt(V[1,1]) 
	local overweight_pca_p_`k' = 2*normal(-abs(b[1,1]/sqrt(V[1,1]))) 
	local obese_pca_`k' = b[1,2]
	local obese_pca_se_`k' = sqrt(V[2,2]) 
	local obese_pca_p_`k' = 2*normal(-abs(b[1,2]/sqrt(V[2,2])))

	restore
	
	dis "End of analysis `k'"
}

*Ns for normal, overweight, obese
foreach i in 2 4 5 {	
	foreach weight in normal overweight obese {
				
		qui sum pca if `weight' == 1 & study == `i' & _mi_m != 0
		local `weight'_ncases_`i' = r(mean)*r(N)/100
		local `weight'_ncontrols_`i' = r(N)/100-r(mean)*r(N)/100
		local `weight'_ntotal_`i' = r(N)/100 
		
		qui sum logpsa if `weight' == 1 & study == `i' & _mi_m != 0
		local `weight'_logpsa_`i' = r(mean)
		local `weight'_logpsa_sd_`i' = r(sd)
	}
}

keep in 1/5

gen studyname = ""
replace studyname = "Krimpen" in 2
replace studyname = "PCPT" in 3
replace studyname = "PLCO" in 4
replace studyname = "ProtecT" in 5

gen studylocation = ""
replace studylocation = "Netherlands" in 2
replace studylocation = "USA" in 4
replace studylocation = "UK" in 5

gen type = "IPD"
drop ethnicity
gen ethnicity = "White" in 1/5
gen overallriskofbias = "Medium"

foreach k in 2 4 5 {
	foreach var in overweight_psa overweight_pca obese_psa obese_pca {
		qui replace `var' = ``var'_`k'' in `k'
		qui replace `var'_se = ``var'_se_`k'' in `k'
		qui replace `var'_p = ``var'_p_`k'' in `k'
	}
}

foreach k in 2 4 5 {
	foreach var in normal overweight obese {
		foreach var2 in ntotal ncases ncontrols logpsa logpsa_sd {
			qui replace `var'_`var2' = ``var'_`var2'_`k'' in `k'
		}
	}
}

keep studyname-overallrisk normal* overweight* obese* study_bmi*
drop in 3
drop in 1

drop normal overweight obese

save "Imputed data (categorical) - working2.dta", replace
use "Imputed data (categorical) - working2.dta", clear

drop *apca*
drop *nacases
drop *nacontrols
drop study_bmi_mean study_bmi_sd
rename obese_pca obese
rename obese_pca_se obese_se
rename obese_pca_p obese_p
rename overweight_pca overweight
rename overweight_pca_se overweight_se
rename overweight_pca_p overweight_p

foreach var in normal overweight obese {
	rename `var'_ncases `var'_cases
	rename `var'_ncontrols `var'_controls
}

gen outcome = "All"
save "Imputed data (categorical) - working3.dta", replace

********************************************************************************

*Advanced
set more off
use "Imputed data (categorical, advanced).dta", clear

foreach var in obese overweight overobese {
	foreach x in apca {
		gen `var'_`x' = .
		gen `var'_`x'_se = .
		gen `var'_`x'_p = .
	}
}

*Ns and mean BMIs
foreach var in normal overweight obese {
	gen `var'_nacases = .
	gen `var'_nacontrols = .
	gen `var'_natotal = .
}
gen or_hr = "OR"

gen nacases_total=.
gen nacontrols_total = .

gen normal = 1 if bmi_cat == 0
replace normal = 0 if bmi_cat == 1 | bmi_cat == 2
gen overweight = 1 if bmi_cat == 1
replace overweight = 0 if bmi_cat == 0 | bmi_cat == 2
gen obese = 1 if bmi_cat == 2
replace obese = 0 if bmi_cat == 0 | bmi_cat == 1

save "Imputed data (categorical, advanced) - working.dta", replace

use "Imputed data (categorical, advanced) - working.dta", clear

*Individual results
foreach k in 2 4 5 {
	dis "Start of analysis `k'"
	preserve
	keep if study == `k'
	
	*Categorical
	
	mi estimate, dots: logit advanced_pca overweight obese age familyhistory if study == `k'
	matrix b = e(b_mi)
	matrix V = e(V_mi)
	local overweight_apca_`k' = b[1,1]
	local overweight_apca_se_`k' = sqrt(V[1,1]) 
	local overweight_apca_p_`k' = 2*normal(-abs(b[1,1]/sqrt(V[1,1]))) 
	local obese_apca_`k' = b[1,2]
	local obese_apca_se_`k' = sqrt(V[2,2]) 
	local obese_apca_p_`k' = 2*normal(-abs(b[1,2]/sqrt(V[2,2])))

	restore
	
	dis "End of analysis `k'"
}

*Ns for normal, overweight, obese
foreach i in 2 4 5 {	
	
	qui sum advanced_pca if normal == 1 & study == `i' & _mi_m != 0
	local normal_nacases_`i' = r(mean)*r(N)/100
	local normal_nacontrols_`i' = r(N)/100-r(mean)*r(N)/100
	local normal_natotal = r(N)/100
	
	qui sum advanced_pca if overweight == 1 & study == `i' & _mi_m != 0
	local overweight_nacases_`i' = r(mean)*r(N)/100
	local overweight_nacontrols_`i' = r(N)/100-r(mean)*r(N)/100
	local overweight_natotal = r(N)/100
	
	qui sum advanced_pca if obese == 1 & study == `i' & _mi_m != 0
	local obese_nacases_`i' = r(mean)*r(N)/100 
	local obese_nacontrols_`i' = r(N)/100-r(mean)*r(N)/100
	local obese_natotal = r(N)/100
	
	qui sum advanced_pca if study == `i' & _mi_m != 0
	local nacases_total = r(N)*r(mean)/100
	local nacontrols_total = (r(N)-r(N)*r(mean))/100

}

keep in 1/5

foreach k in 2 4 5 {
	foreach var in overweight_apca obese_apca {
		qui replace `var' = ``var'_`k'' in `k'
		qui replace `var'_se = ``var'_se_`k'' in `k'
		qui replace `var'_p = ``var'_p_`k'' in `k'
	}
}

foreach k in 2 4 5 {
	foreach var in normal overweight obese {
		foreach var2 in nacases nacontrols {
			qui replace `var'_`var2' = ``var'_`var2'_`k'' in `k'
		}
	}
}

drop in 3
drop in 1
keep *nacases *nacontrols *natotal overweight_apca* obese_apca*

save "Imputed data (categorical, advanced) - working2.dta", replace
use "Imputed data (categorical, advanced) - working2.dta", clear

rename obese_apca obese
rename obese_apca_se obese_se
rename obese_apca_p obese_p
rename overweight_apca overweight
rename overweight_apca_se overweight_se
rename overweight_apca_p overweight_p

foreach var in normal overweight obese {
	rename `var'_nacases `var'_cases
	rename `var'_nacontrols `var'_controls
	rename `var'_natotal `var'_total
}

gen outcome = "Advanced"
gen studyname = ""
replace studyname = "Krimpen" in 1
replace studyname = "PLCO" in 2
replace studyname = "ProtecT" in 3

gen studylocation = ""
replace studylocation = "Netherlands" in 1
replace studylocation = "USA" in 2
replace studylocation = "UK" in 3

gen type = "IPD"
gen ethnicity = "White"
gen overallriskofbias = "Medium"

append using "Imputed data (categorical) - working3.dta"

preserve
drop in 1/3
drop obese-obese_controls
save "IPD PSA - complete (categorical).dta", replace
restore

drop obese_psa-overweight_psa_p
gen whichvariablesareadjustedfor = "Age, ethnicity, log(PSA)"
save "IPD PCa - complete (categorical).dta", replace



}
********************************************************************************

*Graphs & Tables for appendix
{

*Cubic splines (complete case)
*BMI - PSA
use "IPD All clear.dta", clear

sort study bmi
gen logpsa_age = .
gen study_i = .
gen se = .
gen beta = .
foreach i in 2 3 4 5 {
	reg logpsa age if study == `i'
	local beta = _b[age]
	replace logpsa_age = logpsa+logpsa*(age-60)*`beta' if study == `i'
	bspline if study == `i', xvar(bmi) knots(0 24.8 27 30 100) p(3) gen(_study_`i')
	qui reg logpsa_age _study_`i'*, noconstant
	predict study_`i'
	drop _study_`i'*
	qui reg logpsa_age bmi if study == `i'
	qui replace study_i = `i' in `i'
	qui replace beta = _b[bmi] in `i'
	qui replace se = _se[bmi] in `i'	
}

twoway (line study_2 bmi if bmi >= 19.8 & bmi <=34.5) (line study_3 bmi if bmi >= 20.4 & bmi <=41) (line study_4 bmi if bmi >= 19.9 & bmi <=40.5) (line study_5 bmi if bmi >= 19.8 & bmi <=39.3) , ytitle("log PSA") legend(label(1 "Krimpen") label(2 "PCPT") label(3 "PLCO") label(4 "ProtecT")) 
graph export "$cd_graphs\IPD - BMI-logPSA splines.tif", as(tif) replace

*BMI - PCa
use "IPD All clear.dta", clear

sort study bmi
gen logpsa_age = .
gen study_i = .
gen se = .
gen beta = .
foreach i in 2 3 4 5 {
	bspline if study == `i', xvar(bmi) knots(0 24.8 27 30 100) p(3) gen(_study_`i')
	qui logit pca _study_`i'*, noconstant
	predict study_`i'
	qui logit pca bmi if study == `i'
	qui replace study_i = `i' in `i'
	qui replace beta = _b[bmi] in `i'
	qui replace se = _se[bmi] in `i'	
}

twoway (line study_2 bmi if bmi >= 19.8 & bmi <=34.5) (line study_3 bmi if bmi >= 20.4 & bmi <=41) (line study_4 bmi if bmi >= 19.9 & bmi <=40.5) (line study_5 bmi if bmi >= 19.8 & bmi <=39.3) , ytitle("Prostate Cancer Risk") legend(label(1 "Krimpen") label(2 "PCPT") label(3 "PLCO") label(4 "ProtecT")) 
graph export "$cd_graphs\IPD - BMI-PCa splines.tif", as(tif) replace

*PSA - PCa
use "IPD All clear.dta", clear

sort study psa
gen study_i = .
gen se = .
gen beta = .
foreach i in 2 3 4 5 {
	if `i' == 2 {
		bspline if study == `i' & psa >= 0.2 & psa <= 13.3, xvar(logpsa) knots(-6 -0.5 0 0.5 8) p(3) gen(_study_`i')
	}
	if `i' == 3 {
		bspline if study == `i' & psa >= 0.3 & psa <=3, xvar(logpsa) knots(-6 -0.5 0 0.5 8) p(3) gen(_study_`i')
	}	
	if `i' == 4 {
		bspline if study == `i' & psa >= 4, xvar(logpsa) knots(-6 -0.5 0 0.5 8) p(3) gen(_study_`i')
	}
	if `i' == 5 {
		bspline if study == `i' & psa >= 3, xvar(logpsa) knots(-6 -0.5 0 0.5 8) p(3) gen(_study_`i')
	}
	qui logit pca _study_`i'*, noconstant
	predict study_`i'
	drop _study_`i'*
	qui logit pca logpsa if study == `i'
	qui replace study_i = `i' in `i'
	qui replace beta = _b[logpsa] in `i'
	qui replace se = _se[logpsa] in `i'	
}

twoway (line study_2 logpsa if logpsa >= -1.6 & logpsa <= 2.4) (line study_3 logpsa if logpsa >= -1.6 & logpsa <= 2.4) (line study_4 logpsa if psa >= 4 & logpsa <= 2.4) (line study_5 logpsa if psa >= 3 & logpsa <= 2.4), xtitle("Log-PSA") ytitle("Prostate Cancer Risk") legend(label(1 "Krimpen") label(2 "PCPT") label(3 "PLCO") label(4 "ProtecT")) xline(1.099 1.386)
graph export "$cd_graphs\IPD - PSA-PCa splines.tif", as(tif) replace width(1200)

********************************************************************************

*Cubic splines (imputed)
*BMI - PSA
use "Imputed data (linear).dta", clear

keep if _mi_m == 1
replace psa = exp(logpsa)
sort study logpsa
gen study_i = .
gen se = .
gen beta = .
foreach i in 2 3 4 5 {
	if `i' == 2 {
		bspline if study == `i' , xvar(logpsa) knots(-6 -0.5 0 0.5 8) p(3) gen(_study_`i')
	}
	if `i' == 3 {
		bspline if study == `i' , xvar(logpsa) knots(-6 -0.5 0 0.5 8) p(3) gen(_study_`i')
	}	
	if `i' == 4 {
		bspline if study == `i' , xvar(logpsa) knots(-6 -0.5 0 0.5 8) p(3) gen(_study_`i')
	}
	if `i' == 5 {
		bspline if study == `i' , xvar(logpsa) knots(-6 -0.5 0 0.5 8) p(3) gen(_study_`i')
	}
	qui logit pca _study_`i'*, noconstant
	predict study_`i'
	drop _study_`i'*
	qui logit pca logpsa if study == `i'
	qui replace study_i = `i' in `i'
	qui replace beta = _b[logpsa] in `i'
	qui replace se = _se[logpsa] in `i'	
}

twoway (line study_2 logpsa if logpsa >= -1.6 & logpsa <= 2.4) (line study_3 logpsa if logpsa >= -1.2 & logpsa <= 1.1) (line study_4 logpsa if logpsa >= -1.8 & logpsa <= 2.4) (line study_5 logpsa if logpsa >= -1.6 & logpsa <= 2.4), xtitle("Log-PSA") ytitle("Prostate Cancer Risk") legend(label(1 "Krimpen") label(2 "PCPT") label(3 "PLCO") label(4 "ProtecT")) xline(1.099 1.386)
graph export "$cd_graphs\IPD - PSA-PCa splines.tif (IMPUTED).tif", as(tif) replace width(1200)

********************************************************************************

*Forest plots for appendix
use "IPD PCa - complete", clear
replace type = "Imputed"
append using "Complete case - pca"
replace outcome = "All" if outcome == ""
append using "Complete case - advanced pca"
replace outcome = "Advanced" if outcome == ""
drop if study == "PCPT"
replace study = studyname if study == ""
gen type2 = 0
replace type2 = 1 if type == "Imputed"
label define type 0 "Not Imputed" 1 "Imputed"
label values type2 type
sort type2 study

metan md5 md5_se if outcome == "All", by(type2) label(namevar = study) eform nooverall xlabel(0.85, 0.9, 0.95, 1.05) xtitle("OR for PCa per 5 kg/m{superscript:2} increase BMI", size(vsmall)) force xline(-0.0576, lpattern(shortdash) lcolor(maroon) lwidth(thin)) random xline(-0.005, lpattern(shortdash) lcolor(maroon) lwidth(thin))
graph export "$cd_graphs\Appendix IPD PCa.tif", replace as(tif) width(1200)

metan md5 md5_se if outcome == "Advanced", by(type2) label(namevar = study) eform nooverall xlabel(0.8, 0.9, 1.1, 1.2, 1.3) xtitle("OR for advanced PCa per 5 kg/m{superscript:2} increase BMI", size(vsmall)) force xline(-0.0171, lpattern(shortdash) lcolor(maroon) lwidth(thin)) xline(0.0159, lpattern(shortdash) lcolor(maroon) lwidth(thin)) random
graph export "$cd_graphs\Appendix IPD aPCa.tif", replace as(tif) width(1200)

********************************************************************************

*Table A1

use "IPD All clear.dta", clear
gen variable = ""
replace variable = "N" in 2
replace variable = "Missing data" in 3
replace variable = "PCa (N [%])" in 4
replace variable = "BMI (N [%])" in 5
replace variable = "Log-PSA (N [%])" in 6
replace variable = "Family history PCa (N [%])" in 7
replace variable = "Missing data for:" in 8
replace variable = "0 variables" in 9
replace variable = "1 variable" in 10
replace variable = "2 variables" in 11
replace variable = "3 variables" in 12
replace variable = "4 variables" in 13

forvalues i = 1/5 {
	gen v`i' = ""
}
replace v2 = "Krimpen" in 1
replace v3 = "PCPT" in 1
replace v4 = "PLCO" in 1
replace v5 = "ProtecT" in 1
replace v1 = "Total for Analysis" in 1

order variable v2 v4 v5 v1 v3

local i = 2
foreach study in 2 3 4 5 {
	qui count if study == `study'
	local value = r(N)
	local value: dis %6.0fc `value'
	local value = subinstr("`value'"," ","",.)
	qui replace v`study' = "`value'" in 2
}
local value = c(N)
local value: dis %6.0fc `value'
local value = subinstr("`value'"," ","",.)
qui replace v1 = "`value'" in 2

local i = 4
foreach var of varlist pca bmi logpsa familyhistory {
	foreach study in 2 3 4 5 {
		qui count if `var' == . & study == `study'
		local n = r(N)
		qui count if study == `study'
		local N = r(N)
		
		local percent = `n'*100/`N'
		local percent: dis %3.1f `percent'
		
		local n: dis %6.0fc `n'
		local n = subinstr("`n'"," ","",.)
		
		local value = "`n' (`percent')"
		qui replace v`study' = "`value'" in `i'

	}
	qui count if `var' == . & study != 3
	local n = r(N)
	qui count if study != 3
	local N = r(N)
	
	local percent = `n'*100/`N'
	local percent: dis %3.1f `percent'
	
	local n: dis %6.0fc `n'
	local n = subinstr("`n'"," ","",.)
	
	local value = "`n' (`percent')"
	qui replace v1 = "`value'" in `i'	
	
	local i = `i'+1
	gen `var'_miss = 1 if `var' == .
	replace `var'_miss = 0 if `var'_miss == .
}

gen x = pca_miss + bmi_miss + logpsa_miss + familyhistory_miss

local i = 9
forvalues k = 0/4 {
	foreach study in 2 3 4 5 {

		qui count if x == `k' & study == `study'
		local n = r(N)
		qui count if study == `study'
		local N = r(N)
		
		local percent = `n'*100/`N'
		local percent: dis %3.1f `percent'
		
		local n: dis %6.0fc `n'
		local n = subinstr("`n'"," ","",.)
		
		local value = "`n' (`percent')"
		qui replace v`study' = "`value'" in `i'

	}
	qui count if x == `k' & study != 3
	local n = r(N)
	qui count if study != 3
	local N = r(N)
	
	local percent = `n'*100/`N'
	local percent: dis %3.1f `percent'
	
	local n: dis %6.0fc `n'
	local n = subinstr("`n'"," ","",.)
	
	local value = "`n' (`percent')"
	qui replace v1 = "`value'" in `i'
	
	local i = `i'+1
}

keep variable-v3
keep in 1/`i'

save "$cd_tables\Table A1.dta", replace

**********

*Table A3
use "Imputed data (linear).dta", clear
gen variable = ""
replace variable = "N" in 2
replace variable = "Men who were biopsied (prostate cancer status not imputed)" in 3
replace variable = "Participants (% of total)" in 4
replace variable = "PCa (%)" in 5
replace variable = "Advanced PCa (%)*" in 6
replace variable = "No PCa (%)" in 7
replace variable = "Men who were not biopsied (prostate cancer status imputed)" in 8
replace variable = "Participants (% of total)" in 9
replace variable = "PCa (%)" in 10
replace variable = "Advanced PCa (%)*" in 11
replace variable = "No PCa (%)" in 12
replace variable = "All participants (men who were and were not biopsied combined)" in 13
replace variable = "Participants" in 14
replace variable = "PCa (%)" in 15
replace variable = "Advanced PCa (%)*" in 16
replace variable = "No PCa (%)" in 17
replace variable = "Age (mean, [SD])" in 18
replace variable = "BMI (mean, [SD])" in 19
replace variable = "Log-PSA (mean, [SD])" in 20
replace variable = "Family history PCa (%)" in 21


forvalues i = 1/5 {
	gen v`i' = ""
}
replace v2 = "Krimpen" in 1
replace v3 = "PCPT" in 1
replace v4 = "PLCO" in 1
replace v5 = "ProtecT" in 1
replace v1 = "Total for Analysis" in 1

order variable v2 v4 v5 v1 v3

local i = 2
foreach study in 1 2 3 4 5 {
	if `study' == 1 {
		qui count if _mi_m == 0
	}
	else {
		qui count if study == `study' & _mi_m == 0
	}
	local value = r(N)
	local value: dis %6.0fc `value'
	local value = subinstr("`value'"," ","",.)
	qui replace v`study' = "`value'" in 2
}

*Complete case

foreach study in 1 2 3 4 5 {
	local i = 4
	*N
	if `study' == 1 {
		qui count if pca != . & _mi_m == 0 & study != 3
	}
	else {
		qui count if pca != . & study == `study' & _mi_m == 0
	}
	local n = r(N)
	
	if `study' == 1 {
		qui count if _mi_m == 0 & study != 3
	}
	else {
		qui count if study == `study' & _mi_m == 0
	}
	local N = r(N)
	
	local percent = `n'*100/`N'
	local percent: dis %3.1f `percent'
	
	local n: dis %6.0fc `n'
	local n = subinstr("`n'"," ","",.)
	
	local value = "`n' (`percent')"
	qui replace v`study' = "`value'" in `i'
	
	local i = `i' + 1
	
	*PCa
	if `study' == 1 {
		qui count if pca == 1 & _mi_m == 0 & study != 3
	}
	else {
		qui count if pca == 1 & study == `study' & _mi_m == 0
	}
	local n = r(N)
	
	if `study' == 1 {
		qui count if pca != . & _mi_m == 0 & study != 3
	}
	else {
		qui count if pca != . & study == `study' & _mi_m == 0
	}
	local N = r(N)
	
	local percent = `n'*100/`N'
	local percent: dis %3.1f `percent'
	
	local n: dis %6.0fc `n'
	local n = subinstr("`n'"," ","",.)
	
	local value = "`n' (`percent')"
	qui replace v`study' = "`value'" in `i'
	
	local i = `i' + 1
	
	*aPCa
	if `study' == 1 {
		qui count if advanced_pca == 1 & _mi_m == 0 & study != 3
	}
	else {
		qui count if advanced_pca == 1 & study == `study' & _mi_m == 0
	}
	local n = r(N)
	
	if `study' == 1 {
		qui count if advanced_pca != . & _mi_m == 0 & study != 3
	}
	else {
		qui count if advanced_pca != . & study == `study' & _mi_m == 0
	}
	local N = r(N)
	
	local percent = `n'*100/`N'
	local percent: dis %3.1f `percent'
	
	local n: dis %6.0fc `n'
	local n = subinstr("`n'"," ","",.)
	
	local value = "`n' (`percent')"
	qui replace v`study' = "`value'" in `i'
	
	local i = `i' + 1
	
	*No PCa
	if `study' == 1 {
		qui count if pca == 0 & _mi_m == 0 & study != 3
	}
	else {
		qui count if pca == 0 & study == `study' & _mi_m == 0
	}
	local n = r(N)
	
	if `study' == 1 {
		qui count if pca != . & _mi_m == 0 & study != 3
	}
	else {
		qui count if pca != . & study == `study' & _mi_m == 0
	}
	local N = r(N)
	
	local percent = `n'*100/`N'
	local percent: dis %3.1f `percent'
	
	local n: dis %6.0fc `n'
	local n = subinstr("`n'"," ","",.)
	
	local value = "`n' (`percent')"
	qui replace v`study' = "`value'" in `i'
	
	local i = `i' + 1
}

*Imputed only
	
foreach study in 1 2 3 4 5 {
	local i = 9
	*N
	if `study' == 1 {
		qui count if pca != . & _mi_m != 0 & study != 3 & biopsy == 0
	}
	else {
		qui count if pca != . & study == `study' & _mi_m != 0 & biopsy == 0
	}
	local n = r(N)/100
	
	if `study' == 1 {
		qui count if _mi_m != 0 & study != 3 
	}
	else {
		qui count if study == `study' & _mi_m != 0
	}
	local N = r(N)/100
	
	local percent = `n'*100/`N'
	local percent: dis %3.1f `percent'
	
	local n: dis %6.0fc `n'
	local n = subinstr("`n'"," ","",.)
	
	local value = "`n' (`percent')"
	qui replace v`study' = "`value'" in `i'
	
	local i = `i' + 1
	
	*PCa
	if `study' == 1 {
		qui count if pca == 1 & _mi_m != 0 & study != 3 & biopsy == 0
	}
	else {
		qui count if pca == 1 & study == `study' & _mi_m != 0 & biopsy == 0
	}
	local n = r(N)/100
	
	if `study' == 1 {
		qui count if pca != . & _mi_m != 0 & study != 3 & biopsy == 0
	}
	else {
		qui count if pca != . & study == `study' & _mi_m != 0 & biopsy == 0
	}
	local N = r(N)/100
	
	local percent = `n'*100/`N'
	local percent: dis %3.1f `percent'
	
	local n: dis %6.0fc `n'
	local n = subinstr("`n'"," ","",.)
	
	local value = "`n' (`percent')"
	qui replace v`study' = "`value'" in `i'
	
	local i = `i' + 1
	
	*aPCa
	*Needs the advanced imputation data	
	local i = `i' + 1
	
	*No PCa
	if `study' == 1 {
		qui count if pca == 0 & _mi_m != 0 & study != 3 & biopsy == 0
	}
	else {
		qui count if pca == 0 & study == `study' & _mi_m != 0 & biopsy == 0
	}
	local n = r(N)/100
	
	if `study' == 1 {
		qui count if pca != . & _mi_m != 0 & study != 3 & biopsy == 0
	}
	else {
		qui count if pca != . & study == `study' & _mi_m != 0 & biopsy == 0
	}
	local N = r(N)/100
	
	local percent = `n'*100/`N'
	local percent: dis %3.1f `percent'
	
	local n: dis %6.0fc `n'
	local n = subinstr("`n'"," ","",.)
	
	local value = "`n' (`percent')"
	qui replace v`study' = "`value'" in `i'
	
	local i = `i' + 1
}

*Imputed (all)
	
foreach study in 1 2 3 4 5 {
	local i = 14
	*N
	if `study' == 1 {
		qui count if pca != . & _mi_m != 0 & study != 3 
	}
	else {
		qui count if pca != . & study == `study' & _mi_m != 0 
	}
	local n = r(N)/100
	
	if `study' == 1 {
		qui count if _mi_m != 0 & study != 3
	}
	else {
		qui count if study == `study' & _mi_m != 0
	}
	local N = r(N)/100
	
	local percent = `n'*100/`N'
	local percent: dis %3.1f `percent'
	
	local n: dis %6.0fc `n'
	local n = subinstr("`n'"," ","",.)
	
	local value = "`n' (`percent')"
	qui replace v`study' = "`value'" in `i'
	
	local i = `i' + 1
	
	*PCa
	if `study' == 1 {
		qui count if pca == 1 & _mi_m != 0 & study != 3 
	}
	else {
		qui count if pca == 1 & study == `study' & _mi_m != 0 
	}
	local n = r(N)/100
	
	if `study' == 1 {
		qui count if pca != . & _mi_m != 0 & study != 3 
	}
	else {
		qui count if pca != . & study == `study' & _mi_m != 0 
	}
	local N = r(N)/100
	
	local percent = `n'*100/`N'
	local percent: dis %3.1f `percent'
	
	local n: dis %6.0fc `n'
	local n = subinstr("`n'"," ","",.)
	
	local value = "`n' (`percent')"
	qui replace v`study' = "`value'" in `i'
	
	local i = `i' + 1
	
	*aPCa
	*Needs the advanced imputation data	
	local i = `i' + 1
	
	*No PCa
	if `study' == 1 {
		qui count if pca == 0 & _mi_m != 0 & study != 3 
	}
	else {
		qui count if pca == 0 & study == `study' & _mi_m != 0 
	}
	local n = r(N)/100
	
	if `study' == 1 {
		qui count if pca != . & _mi_m != 0 & study != 3 
	}
	else {
		qui count if pca != . & study == `study' & _mi_m != 0 
	}
	local N = r(N)/100
	
	local percent = `n'*100/`N'
	local percent: dis %3.1f `percent'
	
	local n: dis %6.0fc `n'
	local n = subinstr("`n'"," ","",.)
	
	local value = "`n' (`percent')"
	qui replace v`study' = "`value'" in `i'
	
	local i = `i' + 1
}

*Age, BMI, log(PSA)
local i = 18
foreach var of varlist age bmi logpsa {
	foreach study in 1 2 3 4 5 {
		if `study' == 1 {
			qui su `var' if _mi_m != 0 & study != 3 
		}
		else {
			qui su `var' if _mi_m != 0 &study == `study'
		}
		local mean = r(mean)
		local sd = r(sd)

		if "`var'" == "logpsa" {
			local mean: dis %3.2f `mean'
		}
		else {
			local mean: dis %3.1f `mean'
		}		

		if "`var'" == "logpsa" {
			local sd: dis %3.2f `sd'
		}
		else {
			local sd: dis %2.1f `sd'
		}
		
		local value = "`mean' (`sd')"
		qui replace v`study' = "`value'" in `i'
		
	}
	local i = `i'+1
}

*Family history
local i = 21
foreach study in 1 2 3 4 5 {
	if `study' == 1 {
		qui count if familyhistory == 1 & _mi_m != 0 & study != 3 
	}
	else {
		qui count if familyhistory == 1 & study == `study' & _mi_m != 0 
	}
	local n = r(N)/100
	
	if `study' == 1 {
		qui count if _mi_m != 0 & study != 3
	}
	else {
		qui count if study == `study' & _mi_m != 0
	}
	local N = r(N)/100
	
	local percent = `n'*100/`N'
	local percent: dis %3.1f `percent'
	
	local n: dis %6.0fc `n'
	local n = subinstr("`n'"," ","",.)
	
	local value = "`n' (`percent')"
	qui replace v`study' = "`value'" in `i'
}

keep variable - v3
keep in 1/21

save "$cd_tables\Table A3.dta", replace 

*Advanced PCa
use "Imputed data (linear, advanced).dta", clear

forvalues i = 1/5 {
	gen v`i' = ""
}

order v2 v4 v5 v1 v3

forvalues study = 1/5 {
	local i = 1
	*Imputed only
	if `study' == 1 {
		qui count if advanced_pca == 1 & _mi_m != 0 & study != 3 & biopsy == 0
	}
	else {
		qui count if advanced_pca == 1 & study == `study' & _mi_m != 0 & biopsy == 0
	}
	local n = r(N)/100
	
	if `study' == 1 {
		qui count if advanced_pca != . & _mi_m == 0 & study != 3 & biopsy == 0
	}
	else {
		qui count if advanced_pca != . & study == `study' & _mi_m != 0 & biopsy == 0
	}
	local N = r(N)/100
	
	local percent = `n'*100/`N'
	local percent: dis %3.1f `percent'
	
	local n: dis %6.0fc `n'
	local n = subinstr("`n'"," ","",.)
	
	local value = "`n' (`percent')"
	qui replace v`study' = "`value'" in `i'
	
	local i = `i' + 1
	
	*Imputed (all)
	if `study' == 1 {
		qui count if advanced_pca == 1 & _mi_m != 0 & study != 3 
	}
	else {
		qui count if advanced_pca == 1 & study == `study' & _mi_m != 0 
	}
	local n = r(N)/100
	
	if `study' == 1 {
		qui count if advanced_pca != . & _mi_m == 0 & study != 3
	}
	else {
		qui count if advanced_pca != . & study == `study' & _mi_m != 0
	}
	local N = r(N)/100
	
	local percent = `n'*100/`N'
	local percent: dis %3.1f `percent'
	
	local n: dis %6.0fc `n'
	local n = subinstr("`n'"," ","",.)
	
	local value = "`n' (`percent')"
	qui replace v`study' = "`value'" in `i'
	
	local i = `i' + 1
}

keep in 1/2
keep v2 - v3

save "$cd_tables\Table A3 advanced.dta", replace

use "$cd_tables\Table A3.dta", clear 
append using "$cd_tables\Table A3 advanced.dta"
forvalues i = 1/5 {
	qui replace v`i' = v`i'[22] in 11
	qui replace v`i' = v`i'[23] in 16
}
drop in 22/23
save "$cd_tables\Table A3.dta", replace 

********************************************************************************

*IPD analysis complete, ready to add to studies from systematic review

********************************************************************************

}

}
********************************************************************************
********************************************************************************

*PSA Analyses
{
{
set more off
clear
use "Original data\AD\PSA MD.dta", clear

*Type = Beta
gen type = "Beta"

*qui replace the authors, if there is more than 1, with just the primary author
qui egen author_1 = ends(author), punct(" and ") trim head
qui egen author_primary = ends(author_1), punct(,) trim head
qui egen author_primary_2 = ends(author_primary) if strpos(author_primary,".")>0, trim tail
qui replace author_primary = author_primary_2 if author_primary_2 !=""
qui egen author_primary_3 = ends(author_primary) if strpos(author_primary_2,".")>0, trim tail
qui replace author_primary = author_primary_3 if author_primary_3 !=""
qui egen author_primary_4 = ends(author_primary) if strpos(author_primary_3,".")>0, trim tail
qui replace author_primary = author_primary_4 if author_primary_4 !=""
drop author_primary_2 author_primary_3 author_primary_4 author_1 
order author_primary, b(author)
rename author author_all
rename author_primary author
order author_all, last

*PSA should be on the beta scale
qui gen beta = .
qui gen beta_se = .
qui gen study_bmi_mean = .

*Total number of participants
qui gen ntotal = ncases
qui replace ntotal = ncontrols if ncases == .

*All studies reporting SD were actually SE, so change this (not for the studies with controls only, that *is* SD)
qui replace sdcases = sdcases*sqrt(ncases)

*Split 16237 into 3 different studies
qui replace endnoteid = 162371 if strpos(group, "SEARCH")>0
qui replace endnoteid = 162372 if strpos(group, "Duke")>0
qui replace endnoteid = 162373 if strpos(group, "Johns")>0
qui replace group = subinstr(group, "SEARCH ","",.)
qui replace group = subinstr(group, "Duke ","",.)
qui replace group = subinstr(group, "Johns ","",.)
qui replace bmisd = 4.9 if endnoteid == 162371
qui replace bmisd = 4.4 if endnoteid == 162372
qui replace bmisd = 3.3 if endnoteid == 162373
qui replace note = "" if endnoteid == 162371 | endnoteid == 162372 | endnoteid == 162373 
qui replace x = x-4 if endnoteid == 162372
qui replace x = x-8 if endnoteid == 162373
qui replace studyname = "SEARCH" if endnoteid == 162371
qui replace studyname = "Duke" if endnoteid == 162372
qui replace studyname = "Johns Hopkins" if endnoteid == 162373

*Also, no real point in having anything other than "participants" rather than cases/controls in this analysis
qui replace ncases = ncontrols if ncases ==.
qui replace sdcases = sdcontrols if sdcases ==.
qui replace meanmediancases = meanmediancontrols if meanmediancases == .

*Calculate SD from 95% CI
qui egen lcicase = ends(cicases), punct("-") trim head
destring lcicase, replace
qui egen ucicase = ends(cicases), punct("-") trim tail
destring ucicase, replace
qui replace sdcases = sqrt(ncases)*(ucicase-lcicase)/3.92 if sdcases ==.

qui egen lcicontrols = ends(cicontrols), punct("-") trim head
destring lcicontrols, replace
qui egen ucicontrols = ends(cicontrols), punct("-") trim tail
destring ucicontrols, replace
qui replace sdcases = sqrt(ncontrols)*(ucicontrols-lcicontrols)/3.92 if sdcases ==. & endnote != 214346

*13545 needs a different SD as it is geometric
qui replace sdcases = exp(sqrt(ncases)*(ln(ucicase)-ln(lcicase))/3.92) if endnote == 13545

drop cicases cicontrols lcicase ucicase lcicontrols ucicontrols
sort endnote x

*Calculate bmi_sd from quantile data
*Note: this is tricky and from Chene and Thompson
*Easiest way is to probably cycle through studies and exit if not quantile

*Destring the quantiles
qui replace bmisubgroup = group if bmisubgroup == "" & group != "All"
qui egen bmi_lower = ends(bmisubgroup), punct("-") trim head
qui egen bmi_lower2 = ends(bmisubgroup), punct(">") trim tail
replace  bmi_lower = bmi_lower2 if bmi_lower2 !=""
replace bmi_lower = "" if strpos(bmi_lower, "<")>0
replace bmi_lower = subinstr(bmi_lower,"=","",.)
destring bmi_lower, replace
qui egen bmi_upper = ends(bmisubgroup), punct("-") trim tail
qui egen bmi_upper2 = ends(bmisubgroup), punct("<") trim tail
replace  bmi_upper = bmi_upper2 if bmi_upper2 !=""
destring bmi_upper, replace

drop bmi_upper2 bmi_lower2

*Generate useful variables
sort endnoteid x
qui egen study = group(endnoteid outcome)
bysort study: egen count = count(study)
qui replace count = . if count <=2

qui gen Xj = .
qui replace Xj = bmi_upper

qui gen Nj = ntotal
qui gen Pj = .
qui gen Zj = .
qui gen Wj = .
qui gen Mj = .

sort study x
qui sum study
local max = r(max)
forvalues study = 1/`max' {
	qui sum count if study == `study'
	if r(mean) >2 & r(mean) <. {
		qui sum x if study == `study'
		local max_x = r(max)
		
		*Total number of participants in the study
		qui sum Nj if study==`study'
		local sum = r(sum)
		
		forvalues value = 1/`max_x' {
			qui sum Nj if x <= `value' & study==`study'
			qui replace Pj = r(sum)/`sum' if x==`value' & study==`study'
		}

		qui replace Zj = invnormal(Pj)
		qui replace Wj = normalden(Zj,0,1)^2/(Pj*(1-Pj))
	
		qui sum x if study==`study'
		local quantile = r(max)-1
		qui sum Wj if study==`study'
		qui replace Wj = Wj*(`quantile'/r(sum))

		*Regress to get the mean and SD of the study, using the beta of Zj for the mean and the constant for the SD
		qui regress Xj Zj if study==`study'
		local m = _b[_cons]
		qui replace study_bmi_mean = _b[_cons] if study == `study'
		local s = _b[Zj]
		qui replace bmisd = _b[Zj] if study == `study' & bmisd == .
	
		forvalues value = 1/`max_x' {
			if `value' == `max_x' {
				local zj1 = 9999999
			}
			else {
				qui sum Zj if x == `value' & study == `study'
				local zj1 = r(mean)
			}
			if `value' == 1 {
				local zj0 = -999999
			}
			else {
				local value0 = `value'-1
				qui sum Zj if x == `value0' & study == `study'
				local zj0 = r(mean)
			}
			qui replace Mj = `m'+`s'*((normalden(`zj0')-normalden(`zj1'))/(normal(`zj1')-normal(`zj0'))) if study == `study' & x == `value'
		}
	
	*drop if x>1 & study == `study'
	}
}
qui replace bmimean = Mj if bmimean == .
qui drop count-Mj group outcome author_all prospectiveretro* timeof meanyears studytype outcomepsa journal title howare psathres dre bmisdcase bmisdcontrol ncases
order bmisd, a(study_bmi_mean)
order notes, last
rename bmisd study_bmi_sd
rename meanmediancases psa_mean
rename sdcases psa_sd

*Get study BMI for the study with just high/low
qui sum ntotal if endnote == 4480 & x == 1
local n1 = r(mean)
qui sum ntotal if endnote == 4480 & x == 2
local n2 = r(mean)
qui sum bmimean if endnote == 4480 & x == 1
local b1 = r(mean)
qui sum bmimean if endnote == 4480 & x == 2
local b2 = r(mean)
replace study_bmi_mean = (`n1'*`b1'+`n2'*`b2')/(`n1'+`n2') if endnote == 4480

*Generate SE for PSA Mean
qui gen psa_mean_se = psa_sd/sqrt(ntotal), after(psa_mean)

*Clean up BMI_lower
replace bmi_lower = 0 if bmi_lower == .
replace bmi_upper = 100 if bmi_upper == .

*PSA needs to be log-PSA
*First geometric means
gen logpsa = log(psa_mean) if endnote == 13545 | endnote == 36677 | endnote == 13747 | endnote == 214346
gen logpsa_sd = log(psa_sd) if endnote == 13545 
replace logpsa_sd = psa_sd if endnote == 36677 | endnote == 13747 
replace logpsa_sd = psa_sd*sqrt(ntotal) if endnote == 214346

*Turn the IQR into an SD - 30553 & 209795
qui egen iqr_lower = ends(iqrcontrols), punct("-") trim head
destring iqr_lower, replace
qui egen iqr_upper = ends(iqrcontrols), punct("-") trim tail
destring iqr_upper, replace
replace psa_sd = (iqr_upper-iqr_lower)/1.35 if endnote == 30553 | endnote == 209795
drop iqr_lower-iqr_upper iqrcases-meandiffere
order logpsa-logpsa_sd, a(type)

*This is also the only median, so get logpsa for 30553
replace logpsa = log(psa_mean) if endnote == 30553
replace logpsa_sd = sqrt(log((1+sqrt(1+4*psa_sd^2/psa_mean^2))/2)) if endnote == 30553

*Everything else should be mean and SD
replace logpsa = log(psa_mean^2/sqrt(psa_sd^2+psa_mean^2)) if logpsa ==.
replace logpsa_sd = sqrt(log(psa_sd^2/psa_mean^2+1)) if logpsa_sd ==.
gen logpsa_se = logpsa_sd/sqrt(ntotal)
order logpsa_se, a(logpsa_sd)

*VWLS - linear
sort study x
qui sum study
local max = r(max)
forvalues study = 1/`max' {
	dis "`study'"
	vwls logpsa bmimean if study == `study', sd(logpsa_se)
	qui replace beta = _b[bmimean] if study == `study'
	qui replace beta_se = _se[bmimean] if study == `study'
	qui replace type = "Categories" if study == `study'
}

*Categorise into normal, overweight and obese
gen bmi_group = 1 if bmi_upper <= 25, a(bmi_upper)
replace bmi_group = 2 if bmi_upper <=30 & bmi_group == .
replace bmi_group = 3 if bmi_group == .
label define bmi_group 1 "Normal" 2 "Overweight" 3 "Obese"
label values bmi_group bmi_group

*JUST look at normal = 18.5-25, overweight = 25-30 and obese = 30-35
*If normal has two categories, merge them (assume similar means)
*Studies with 2 lower categories
qui sum logpsa if bmi_lower == 0 & endnote == 13545
local mean1 = r(mean)
qui sum logpsa_sd if bmi_lower == 0 & endnote == 13545
local sd1 = r(mean)
qui sum ntotal if bmi_lower == 0 & endnote == 13545
local n1 = r(mean)
qui sum bmimean if bmi_lower == 0 & endnote == 13545
local bmi1 = r(mean)
qui sum logpsa if bmi_lower == 23 & endnote == 13545
local mean2 = r(mean)
qui sum logpsa_sd if bmi_lower == 23 & endnote == 13545
local sd2 = r(mean)
qui sum ntotal if bmi_lower == 23 & endnote == 13545
local n2 = r(mean)
qui sum bmimean if bmi_lower == 23 & endnote == 13545
local bmi2 = r(mean)

local mean = (`mean1'*`n1'+`mean2'*`n2')/(`n1'+`n2')
local sd = sqrt(((`n1'-1)*`sd1'^2+(`n2'-1)*`sd2'^2)/(`n1'+`n2'-2))
local bmi = (`bmi1'*`n1'+`bmi2'*`n2')/(`n1'+`n2')

qui replace logpsa = `mean' if bmi_lower == 23 & endnote == 13545
qui replace logpsa_sd = `sd' if bmi_lower == 23 & endnote == 13545
qui replace ntotal = `n1'+`n2' if bmi_lower == 23 & endnote == 13545
qui replace bmimean = `bmi' if bmi_lower == 23 & endnote == 13545

*Studies with 4 BMI categories (BMI>30 & BMI>35)
foreach endnote in 14033 15526 18448 162371 162372 162373 {

	qui sum logpsa if bmi_lower == 30 & endnote == `endnote'
	local mean1 = r(mean)
	qui sum logpsa_sd if bmi_lower == 30 & endnote == `endnote'
	local sd1 = r(mean)
	qui sum ntotal if bmi_lower == 30 & endnote == `endnote'
	local n1 = r(mean)
	qui sum bmimean if bmi_lower == 30 & endnote == `endnote'
	local bmi1 = r(mean)
	qui sum logpsa if bmi_lower == 35 & endnote == `endnote' 
	local mean2 = r(mean)
	qui sum logpsa_sd if bmi_lower == 35 & endnote == `endnote' 
	local sd2 = r(mean)
	qui sum ntotal if bmi_lower == 35 & endnote == `endnote' 
	local n2 = r(mean)
	qui sum bmimean if bmi_lower == 35 & endnote == `endnote'
	local bmi2 = r(mean)

	local mean = (`mean1'*`n1'+`mean2'*`n2')/(`n1'+`n2')
	local sd = sqrt(((`n1'-1)*`sd1'^2+(`n2'-1)*`sd2'^2)/(`n1'+`n2'-2))
	local bmi = (`bmi1'*`n1'+`bmi2'*`n2')/(`n1'+`n2')

	qui replace logpsa = `mean' if bmi_lower == 30 & endnote == `endnote' 
	qui replace logpsa_sd = `sd' if bmi_lower == 30 & endnote == `endnote' 
	qui replace ntotal = `n1'+`n2' if bmi_lower == 30 & endnote == `endnote'
	qui replace bmimean = `bmi' if bmi_lower == 30 & endnote == `endnote'
	
	drop if bmi_lower == 35 & endnote == `endnote'
}

*And 5 categories
qui sum logpsa if bmi_lower == 30 & endnote == 19431
local mean1 = r(mean)
qui sum logpsa_sd if bmi_lower == 30 & endnote == 19431
local sd1 = r(mean)
qui sum ntotal if bmi_lower == 30 & endnote == 19431 
local n1 = r(mean)
qui sum bmimean if bmi_lower == 30 & endnote == 19431
local bmi1 = r(mean)

qui sum logpsa if bmi_lower == 35 & endnote == 19431
local mean2 = r(mean)
qui sum logpsa_sd if bmi_lower == 35 & endnote == 19431
local sd2 = r(mean)
qui sum ntotal if bmi_lower == 35 & endnote == 19431
local n2 = r(mean)
qui sum bmimean if bmi_lower == 35 & endnote == 19431
local bmi2 = r(mean)
	
qui sum logpsa if bmi_lower == 40 & endnote == 19431 
local mean3 = r(mean)
qui sum logpsa_sd if bmi_lower == 40 & endnote == 19431
local sd3 = r(mean)
qui sum ntotal if bmi_lower == 40 & endnote == 19431
local n3 = r(mean)
qui sum bmimean if bmi_lower == 40 & endnote == 19431
local bmi3 = r(mean)

local mean = (`mean1'*`n1'+`mean2'*`n2'+`mean3'*`n3')/(`n1'+`n2'+`n3')
local sd = sqrt(((`n1'-1)*`sd1'^2+(`n2'-1)*`sd2'^2+(`n3'-1)*`sd3'^2)/(`n1'+`n2'+`n3'-2))
local bmi = (`bmi1'*`n1'+`bmi2'*`n2'+`bmi3'*`n3')/(`n1'+`n2'+`n3')

qui replace logpsa = `mean' if bmi_lower == 30 & endnote == 19431
qui replace logpsa_sd = `sd' if bmi_lower == 30 & endnote == 19431
qui replace ntotal = `n1'+`n2'+`n3' if bmi_lower == 30 & endnote == 19431
qui replace bmimean = `bmi' if bmi_lower == 30 & endnote == 19431

drop if bmi_lower >= 35 & endnote == 19431

*And studies with >1 "Normal" BMI
drop if bmi_lower == 0 & endnote == 13545

foreach endnote in 211732 212495 214528 {

	qui sum logpsa if x==1 & endnote == `endnote'
	local mean1 = r(mean)
	qui sum logpsa_sd if x==1 & endnote == `endnote'
	local sd1 = r(mean)
	qui sum ntotal if x==1 & endnote == `endnote'
	local n1 = r(mean)
	qui sum bmimean if x==1 & endnote == `endnote'
	local bmi1 = r(mean)
	qui sum logpsa if x==2 & endnote == `endnote' 
	local mean2 = r(mean)
	qui sum logpsa_sd if x==2 & endnote == `endnote' 
	local sd2 = r(mean)
	qui sum ntotal if x==2 & endnote == `endnote' 
	local n2 = r(mean)
	qui sum bmimean if x==2 & endnote == `endnote'
	local bmi2 = r(mean)

	local mean = (`mean1'*`n1'+`mean2'*`n2')/(`n1'+`n2')
	local sd = sqrt(((`n1'-1)*`sd1'^2+(`n2'-1)*`sd2'^2)/(`n1'+`n2'-2))
	local bmi = (`bmi1'*`n1'+`bmi2'*`n2')/(`n1'+`n2')

	qui replace logpsa = `mean' if x==2 & endnote == `endnote' 
	qui replace logpsa_sd = `sd' if x==2 & endnote == `endnote' 
	qui replace ntotal = `n1'+`n2' if x==2 & endnote == `endnote'
	qui replace bmimean = `bmi' if x==2 & endnote == `endnote'
	
	drop if x==1 & endnote == `endnote'
	replace x = x-1 if endnote == `endnote'
}

foreach endnote in 214346 {

	qui sum logpsa if x==1 & endnote == `endnote'
	local mean1 = r(mean)
	qui sum logpsa_sd if x==1 & endnote == `endnote'
	local sd1 = r(mean)
	qui sum ntotal if x==1 & endnote == `endnote'
	local n1 = r(mean)
	qui sum bmimean if x==1 & endnote == `endnote'
	local bmi1 = r(mean)
	qui sum logpsa if x==2 & endnote == `endnote' 
	local mean2 = r(mean)
	qui sum logpsa_sd if x==2 & endnote == `endnote' 
	local sd2 = r(mean)
	qui sum ntotal if x==2 & endnote == `endnote' 
	local n2 = r(mean)
	qui sum bmimean if x==2 & endnote == `endnote'
	local bmi2 = r(mean)
	qui sum logpsa if x==3 & endnote == `endnote' 
	local mean3 = r(mean)
	qui sum logpsa_sd if x==3 & endnote == `endnote' 
	local sd3 = r(mean)
	qui sum ntotal if x==3 & endnote == `endnote' 
	local n3 = r(mean)
	qui sum bmimean if x==3 & endnote == `endnote'
	local bmi3 = r(mean)	

	local mean = (`mean1'*`n1'+`mean2'*`n2'+`mean3'*`n3')/(`n1'+`n2'+`n3')
	local sd = sqrt(((`n1'-1)*`sd1'^2+(`n2'-1)*`sd2'^2+(`n3'-1)*`sd3'^2)/(`n1'+`n2'+`n3'-2))
	local bmi = (`bmi1'*`n1'+`bmi2'*`n2'+`bmi3'*`n3')/(`n1'+`n2'+`n3')

	qui replace logpsa = `mean' if x==3 & endnote == `endnote' 
	qui replace logpsa_sd = `sd' if x==3 & endnote == `endnote' 
	qui replace ntotal = `n1'+`n2' if x==3 & endnote == `endnote'
	qui replace bmimean = `bmi' if x==3 & endnote == `endnote'
	
	drop if x==1 & endnote == `endnote'
	drop if x==2 & endnote == `endnote'
	replace x = x-2 if endnote == `endnote'
}


*Now compare normal:overweight/obese - MD is probably best
gen md_overweight = .
gen md_se_overweight = .
gen md_obese = .
gen md_se_obese = .
gen md_overobese = .
gen md_se_overobese = .
gen n1 = .
gen n2 = .
gen n3 = .
gen logpsa_1 = .
gen logpsa_sd_1 = .
gen logpsa_2 = .
gen logpsa_sd_2 = .
gen logpsa_3 = .
gen logpsa_sd_3 = .
gen bmi_1 = .
gen bmi_2 = .
gen bmi_3 = .

qui sum study
local max = r(max)
forvalues i = 1/`max' {
	qui sum logpsa_sd if study == `i' & bmi_group == 1
	local sd1 = r(mean)
	replace logpsa_sd_1 = `sd1' if study == `i'
	qui sum logpsa_sd if study == `i' & bmi_group == 2
	local sd2 = r(mean)
	replace logpsa_sd_2 = `sd2' if study == `i'
	qui sum logpsa_sd if study == `i' & bmi_group == 3
	local sd3 = r(mean)
	replace logpsa_sd_3 = `sd3' if study == `i'
	
	qui sum logpsa if study == `i' & bmi_group == 1
	local u1 = r(mean)
	replace logpsa_1 = `u1' if study == `i'
	qui sum logpsa if study == `i' & bmi_group == 2
	local u2 = r(mean)
	replace logpsa_2 = `u2' if study == `i'
	qui sum logpsa if study == `i' & bmi_group == 3
	local u3 = r(mean)	
	replace logpsa_3 = `u3' if study == `i'

	qui sum ntotal if study == `i' & bmi_group == 1
	local n1 = r(mean)
	qui sum ntotal if study == `i' & bmi_group == 2
	local n2 = r(mean)
	qui sum ntotal if study == `i' & bmi_group == 3
	local n3 = r(mean)
	
	qui sum bmimean if study == `i' & bmi_group == 1
	local bmi1 = r(mean)
	qui sum bmimean if study == `i' & bmi_group == 2
	local bmi2 = r(mean)
	qui sum bmimean if study == `i' & bmi_group == 3
	local bmi3 = r(mean)
	
	local sdp1 = sqrt(((`n1'-1)*`sd1'^2+(`n2'-1)*`sd2'^2)/(`n1'+`n2'-2))
	local sdp2 = sqrt(((`n1'-1)*`sd1'^2+(`n3'-1)*`sd3'^2)/(`n1'+`n3'-2))
	local sdp3 = sqrt(((`n2'-1)*`sd2'^2+(`n3'-1)*`sd3'^2)/(`n2'+`n3'-2))	
	
	local md1 = (`u2'-`u1') 
	local md2 = (`u3'-`u1')	
	local md3 = (`u3'-`u2')		
	
	replace md_overweight = (`u2'-`u1') if study == `i'
	replace md_obese = (`u3'-`u1') if study == `i'
	replace md_overobese = (`u3'-`u2') if study == `i'	
	
	replace md_se_overweight = sqrt((`n1'+`n2')*`sdp1'^2/(`n1'*`n2')) if study == `i'
	replace md_se_obese = sqrt((`n1'+`n3')*`sdp2'^2/(`n1'*`n3')) if study == `i'
	replace md_se_overobese = sqrt((`n2'+`n3')*`sdp3'^2/(`n2'*`n3')) if study == `i'	
	
	forvalues j = 1/3 {
		replace n`j' = `n`j'' if study == `i'
		replace bmi_`j' = `bmi`j'' if study == `i'
	}
}

*Replace ntotal with actual totals
qui sum study
local max = r(max)

forvalues i = 1/`max' {
	qui sum ntotal if study == `i'
	qui replace ntotal = r(sum) if study == `i'
}

*Make sure Yang's results are deleted (4480)
foreach var of varlist beta* md* n1 n3 logpsa_1 logpsa_sd_1 logpsa_3 logpsa_sd_3 bmi_1 bmi_3 {
	qui replace `var' = . if endnote == 4480
}
qui replace albatross = 1 if endnote == 4480

*Get to 1 record per study
keep if bmi_group == 1

*Put everything together
append using "Original data\AD\PSA categories.dta"

*Drop 36677 & 13747 since they are already in the MD
drop if (endnote == 36677 | endnote == 13747) & controlsn != .

*Deal with the continuous one first (5671)
qui replace beta = adjor if endnote == 5671
qui replace beta_se = adjse if endnote == 5671
qui replace type = "Beta" if endnote == 5671

*And now 12998
qui replace beta = adjor if endnote == 12998
qui replace beta_se = adjse if endnote == 12998
qui replace type = "Beta" if endnote == 12998

*So the ORs will just be used in the AP, therefore I'm going to select PSA>4 as the outcome of interest (comparable between the two studies, and useful point).
*Delete data that can't be used
drop if outcome == "PSA >= 2.5" | outcome == "PSA >= 10.0" | outcome == "PSA>2.5"

*Replace ntotals
qui sum controlsn if endnote == 8796 & (bmisubgroup == "Obese" | bmisubgroup == "Normal")
qui replace ntotal = r(sum) if endnote == 8796
qui replace ntotal = casesn if endnote == 5671 
qui replace ntotal = controlsn if endnote == 12998 | endnote == 36480 

*First keep only the obese for 8796 P values, as this will be the comparison (obese versus normal, rather than P value for beta)
drop if endnote==8796 & bmisubgroup!="Obese"
replace type = "OR" if endnote == 8796

*Do things with the new additions (200000+)
replace beta = adjor if endnote == 208273
replace beta_se = adjse if endnote == 208273

replace beta = adjor if endnote == 209405
replace beta_se = adjse if endnote == 209405

replace beta = adjor/5 if endnote == 212003
replace beta_se = adjse/5 if endnote == 212003
replace type = "Continuous" if endnote == 208273 | endnote == 209405 | endnote == 212003

qui sum adjor if endnote == 209405 & bmirange == "25-30"
replace md_overweight = r(mean) if endnote == 209405
qui sum adjse if endnote == 209405 & bmirange == "25-30"
replace md_se_overweight = r(mean) if endnote == 209405
qui sum controlsn if endnote == 209405 & bmirange == "<25"
replace n1 = r(mean) if endnote == 209405
qui sum controlsn if endnote == 209405 & bmirange == "25-30"
replace n2 = r(mean) if endnote == 209405
qui sum controlsn if endnote == 209405 & bmirange == "30-35"
replace n3 = r(mean) if endnote == 209405
qui sum adjor if endnote == 209405 & bmirange == "30-35"
replace md_obese = r(mean) if endnote == 209405
qui sum adjse if endnote == 209405 & bmirange == "30-35"
replace md_se_obese = r(mean) if endnote == 209405

*Note, I did this manually. Need to write a program to do this really.
replace bmi_1 = 21.75 if endnote == 209405
replace bmi_2 = 27.46 if endnote == 209405
replace bmi_3 = 32.22 if endnote == 209405
replace study_bmi_mean = 26.04 if endnote == 209405
replace study_bmi_sd = 4.26 if endnote == 209405

*Drop parts of 209405 that can't be used
drop if endnote == 209405 & bmirange != "1 kg/m2"

drop if endnote == 208801 //ProtecT study, use IPD

replace ntotal = controlsn if ntotal == .
replace ntotal = casesn if ntotal == .

replace beta = adjor if albatross == 1 & beta == .

replace ethnicity = group if endnote == 208273

*Replace P values
replace pvalue = adjpvalue if pvalue ==.
drop othervariables x bmi_lower bmi_upper study notes bmisd-yearsfu bmisubgroup

*qui replace the authors, if there is more than 1, with just the primary author
qui egen author_1 = ends(author), punct(" and ") trim head
qui egen author_primary = ends(author_1), punct(,) trim head
qui egen author_primary_2 = ends(author_primary) if strpos(author_primary,".")>0, trim tail
qui replace author_primary = author_primary_2 if author_primary_2 !=""
qui egen author_primary_3 = ends(author_primary) if strpos(author_primary_2,".")>0, trim tail
qui replace author_primary = author_primary_3 if author_primary_3 !=""
qui egen author_primary_4 = ends(author_primary) if strpos(author_primary_3,".")>0, trim tail
qui replace author_primary = author_primary_4 if author_primary_4 !=""
drop author_primary_2 author_primary_3 author_primary_4 author_1
order author_primary, b(author)
rename author author_all
rename author_primary author
order author_all, last

*Culp renaming for ethnicity
replace author = "Culp (Whites)" if author == "Culp" & ethnicity == "Whites"
replace author = "Culp (African American)" if author == "Culp" & ethnicity == "AA"
replace author = "Culp (Hispanic)" if author == "Culp" & ethnicity == "Hispanic"

*And Banez for studies
replace author = "Banez (SEARCH)" if author == "Banez" & studyname == "SEARCH"
replace author = "Banez (Duke)" if author == "Banez" & studyname == "Duke"
replace author = "Banez (Johns Hopkins)" if author == "Banez" & studyname == "Johns Hopkins"

*And Waters for ethnicity
gen x1 = " "
gen x2 = "("
gen x3 = ")"
egen author2 = concat(author x1 x2 ethnicity x3)
replace author = author2 if endnote == 208273
drop x1-author2

*NHANES has 3 studies, 209795 (use for cat), 212003 (use for cont), and 12998 (now removed at database level)
replace beta = . if endnote == 209795
replace beta_se = . if endnote == 209795

*Standardised beta = %change in PSA per 5 kg/m2 increase in BMI
gen beta5 = (beta*5)
gen beta5_se = (beta_se*5)

*Drop the PLCO study
drop if endnote == 14033
sort albatross year author

rename ntotal N
qui replace pvalue = 2*normal(-abs(beta/beta_se)) if pvalue == .

qui gen p=pvalue
qui gen e=beta
qui gen n=N 

qui gen included = "Included" if albatross == .
qui replace included = "Excluded" if included == ""

qui replace e = 1 if endnote == 8796
qui replace e = -1 if endnote == 36480 | endnote == 4480 | endnote == 214367 | endnote == 208597 | endnote == 5671 | endnote == 210783

drop psaassay bmiascer screening psa_mean-bmimean logpsa-logpsa_se bmi_group

qui gen ethnic_bin = 1
qui replace ethnic_bin = 2 if ethnicity == "Chinese" | ethnicity == "Japan" | ethnicity == "Japanese" | ethnicity == "Japanese (100)" | ethnicity == "Korean" | ethnicity == "Black" | ethnicity=="Multiethnic" | author == "Whittemore (Blacks)" | ethnicity == "Nigerian" | ethnicity == "African American" | ethnicity == "Japanese American" | ethnicity == "Native Hawaiian"
capture label define ethnic_bin 1 "Caucasian" 2 "Non-Caucasian" 
label values ethnic_bin ethnic_bin
encode overallrisk, gen(rob)
order rob, a(overallrisk)

save "PSA AD.dta", replace

use "PSA AD.dta", clear

keep endnote study_bmi_mean study_bmi_sd N md* beta5 beta5_se p e albatross type

save "PSA AD Table.dta", replace
}

*IPD results

{
use "IPD PSA - complete.dta", clear
merge 1:1 _n using "IPD PSA - complete (categorical).dta"
drop _merge
rename pvalue p
drop ncases ncontrols bmi_cases* bmi_controls* outcome
rename ntotal N
gen n = N
rename bmi_all study_bmi_mean
rename bmi_all_sd study_bmi_sd
rename overallriskofbias overallriskofbiaspsa
gen rob = 2
rename whichvariablesare variablesadjustedforpsa
gen e = beta5
gen included = "Included"
gen ethnic_bin = 1
gen author = studyname

foreach var in obese overweight {
	rename `var'_psa md_`var'
	rename `var'_psa_se md_se_`var'
	drop `var'_psa_p 
}

rename normal_ntotal n1
rename overweight_ntotal n2
rename obese_ntotal n3

rename normal_bmi bmi_1
rename overweight_bmi bmi_2
rename obese_bmi bmi_3

save "IPD PSA - to add", replace
}

*Add in IPD to AD
use "PSA AD.dta", clear
append using "IPD PSA - to add"
gen type2 = "IPD" if type == "IPD"
replace type2 = "AD" if type != "IPD"

save "PSA pre-graphs.dta", replace

use "PSA pre-graphs.dta", clear
local obs = c(N)
local obs5 = `obs'+1
set obs `obs5'
replace author = "Summary" in `obs5'
replace studyname = "Summary" in `obs5'

foreach var in overweight obese {
	gen md_p_`var' = 2*normal(-abs(md_`var'/md_se_`var')), a(md_se_`var')
	gen md_`var'_fixed = ., a(md_p_`var')
	gen md_se_`var'_fixed = ., a(md_`var'_fixed)
	gen md_p_`var'_fixed = ., a(md_se_`var'_fixed)
}

gen beta5_p = p, b(p)
gen beta5_fixed = ., a(beta5_p)
gen beta5_se_fixed = ., a(beta5_fixed)
gen beta5_p_fixed = ., a(beta5_se_fixed)

rename year year_num
tostring year_num, gen(year)
replace year = "IPD" if type == "IPD"

*All
metan beta5 beta5_se, eform xlabel(0.8, 0.85, 0.9, 0.95, 1.05) xtitle("Percentage change in PSA for a 5 kg/m{superscript:2} increase in BMI",size(vsmall)) rcols(N) label(namevar=author, yearvar=year) by(type2) force textsize(90) random second(fixed)
graph save "Graphs\BMI-PSA\PSA.gph", replace
metan beta5 beta5_se, eform random second(fixed) nograph
qui replace beta5 = ln(r(ES)) in `obs5'
qui replace beta5_se = r(selogES) in `obs5'
qui replace beta5_p = 2*normal(-r(z)) in `obs5'
qui replace beta5_fixed = ln(r(ES_2)) in `obs5'
qui replace beta5_se_fixed = (r(selogES_2)) in `obs5'
qui replace beta5_p_fixed = 2*normal(-abs(ln(r(ES_2))/r(selogES_2))) in `obs5'

metan md_overweight md_se_overweight, eform xlabel(0.6, 0.7, 0.8, 0.9, 1.1, 1.2) xtitle("Percentage change in PSA, overweight versus normal weight men",size(vsmall)) rcols(n1 n2) label(namevar=author, yearvar=year) by(type2) force textsize(90) random second(fixed)
graph save "Graphs\BMI-PSA\Overweight PSA.gph", replace
metan md_obese md_se_obese, eform xlabel(0.6, 0.7, 0.8, 0.9, 1.1, 1.2) xtitle("Percentage change in PSA, obese versus normal weight men",size(vsmall)) rcols(n1 n3) label(namevar=author, yearvar=year) by(type2) force textsize(90) random second(fixed)
graph save "Graphs\BMI-PSA\Obese PSA.gph", replace

foreach var in overweight obese {
	metan md_`var' md_se_`var' if type != "", eform nograph random second(fixed)
	qui replace md_`var' = ln(r(ES)) in `obs5'
	qui replace md_se_`var' = r(selogES) in `obs5'
	qui replace md_p_`var' = 2*normal(-r(z)) in `obs5'
	qui replace md_`var'_fixed = ln(r(ES_2)) in `obs5'
	qui replace md_se_`var'_fixed = (r(selogES_2)) in `obs5'
	qui replace md_p_`var'_fixed = 2*normal(-abs(ln(r(ES_2))/r(selogES_2))) in `obs5'
}

metafunnel beta5 beta5_se if type != "", eform xtitle("Percentage change in PSA") ytitle("SE of beta coefficient")
graph save "Graphs\BMI-PSA\PSA Funnel.gph", replace

qui sum N if type != ""
qui replace N = r(sum) in `obs5'

*Need a weighted estimate of the mean BMI in each group
forvalues i = 1/3 {
	
	gen x`i' = bmi_`i'*n`i' if type != "" & albatross != 1
	qui sum x`i' 
	local xsum`i' = r(sum)
	qui sum n`i'
	local nsum`i' = r(sum)
	local mean = round(`xsum`i''/`nsum`i'',0.001)
	qui replace bmi_`i' = `mean' in `obs5'
	qui sum n`i' if type != ""
	qui replace n`i' = r(sum) in `obs5'
	drop x`i'
	
	gen x`i' = logpsa_`i'*n`i' if type != "" & albatross != 1
	qui sum x`i'
	local xsum`i' = r(sum)
	qui sum n`i'
	local nsum`i' = r(sum)
	local mean = `xsum`i''/`nsum`i''
	replace logpsa_`i' = `mean' in `obs5'
	drop x`i'
	
}

gen x = study_bmi_mean*N if type != "" & albatross != 1
qui sum x 
local xsum = r(sum)
qui sum N if type != "" & albatross != 1
local nsum = r(sum)
local mean = `xsum'/`nsum'
replace study_bmi_mean = `mean' in `obs5'
drop x

foreach var of varlist study_bmi_sd logpsa_sd_* {
	gen x = `var'^2*(N-1) if type != "" & albatross != 1
	qui sum x
	local xsum = r(sum)
	qui sum N if type != "" & albatross != 1
	local nsum = r(sum)-12
	local sd = sqrt(`xsum'/`nsum')
	replace `var' = `sd' in `obs5'
	drop x
}

gen Included = included
replace Included = "Included - IPD" if type == "IPD"
replace Included = "Included - AD" if Included == "Included"

*Add in Loeb - 31/01/2019
local obs = c(N) +1
set obs `obs'
replace endnoteid = 213859 in `obs'
replace albatross = 1 in `obs'
replace overall = "Medium" in `obs'
replace rob = 2 in `obs'
replace variablesadjustedforpsa = "Age" in `obs'
replace author = "Loeb" in `obs'
replace year = 2009 in `obs'
replace studyname = "BLSA" in `obs'
replace studylocation = "USA" in `obs'
replace midyear = 1958 in `obs'
replace ethnicity = "Caucasian" in `obs'
replace type = "Continuous" in `obs'
replace e = -1 in `obs'
replace n = 994 in `obs'
replace p = 0.06 in `obs'
replace included = "Excluded" in `obs'
replace Included = "Excluded" in `obs'
replace ethnic_bin = 1 in `obs'

albatross n p e, type(beta) by(Included) color con(0.05 0.1 0.15)
graph save "Graphs\BMI-PSA\PSA Albatross.gph", replace
graph export "Graphs\Final graphs\For processing\PSA albatross.tif", replace as(tif) width(1200)

save "PSA complete.dta", replace

*Meta-regression
use "PSA complete.dta", clear
replace midyear = 2004 if author == "ProtecT"
replace midyear = 1998 if author == "PLCO"
replace midyear = 1995 if author == "Krimpen"

metareg beta5 ethnic_bin rob midyear if albatross != 1, wsse(beta5_se)
metareg beta5 study_bmi_mean if albatross != 1, wsse(beta5_se)

********************************************************************************

*Linear Table
use "PSA complete.dta", clear
keep endnote author year studyname studylocation ethnicity midyear beta5 beta5_se N p type albatross e
order endnote author year studyname studylocation ethnicity midyear beta5 N p type
drop if author == "Summary"
replace albatross = 0 if albatross == .
replace endnote = -1 if author == "ProtecT"
replace endnote = -2 if author == "PLCO"
replace endnote = -3 if author == "Krimpen"
sort albatross year author
gen or = round(100*(exp(beta5)-1),0.01)
gen l = round((100*(exp(beta5-1.96*beta5_se)-1)),0.01)
gen u = round((100*(exp(beta5+1.96*beta5_se)-1)),0.01)
gen x = "("
gen y = ")"
gen z = " "
gen a = " to "
gen b = "%"
egen ci = concat(or b z x l b a u b y) if beta5 != .
order ci, a(beta5)

replace ci = subinstr(ci,"-.","-0.",.)
replace ci = subinstr(ci," ."," 0.",.)

forvalues i = 1/9 {
	replace ci = subinstr(ci,".`i'%",".`i'0%",.)
}

qui replace ci = "Positive" if e > 0 & albatross == 1
qui replace ci = "Negative" if e < 0 & albatross == 1

*Not convinced by the P value in 19431 (0.0001 when 95% CI crosses 0)
replace p = 2*normal(-abs(beta5/beta5_se)) if endnote == 19431

drop beta5 beta5_se l-a or e

qui replace ethnicity = "" if ethnicity == studylocation | ethnicity == "Australian" | ethnicity == "Canadian" | ethnicity == "Chinese" | ethnicity == "Czech" | ethnicity == "Finnish" | ///
ethnicity == "Greece" | ethnicity == "Greek" | ethnicity == "Italian" | ethnicity == "Japanese" | ethnicity == "Korean" | ethnicity == "Norwegian" | ethnicity == "Swedish" | ethnicity == "German"

replace type = "" if albatross == 1
*209795 was just for categories, not main MA
drop if endnote == 209795

replace endnote = 208273 if author == "Waters (African American)"
replace endnote = 208274 if author == "Waters (European American)"
replace endnote = 208275 if author == "Waters (Japanese American)"
replace endnote = 208276 if author == "Waters (Latino)"
replace endnote = 208277 if author == "Waters (Native Hawaiian)"

save "Tables\PSA - Supplementary table.dta", replace 

********************************************************************************

*Catgeorical table
use "PSA complete.dta", clear
replace logpsa_1 = normal_logpsa if normal_logpsa != .
replace logpsa_sd_1 = normal_logpsa_sd if normal_logpsa_sd != .
replace logpsa_2 = overweight_logpsa if overweight_logpsa != .
replace logpsa_sd_2 = overweight_logpsa_sd if overweight_logpsa_sd != .
replace logpsa_3 = obese_logpsa if obese_logpsa != .
replace logpsa_sd_3 = obese_logpsa_sd if obese_logpsa_sd != .
keep endnote studyname author md_overweight md_se_overweight md_obese md_se_obese n1 n2 n3 logpsa_1 logpsa_sd_1 logpsa_2 logpsa_sd_2 logpsa_3 logpsa_sd_3 bmi_1 bmi_2 bmi_3 year type
keep if bmi_1 != . | logpsa_1 != .
drop if author == "Summary"
sort type year author

gen s = " "
gen b1 = "("
gen b2 = ")"
gen to = " to "
gen p = "%"
foreach var of varlist logpsa_1-logpsa_sd_3 {
	replace `var' = round(`var',0.01)
}
foreach var of varlist bmi_1-bmi_3 {
	replace `var' = round(`var',0.1)
}
foreach var of varlist n1-n3 {
	replace `var' = round(`var',1)
}
egen psa_1 = concat(logpsa_1 s b1 logpsa_sd_1 b2)
egen psa_2 = concat(logpsa_2 s b1 logpsa_sd_2 b2)
egen psa_3 = concat(logpsa_3 s b1 logpsa_sd_3 b2)

order year type endnote studyname author bmi_1 psa_1 n1 bmi_2 psa_2 n2 bmi_3 psa_3 n3 

foreach i in 1 2 3 {
	replace psa_`i' = subinstr(psa_`i',".","0.",.) if logpsa_`i' > -1 & logpsa_`i' < 1 & logpsa_sd_`i' < 1
	replace psa_`i' = subinstr(psa_`i',"(.","(0.",.)
	replace psa_`i' = subinstr(psa_`i',"-.","-0.",.)
	
	forvalues j = 1/9 {
		replace psa_`i' = subinstr(psa_`i',".`j')",".`j'0)",.)
		replace psa_`i' = subinstr(psa_`i',".`j' ",".`j'0 ",.)
	}
}

gen md1 = round(100*(exp(md_overweight)-1),0.01)
gen ci11 = round(100*(exp(md_overweight-1.96*md_se_overweight)-1),0.01)
gen ci12 = round(100*(exp(md_overweight+1.96*md_se_overweight)-1),0.01)
gen md2 = round(100*(exp(md_obese)-1),0.01)
gen ci21 = round(100*(exp(md_obese-1.96*md_se_obese)-1),0.01)
gen ci22 = round(100*(exp(md_obese+1.96*md_se_obese)-1),0.01)

egen meandif1 = concat(md1 p s b1 ci11 p to ci12 p b2)
egen meandif2 = concat(md2 p s b1 ci21 p to ci22 p b2)

foreach i in 1 2 {
	forvalues j = 1/9 {
		replace meandif`i' = subinstr(meandif`i',".`j'%",".`j'0%",.)
		replace meandif`i' = subinstr(meandif`i',"-`j'%","-`j'.00%",.)
		replace meandif`i' = subinstr(meandif`i',"-20%","-20.00%",.)
	}
	replace meandif`i' = subinstr(meandif`i',".","0.",1) if md`i' > -1 & md`i' < 1
}

save "Tables\PSA - Supplementary table (categorical).dta", replace 




********************************************************************************

*PSA analyses complete

********************************************************************************

}


********************************************************************************
********************************************************************************

*Prostate Cancer analyses
{

*MD
{
set more off
clear

********************************************************************************
**Prostate cancer mean differences**
********************************************************************************

use "Original data\AD\PCa MD.dta", clear

*Type = MD
gen type = "MD"

*qui replace the authors, if there is more than 1, with just the primary author
qui egen author_1 = ends(author), punct(" and ") trim head
qui egen author_primary = ends(author_1), punct(,) trim head
qui egen author_primary_2 = ends(author_primary) if strpos(author_primary,".")>0, trim tail
qui replace author_primary = author_primary_2 if author_primary_2 !=""
qui egen author_primary_3 = ends(author_primary) if strpos(author_primary_2,".")>0, trim tail
qui replace author_primary = author_primary_3 if author_primary_3 !=""
qui egen author_primary_4 = ends(author_primary) if strpos(author_primary_3,".")>0, trim tail
qui replace author_primary = author_primary_4 if author_primary_4 !=""
drop author_primary_2 author_primary_3 author_primary_4 author_1
order author_primary, b(author)
rename author author_all
rename author_primary author
order author_all, last

*Mean differences require an SMD, therefore the pooled SD must be calculated
qui gen smd = .
qui gen smd_se = .
qui gen study_bmi_mean = .

*Total number of participants
qui gen ntotal = ncases+ncontrols

*Three studies reported SE, not SD, so replacing this
qui replace sdcases = sdcases*sqrt(ncases) if endnoteid == 24897 | endnoteid == 44138 | endnote == 207948
qui replace sdcontrols = sdcontrols*sqrt(ncontrols) if endnoteid == 24897 | endnoteid == 44138 | endnote == 207948

*Calculate SD from 95% CI
qui egen lcicase = ends(cicases), punct("-") trim head
destring lcicase, replace
qui egen ucicase = ends(cicases), punct("-") trim tail
destring ucicase, replace
qui replace sdcases = sqrt(ncases)*(ucicase-lcicase)/3.92 if sdcases ==.

qui egen lcicontrol = ends(cicontrols), punct("-") trim head
destring lcicontrol, replace
qui egen ucicontrol = ends(cicontrols), punct("-") trim tail
destring ucicontrol, replace
qui replace sdcontrols = sqrt(ncontrols)*(ucicontrol-lcicontrol)/3.92 if sdcontrols==.

*Calculate SD from IQR
qui egen liqrcase = ends(iqrcases), punct("-") trim head
destring liqrcase, replace
qui egen uiqrcase = ends(iqrcases), punct("-") trim tail
destring uiqrcase, replace
qui replace sdcases = (uiqrcase-liqrcase)/1.35 if sdcases ==.

qui egen liqrcontrol = ends(iqrcontrols), punct("-") trim head
destring liqrcontrol, replace
qui egen uiqrcontrol = ends(iqrcontrols), punct("-") trim tail
destring uiqrcontrol, replace
qui replace sdcontrols = (uiqrcontrol-liqrcontrol)/1.35 if sdcontrols==.

*Kopp (5189) has 5-95% percentiles, not 25-75%
qui replace sdcontrols = (uiqrcontrol-liqrcontrol)/3.29 if endnote == 5189
qui replace sdcases = (uiqrcase-liqrcase)/3.29 if endnote == 5189

*Also 209072
qui replace sdcontrols = (uiqrcontrol-liqrcontrol)/3.29 if endnote == 209072
qui replace sdcases = (uiqrcase-liqrcase)/3.29 if endnote == 209072
replace studyname = "danish diet, cancer, and health (dch)" if endnote == 209072

*Calculate SD_pooled
qui gen sd_pooled = sqrt(((ncases-1)*sdcases^2+(ncontrols-1)*sdcontrols^2)/(ntotal-2))
qui replace sd_pooled = bmisd if bmisd !=.

*SD_pooled from P value alone (1 study)
*Assume same SD in cases and controls
qui replace sd_pooled = sqrt(abs((meanmediancases-meanmediancontrols)^2*ncases*ncontrols/(ntotal*invnormal(pvalue/2)^2))) if pvalue!=. & sd_pooled ==.

*Drop variables that can't be used (now used)
drop lcicase ucicase lcicontrol ucicontrol liqrcase uiqrcase liqrcontrol uiqrcontrol

*Calculate SD_pooled from binary data (1 study)
sort endnoteid x
qui egen study = group(endnoteid outcome group)
bysort study: egen count = count(study)
qui replace count = . if count == 1
*Assume if count == 2 that the data has been made binary
qui sum study
local max = r(max)
forvalues i = 1 (1) `max' {
	qui sum count if study == `i'
	if r(mean) ==2 {
		qui sum ncases if study == `i' & x == 1
		local a = r(mean)
		qui sum ncontrols if study == `i' & x == 1
		local b = r(mean)
		qui sum ncases if study == `i' & x == 2
		local c = r(mean)
		qui sum ncontrols if study == `i' & x == 2
		local d = r(mean)
		local or = `a'*`d'/(`b'*`c')
		local pi1 = `a'/(`a'+`b')
		local pi2 = `c'/(`c'+`d')
		qui replace smd = invnormal(`pi1')-invnormal(`pi2') if study == `i'
		drop if study == `i' & x == 2
		qui replace ncases = `a'+`c' if study == `i'
		qui replace ncontrols = `b'+`d' if study == `i'
		qui replace ntotal = ncases+ncontrols if study == `i'
		qui replace smd_se = sqrt(1/(`a'+`c')+1/(`b'+`d')+(invnormal(`pi1')-invnormal(`pi2'))^2/(2*(`a'+`b'+`c'+`d'))) if study == `i'
		dis sqrt((`pi1'*(1-`pi1'))^-2*(`pi1'*(1-`pi1')/(`a'+`c'))+(`pi2'*(1-`pi2'))^-2*(`pi2'*(1-`pi2')/(`b'+`d')))
	}
}

*Calculate SD_pooled from quantile data (2 studies)
*Note: this is for studies that had to be collapsed as the quantiles broke the age-matching, and required estimates of the SD
*Note: this is tricky and from Chene and Thompson
*Easiest way is to probably cycle through studies and exit if not quantile

*Destring the quantiles
qui egen bmi_lower = ends(bmisubgroup), punct("-") trim head
qui egen bmi_lower2 = ends(bmisubgroup), punct(">") trim tail
replace  bmi_lower = bmi_lower2 if bmi_lower2 !=""
replace bmi_lower = "" if strpos(bmi_lower, "<")>0
destring bmi_lower, replace
qui egen bmi_upper = ends(bmisubgroup), punct("-") trim tail
qui egen bmi_upper2 = ends(bmisubgroup), punct("<") trim tail
replace  bmi_upper = bmi_upper2 if bmi_upper2 !=""
destring bmi_upper, replace

drop bmi_upper2 bmi_lower2

qui gen Xj = .
qui replace Xj = bmi_upper

qui gen Nj = ntotal
qui gen Pj = .
qui gen Zj = .
qui gen Wj = .
qui gen Mj = .
qui gen Mj_cases = .
qui gen Mj_controls = .

sort study x
forvalues study = 1/`max' {
	qui sum count if study == `study'
	if r(mean) >2 & r(mean) <. {
		qui sum x if study == `study'
		local max_x = r(max)
		
		*Total number of participants in the study
		qui sum Nj if study==`study'
		local sum = r(sum)
		
		forvalues value = 1/`max_x' {
			qui sum Nj if x <= `value' & study==`study'
			qui replace Pj = r(sum)/`sum' if x==`value' & study==`study'
		}

		qui replace Zj = invnormal(Pj)
		qui replace Wj = normalden(Zj,0,1)^2/(Pj*(1-Pj))
	
		qui sum x if study==`study'
		local quantile = r(max)-1
		qui sum Wj if study==`study'
		qui replace Wj = Wj*(`quantile'/r(sum))

		*Regress to get the mean and SD of the study, using the beta of Zj for the mean and the constant for the SD
		qui regress Xj Zj if study==`study'
		local m = _b[_cons]
		qui replace study_bmi_mean = _b[_cons] if study == `study'
		local s = _b[Zj]
		qui replace sd_pooled = _b[Zj] if study == `study'
		
		forvalues value = 1/`max_x' {
			if `value' == `max_x' {
				local zj1 = 9999999
			}
			else {
				qui sum Zj if x == `value' & study == `study'
				local zj1 = r(mean)
			}
			if `value' == 1 {
				local zj0 = -999999
			}
			else {
				local value0 = `value'-1
				qui sum Zj if x == `value0' & study == `study'
				local zj0 = r(mean)
			}
			replace Mj = `m'+`s'*((normalden(`zj0')-normalden(`zj1'))/(normal(`zj1')-normal(`zj0'))) if study == `study' & x == `value'
		}
	
		*Now need to get the means for cases and controls seperately
		foreach var in cases controls {
			qui replace Nj = n`var'
			qui replace Pj = .
			qui replace Zj = .
			qui replace Wj = .
			qui replace Mj = .
		
			*Total number of participants in the study
			qui sum Nj if study==`study'
			local sum = r(sum)
		
			forvalues value = 1/`max_x' {
				qui sum Nj if x <= `value' & study==`study'
				qui replace Pj = r(sum)/`sum' if x==`value' & study==`study'
			}
		
			qui replace Zj = invnormal(Pj)
			qui replace Wj = normalden(Zj,0,1)^2/(Pj*(1-Pj))
	
			qui sum x if study==`study'
			local quantile = r(max)-1
			qui sum Wj if study==`study'
			qui replace Wj = Wj*(`quantile'/r(sum)) if study == `study'
	
			*Regress to get the mean and SD of the study, using the beta of Zj for the mean and the constant for the SD
			qui regress Xj Zj [aweight = Wj] if study==`study' 
			local m = _b[_cons]
			qui replace meanmedian`var' = _b[_cons] if study == `study'
			local s = _b[Zj]
			qui replace sd`var' = _b[Zj] if study == `study'
		
			forvalues value = 1/`max_x' {
				if `value' == `max_x' {
					local zj1 = 9999999
				}
				else {
					qui sum Zj if x == `value' & study == `study'
					local zj1 = r(mean)
				}
				if `value' == 1 {
					local zj0 = -999999
				}
				else {
					local value0 = `value'-1
					qui sum Zj if x == `value0' & study == `study'
					local zj0 = r(mean)
				}
				replace Mj_`var' = `m'+`s'*((normalden(`zj0')-normalden(`zj1'))/(normal(`zj1')-normal(`zj0'))) if study == `study' & x == `value'
			}
	
		}
	*Replace ncases/controls in x==1 with total cases/controls
	qui sum ncases if study==`study'
	qui replace ncases = r(sum) if study == `study'
	qui sum ncontrols if study==`study'
	qui replace ncontrols = r(sum) if study == `study'	
	qui replace ntotal = ncases+ncontrols if study == `study'
	drop if x>1 & study == `study'
	}
}
qui drop bmi_lower-Wj count Mj-Mj_controls

*Generate SMD
qui replace smd = (meanmediancases-meanmediancontrols)/sd_pooled if smd == .
qui replace smd_se = sqrt(1/ncases+1/ncontrols+smd^2/(2*ntotal)) if smd_se == .
qui replace study_bmi_mean = (meanmediancases*ncases+meanmediancontrols*ncontrols)/ntotal if study_bmi_mean == .

*Individual corrections
replace studyname = "northern sweden health and disease cohort (nshdc)" if endnote == 44138 | endnote == 208273
drop if endnote == 15307 //PLCO study
qui replace studyname = "health professionals follow-up study (hpfs)" if endnote == 3101
qui replace studyname= "malmo diet and cancer (mdc)" if endnote == 6846
replace ethnicity = group if endnote == 24897

*If SMD=ln(OR) per SD, then divide by SD and multiply by 5 to get ln(OR) per 5 kg/m2
qui gen smd5 = smd*5/sd_pooled
qui gen smd_se5 = smd_se*5/sd_pooled
qui replace pvalue = round(2*normal(-abs(smd5/smd_se5)),0.01) if pvalue == .

*Get rid of duplicate studies (with the highest SEs)
qui replace studyname = lower(studyname)
sort studyname group outcome
qui egen study2 = group(studyname outcome group)
bysort study2: egen count = count(study2)
qui replace study2 = . if count == 1
qui replace count = . if count == 1 | count == 0
sort count studyname outcome group 
qui gen include = .
qui replace include = 1 if count == .
bysort study2: egen rank = rank(smd_se5) if study2 !=., unique
qui replace include = 1 if count !=. & rank == 1

*Also, need to keep those studies with the closest BMI measurement to diagnosis (if multiple measurements exist)
qui replace include = . if endnoteid == 42897 | endnoteid == 38107

*Keep only the useful variables
keep endnote albatross author year studyname studylocation overallrisk midyear ethnicity timeof meanyears meanage whichvariables* group outcome ncases ncontrols ///
meanmedian* sdcases sdcontrol pvalue smd smd_se study_bmi ntotal sd_pooled smd5 smd_se5 include type

*Retitle Whittemore for each ethnicity
qui replace author = "Whittemore (Blacks)" if author == "Whittemore" & group == "Blacks"
qui replace author = "Whittemore (Japanese-American)" if author == "Whittemore" & group == "Japanese-American"
qui replace author = "Whittemore (Whites)" if author == "Whittemore" & group == "Whites"
qui replace author = "Whittemore (Chinese-American)" if author == "Whittemore" & group == "Chinese-American"

*More tidying
drop if smd5 == . & albatross == .
drop if outcome == "Gleason 6" | outcome == "Gleason 7" | outcome == "Gleason 8-10"

keep if include == 1
sort endnote

rename meanmediancases bmi_cases
rename meanmediancontrols bmi_controls
rename sdcases bmi_cases_sd
rename sdcontrols bmi_controls_sd
rename study_bmi_mean bmi_all
rename sd_pooled bmi_all_sd
rename smd5 md5
rename smd_se5 md5_se

gen or_hr = "OR"

order endnote albatross outcome md5 md5_se pvalue ncases ncontrols ntotal bmi_cases bmi_cases_sd bmi_controls bmi_controls_sd bmi_all bmi_all_sd studyname

save "PCa MD complete.dta", replace

}

*Categorical
{
************************************************************************************************************************************************************************************************************************************************
**Prostate cancer categories**
************************************************************************************************************************************************************************************************************************************************

use "Original data\AD\PCa categories.dta", clear

*Rename casesn
rename casesn ncases
rename controlsn ncontrols

*qui replace the authors, if there is more than 1, with just the primary author
qui egen author_1 = ends(author), punct(" and ") trim head
qui egen author_primary = ends(author_1), punct(,) trim head
qui egen author_primary_2 = ends(author_primary) if strpos(author_primary,".")>0, trim tail
qui replace author_primary = author_primary_2 if author_primary_2 !=""
qui egen author_primary_3 = ends(author_primary) if strpos(author_primary_2,".")>0, trim tail
qui replace author_primary = author_primary_3 if author_primary_3 !=""
qui egen author_primary_4 = ends(author_primary) if strpos(author_primary_3,".")>0, trim tail
qui replace author_primary = author_primary_4 if author_primary_4 !=""
drop author_primary_2 author_primary_3 author_primary_4 author_1
order author_primary, b(author)
rename author author_all
rename author_primary author
order author_all, last

*Destring the quantiles/categories
qui egen bmi_lower = ends(bmirange), punct("-") trim head
qui egen bmi_lower2 = ends(bmirange), punct(">") trim tail
replace  bmi_lower = bmi_lower2 if bmi_lower2 !=""
replace bmi_lower = "" if strpos(bmi_lower, "<")>0
destring bmi_lower, replace
qui egen bmi_upper = ends(bmirange), punct("-") trim tail
qui egen bmi_upper2 = ends(bmirange), punct("<") trim tail
replace  bmi_upper = bmi_upper2 if bmi_upper2 !=""
destring bmi_upper, replace

drop bmi_upper2 bmi_lower2

*Take out AP, they'll only get in the way
preserve
keep if albatross == 1
save "Albatross - PCa.dta", replace
restore
drop if albatross == 1

*Also need to remove continuous data

qui gen continuous = 1 if strpos(bmirange, "unit") >0 | strpos(bmirange, "SD") >0 | strpos(bmisubgroup, "unit") >0 | strpos(bmisubgroup, "kg/m2") >0 
preserve
drop if continuous == .
gen type = "continuous"
save "PCa continuous.dta", replace
restore
drop if continuous == 1
drop continuous
gen type = "categorical"

*Also need to condense 7584, 8774, 39497, 44677, 19122, 209177 (but remove the "trend" in bmirange first)
*These will be treated like SMDs once the mean/SD for cases/controls has been calculated
drop if bmirange == "trend" & endnote == 44677
qui gen condense = 1 if strpos(notes, "condense")>0
qui replace condense = 1 if strpos(notes, "Condense")>0

preserve

keep if condense == 1
qui gen Xj = .
qui replace Xj = bmi_upper

sort endnoteid x
qui egen study = group(endnoteid outcome group)
bysort study: egen count = count(study)
qui replace count = . if count == 1

qui gen meanmediancases = .
qui gen meanmediancontrols = .
qui gen sdcases = .
qui gen sdcontrols = .
qui gen sd_pooled = .
qui gen ntotal = .
qui gen Nj = .
qui gen Pj = .
qui gen Zj = .
qui gen Wj = .
qui gen Mj = .
qui gen Mj_cases = .
qui gen Mj_controls = .

qui sum study
local max = r(max)
sort study x
forvalues study = 1/`max' {
	*Now need to get the means for cases and controls seperately
	qui sum x
	local max_x = r(max)
	foreach var in cases controls {
		qui replace Nj = n`var'
		qui replace Pj = .
		qui replace Zj = .
		qui replace Wj = .
		qui replace Mj = .
	
		*Total number of participants in the study
		qui sum Nj if study==`study'
		local sum = r(sum)
	
		forvalues value = 1/`max_x' {
			qui sum Nj if x <= `value' & study==`study'
			qui replace Pj = r(sum)/`sum' if x==`value' & study==`study'
		}
	
		qui replace Zj = invnormal(Pj)
		qui replace Wj = normalden(Zj,0,1)^2/(Pj*(1-Pj))

		qui sum x if study==`study'
		local quantile = r(max)-1
		qui sum Wj if study==`study'
		qui replace Wj = Wj*(`quantile'/r(sum)) if study == `study'

		*Regress to get the mean and SD of the study, using the beta of Zj for the mean and the constant for the SD
		qui regress Xj Zj [aweight = Wj] if study==`study' 
		local m = _b[_cons]
		qui replace meanmedian`var' = _b[_cons] if study == `study'
		local s = _b[Zj]
		qui replace sd`var' = _b[Zj] if study == `study'
	
		forvalues value = 1/`max_x' {
			if `value' == `max_x' {
				local zj1 = 9999999
			}
			else {
				qui sum Zj if x == `value' & study == `study'
				local zj1 = r(mean)
			}
			if `value' == 1 {
				local zj0 = -999999
			}
			else {
				local value0 = `value'-1
				qui sum Zj if x == `value0' & study == `study'
				local zj0 = r(mean)
			}
			replace Mj_`var' = `m'+`s'*((normalden(`zj0')-normalden(`zj1'))/(normal(`zj1')-normal(`zj0'))) if study == `study' & x == `value'
		}

	}
	*Replace ncases/controls in x==1 with total cases/controls
	qui sum ncases if study==`study'
	qui replace ncases = r(sum) if study == `study'
	qui sum ncontrols if study==`study'
	qui replace ncontrols = r(sum) if study == `study'
	qui sum yearsfu if study == `study'
	qui replace yearsfu = r(sum) if study == `study'
	qui replace ntotal = ncases+ncontrols if study == `study'
	drop if x>1 & study == `study'
}

*Generate SMD
qui replace sd_pooled = sqrt(((ncases-1)*sdcases^2+(ncontrols-1)*sdcontrols^2)/(ntotal-2))

qui gen smd = (meanmediancases-meanmediancontrols)/sd_pooled
qui gen smd_se = sqrt(1/ncases+1/ncontrols+smd^2/(2*ntotal))
qui gen study_bmi_mean = (meanmediancases*ncases+meanmediancontrols*ncontrols)/ntotal
qui gen pvalue = 2*normal(-abs(smd/smd_se))

drop Nj-Mj_controls count condense unadjor-x
replace type = "MD"
gen md5 = smd*5/sd_pooled
gen md5_se = smd_se*5/sd_pooled

save "PCa MD from Cat.dta", replace
restore
drop if condense == 1
drop condense


*Extensive data cleaning to do here, study by study
*Generating study-specific IDs is probably easiest

sort endnoteid x
qui egen study = group(endnoteid outcome group) if endnote != 209014 //209014 needs to be in the categorical analysis, but is albatross otherwise
bysort study: egen count = count(study)
qui replace count = . if count == 1

*Total number of participants
qui gen ntotal = ncases+ncontrols
qui replace ntotal = ncases if ntotal ==.

*Calculate SE from 95% CI
qui egen lciadj = ends(adj95ci), punct("-") trim head
destring lciadj, replace
qui egen uciadj = ends(adj95ci), punct("-") trim tail
destring uciadj, replace

qui egen lciunadj = ends(unadj95ci), punct("-") trim head
destring lciunadj, replace
qui egen uciunadj = ends(unadj95ci), punct("-") trim tail
destring uciunadj, replace

qui replace adjse = (ln(uciadj)-ln(lciadj))/3.92
qui replace unadjse = (ln(uciunadj)-ln(lciunadj))/3.92

*Get all the x values to 1-`max' for each study
qui sum study
local max = r(max)
forvalues study = 1/`max' {
	qui sum x if study == `study'
	local dif = r(min)-1
	qui replace x = x-`dif' if study == `study'
}

*Calculate SD from 95% CI for 28720 (written in notes)
local sd1 = sqrt(6214)*(0.4)/3.92 
local sd2 = sqrt(6930)*(0.3)/3.92 
qui replace bmisd = sqrt((6214*`sd1'^2+6930*`sd2')/(6214+6930)) if endnote == 28720

*Mean years from bmi to outcome need replacing for multi-time-point studies 
qui replace meanyearsbetween = 47 if endnote == 4352 & group == "All - 20 years"
qui replace meanyearsbetween = 37 if endnote == 4352 & group == "All - 30 years"
qui replace meanyearsbetween = 27 if endnote == 4352 & group == "All - 40 years"
qui replace meanyearsbetween = 17 if endnote == 4352 & group == "All - 50 years"
qui replace meanyearsbetween = 7 if endnote == 4352 & group == "All - 60 years"
qui replace meanyearsbetween = 20 if endnote == 16570 & group == "Younger"
qui replace meanyearsbetween = 46 if endnote == 16801 & group == "BMI at age 18"
qui replace bmisd = 3.3 if endnote == 16801 & group == "BMI at age 18"
qui replace meanyearsbetween = 34 if endnote == 16801 & group == "BMI at age 30"
qui replace meanyearsbetween = 19 if endnote == 16801 & group == "BMI at age 45"
qui replace meanyearsbetween = 40 if endnote == 19122
qui replace whichvariablesare = "Age, ethnicity, family history, fat intake" if endnote == 19122
qui replace meanyearsbetween = 54 if endnote == 24392 & group == "All - age 25"
qui replace meanyearsbetween = 29 if endnote == 24392 & group == "All - age 50"
qui replace meanyearsbetween = 39 if endnote == 36541 & group == "All - BMI 30 years"
qui replace meanyearsbetween = 54 if endnote == 39525 & group == "All - 21 years"
qui replace meanyearsbetween = 33 if endnote == 42842 & group == "All - 30 years"
qui replace meanyearsbetween = 40 if endnote == 43538
qui replace meanyearsbetween = 48 if endnote == 44248 & group == "All - 25 years"
qui replace meanyearsbetween = 28 if endnote == 44248 & group == "All - 45 years"
qui replace meanyearsbetween = 8 if endnote == 44248 & group == "All - 65 years"
qui replace meanyearsbetween = 48 if endnote == 44358 & group == "All - 20 years"

order title-bmiascer notes, last
order timeofbmi-overallrisk, last
order bmisub, last
order bmimean bmimedian, a(bmisdcontrol)
order ntotal, a(ncontrol)
order study, a(endnote)

*Gen ntotal, total number of participants in a study
*Note that n needs to change to sum of ncases & ncontrols, as they may not be the same
gen ncases_study = .
gen ncontrols_study = .
qui sum study
local study_max = r(max)
forvalues i = 1/`study_max' {
	qui sum ncases if study == `i'
	qui replace ncases_study = r(sum) if study == `i'
	qui sum ncontrols if study == `i'
	qui replace ncontrols_study = r(sum) if study == `i'
}

*Some studies had f/u time and number of cases/controls estimated - mark these
qui gen fuestimated = 1 if strpos(notes, "f/u estimate")>0
qui replace fuestimated = 1 if strpos(notes, "fu estimate")>0
qui replace fuestimated = 1 if strpos(notes, "fuestimate")>0

qui gen nccestimated = 1 if strpos(notes, "ncc estimated")>0

sort study x

*Need BMI SD
*Easy - Pooled SD from SD1 SD2
qui replace bmisd = sqrt((bmisdcase^2*(allcancern-1)+bmisdcontrol^2*(nocancern-1))/(allcancern-2)) if bmisd ==.

*Hard - Pooled SD from bmirange
*Calculate case and control SD from quantile data

qui gen Xj = .
qui replace Xj = bmi_upper

qui gen Nj = ntotal
qui gen Pj = .
qui gen Zj = .
qui gen Wj = .
qui gen Mj = .
qui gen Mj_cases = .
qui gen Mj_controls = .

sort study x
qui gen lnadjor = ln(adjor)
qui gen lnunadjor = ln(unadjor)
qui gen lnadjhr = ln(adjhr)
qui gen lnunadjhr = ln(unadjhr)
qui gen lnor_a = .
qui gen lnor_a_se = .
qui gen lnor_u = .
qui gen lnor_u_se = .
qui gen lnhr_a = .
qui gen lnhr_a_se = .
qui gen lnhr_u = .
qui gen lnhr_u_se = .

gen overweight = .
gen overweight_se = .
gen obese = .
gen obese_se = .

qui gen bmimeancases = .
qui gen bmimeancontrols = .
qui gen study_bmi_mean = .
qui gen sd_pooled = .
rename bmisdcase bmisdcases
rename bmisdcontrol bmisdcontrols
qui gen ncases_total = .
qui gen ncontrols_total = .

*Need all studies to have the referent exposure level at x=1 but only AFTER the SD has been calculated
qui gen referent = .
qui replace referent = 1 if (unadjor != . | unadjhr != .) & unadjse == . 
qui replace referent = 1 if (adjor != . | adjhr != . ) & adjse == .

*Need all non-referent exposure levels to be slightly different to 1
qui replace lnadjor = 0.00001 if lnadjor == 0 & referent!=1
qui replace lnunadjor = 0.00001 if lnunadjor == 0 & referent!=1
qui replace lnadjhr = 0.00001 if lnadjhr == 0 & referent!=1
qui replace lnunadjhr = 0.00001 if lnunadjhr == 0 & referent!=1

qui sum study
local study_max = r(max)
dis "`study_max'"
forvalues study = 1/`study_max' {
	dis _newline "study = `study'"
	qui sum count if study == `study'
	if r(mean) >2 & r(mean) <. {
		qui sum x if study == `study'
		local max_x = r(max)
		
		*Total number of participants in the study
		qui sum Nj if study==`study'
		local sum = r(sum)
		
		forvalues value = 1/`max_x' {
			qui sum Nj if x <= `value' & study==`study'
			qui replace Pj = r(sum)/`sum' if x==`value' & study==`study'
		}

		qui replace Zj = invnormal(Pj)
		qui replace Wj = normalden(Zj,0,1)^2/(Pj*(1-Pj))
	
		qui sum x if study==`study'
		local quantile = r(max)-1
		qui sum Wj if study==`study'
		qui replace Wj = Wj*(`quantile'/r(sum))

		*Regress to get the mean and SD of the study, using the beta of Zj for the mean and the constant for the SD
		qui regress Xj Zj if study==`study'
		local m = _b[_cons]
		qui replace study_bmi_mean = _b[_cons] if study == `study'
		local s = _b[Zj]
		qui replace sd_pooled = _b[Zj] if study == `study'
		
		forvalues value = 1/`max_x' {
			if `value' == `max_x' {
				local zj1 = 9999999
			}
			else {
				qui sum Zj if x == `value' & study == `study'
				local zj1 = r(mean)
			}
			if `value' == 1 {
				local zj0 = -999999
			}
			else {
				local value0 = `value'-1
				qui sum Zj if x == `value0' & study == `study'
				local zj0 = r(mean)
			}
			qui replace Mj = `m'+`s'*((normalden(`zj0')-normalden(`zj1'))/(normal(`zj1')-normal(`zj0'))) if study == `study' & x == `value'
		}
		
		*Then for cases and controls (Only for BMI SD by case/control status)
		*IRRs don't have to have ncontrols, so replace with ncases as ncontrols don't exist
		foreach var in cases controls {
			if "`var'" == "controls" {
				qui sum ncontrols if study == `study'
				local controls_n_escape = r(N)
				local controls_n_escape_2 = r(min)
				if `controls_n_escape' == 0 | `controls_n_escape_2'==0 {
					local var = "cases"
					dis _newline "controls -> cases"
				}
			}
			qui replace Nj = n`var' if study == `study'
			qui replace Pj = .
			qui replace Zj = .
			qui replace Wj = .
		
			*Total number of participants in the study
			qui sum Nj if study==`study'
			local sum = r(sum)
		
			forvalues value = 1/`max_x' {
				qui sum Nj if x <= `value' & study==`study'
				qui replace Pj = r(sum)/`sum' if x==`value' & study==`study'
			}
		
			qui replace Zj = invnormal(Pj)
			qui replace Wj = normalden(Zj,0,1)^2/(Pj*(1-Pj))
	
			qui sum x if study==`study'
			local quantile = r(max)-1
			qui sum Wj if study==`study'
			qui replace Wj = Wj*(`quantile'/r(sum))
	
			*Regress to get the mean and SD of the study, using the beta of Zj for the mean and the constant for the SD
			qui regress Xj Zj if study==`study'
			local m = _b[_cons]
			qui replace bmimean`var' = _b[_cons] if study == `study'
			local s = _b[Zj]
			qui replace bmisd`var' = _b[Zj] if study == `study'
		
			forvalues value = 1/`max_x' {
				if `value' == `max_x' {
					local zj1 = 9999999
				}
				else {
					qui sum Zj if x == `value' & study == `study'
					local zj1 = r(mean)
				}
				if `value' == 1 {
					local zj0 = -999999
				}
				else {
					local value0 = `value'-1
					qui sum Zj if x == `value0' & study == `study'
					local zj0 = r(mean)
				}
				replace Mj_`var' = `m'+`s'*((normalden(`zj0')-normalden(`zj1'))/(normal(`zj1')-normal(`zj0'))) if study == `study' & x == `value'
			}
	
		}
	}
	replace Mj_cases = bmimean if endnote == 33625 & bmimean != .
	qui sum lnadjor if study == `study'
	local adj = r(N)
	qui sum yearsfu if study == `study'
	local fu = r(N)
	qui sum ncontrols if study == `study'
	local sum_con = r(N)
	qui sum ncases if study == `study'
	local sum_cases = r(N)
	qui sum lnunadjor if study == `study'
	local unadj = r(N)
	qui sum lnunadjhr if study == `study'
	local unadjhr = r(N)	
	qui sum lnadjhr if study == `study'
	local adjhr = r(N)
	
	*Switch the referent to be first
	qui sum x if referent == 1 & study == `study'
	local ref = r(sum)
	if `ref' !=1 {
		qui replace x = 0 if referent == 1 & study == `study'
		qui replace x = x+1 if x < `ref' & study == `study'
	}
	sort study x
	
	*GLST needs the reference group to be dose == 0
	foreach k in Mj Mj_cases Mj_controls {
		qui su `k' if x == 1 & study == `study'
		qui replace `k' = `k'-r(mean) if study == `study'
		
	}
	
	*Need to make a category for years f/u (ir) instead of cc
	if `adj' > 0 & `fu' == 0 {
		glst lnadjor Mj if study == `study', se(adjse) cov(ntotal ncases) cc
		qui replace lnor_a = _b[Mj] if study == `study'
		qui replace lnor_a_se = _se[Mj] if study == `study'
	}

	if `unadj' > 0 & `fu' == 0 {
		glst lnunadjor Mj if study == `study', se(unadjse) cov(ntotal ncases) cc
		qui replace lnor_u = _b[Mj] if study == `study'
		qui replace lnor_u_se = _se[Mj] if study == `study'
	}
	if `fu' > 0 & `unadjhr' > 0 {
		glst lnunadjhr Mj_cases if study == `study', se(unadjse) cov(yearsfu ncases) ir
		qui replace lnhr_u = _b[Mj] if study == `study'
		qui replace lnhr_u_se = _se[Mj] if study == `study'
	}
	if `fu' > 0 & `adjhr' > 0 {
		glst lnadjhr Mj_cases if study == `study', se(adjse) cov(yearsfu ncases) ir
		qui replace lnhr_a = _b[Mj] if study == `study'
		qui replace lnhr_a_se = _se[Mj] if study == `study'
	}	

	*Replace ncases/controls in x==1 with total cases/controls
	qui sum ncases if study==`study'
	qui replace ncases_total = r(sum) if study == `study'
	qui sum ncontrols if study==`study'
	qui replace ncontrols_total = r(sum) if study == `study'	
	qui replace ntotal = ncases_total+ncontrols_total if study == `study'
}

*Binary outcomes
destring bmi_lower, gen(bmi_l)
replace bmi_l = int(bmi_l+0.999)
gen bmi_u = int(bmi_upper+0.999)
gen bmi_cat = 0 if bmi_u == 25 & (bmi_l == . | bmi_l <=20)
replace bmi_cat = 1 if bmi_u == 30 & bmi_l == 25
replace bmi_cat = 2 if bmi_l == 30 
drop bmi_l bmi_u

*N in each BMI group
foreach var in normal overweight obese {
	gen `var'_ncases = .
	gen `var'_ncontrols = .
	gen `var'_ntotal = .
}
gen or_hr = ""

gen normal_bmi = .
gen overweight_bmi = .
gen obese_bmi = .

*Exception for 209014
qui sum study
replace study = r(max)+1 if endnote == 209014

qui sum study
local max=  r(max)
forvalues study = 1/`max' {
	foreach var of varlist adjor unadjor {
		qui sum `var' if study == `study'
		replace or_hr = "OR" if r(mean) != . & study == `study'
	}
	foreach var of varlist adjhr unadjhr {
		qui sum `var' if study == `study'
		replace or_hr = "HR" if r(mean) != . & or_hr == "" & study == `study'
	}	
	*Check if normal BMI cat has an OR of 1
	foreach var of varlist adjor adjhr unadjor unadjhr {
		local cont = 0
		qui sum `var' if study == `study' & bmi_cat == 0
		if r(mean) == 1 {
			local cont = 1
		}
		
		if `cont' == 1 {
			
			*Get Ns
			foreach n in ncases ncontrols {
				qui sum `n' if study == `study' & bmi_cat == 0
				qui replace normal_`n' = r(mean) if study == `study'
				qui sum `n' if study == `study' & bmi_cat == 1
				qui replace overweight_`n' = r(mean) if study == `study'				
				qui sum `n' if study == `study' & bmi_cat == 2
				qui replace obese_`n' = r(mean) if study == `study'
			}
			
			foreach xx in normal overweight obese {
				qui replace `xx'_ntotal = `xx'_ncases+`xx'_ncontrols
			}
			
			*Mean BMI in each group
			qui sum Mj if study == `study' & bmi_cat == 0 & or_hr == "OR"
			qui replace normal_bmi = r(mean) if study == `study' & or_hr == "OR"
			qui sum Mj if study == `study' & bmi_cat == 1 & or_hr == "OR"
			qui replace overweight_bmi = r(mean) if study == `study' & or_hr == "OR"
			qui sum Mj if study == `study' & bmi_cat == 2 & or_hr == "OR"
			qui replace obese_bmi = r(mean) if study == `study' & or_hr == "OR"
			
			qui sum Mj_cases if study == `study' & bmi_cat == 0 & or_hr == "HR"
			qui replace normal_bmi = r(mean) if study == `study' & or_hr == "HR"
			qui sum Mj_cases if study == `study' & bmi_cat == 1 & or_hr == "HR"
			qui replace overweight_bmi = r(mean) if study == `study' & or_hr == "HR"
			qui sum Mj_cases if study == `study' & bmi_cat == 2 & or_hr == "HR"
			qui replace obese_bmi = r(mean) if study == `study'	 & or_hr == "HR"		
			
			qui sum `var' if study == `study' & bmi_cat == 1
			qui replace overweight = r(mean) if study == `study' & overweight == .
			if "`var'" == "adjor" | "`var'" == "adjhr" {
				qui sum adjse if study == `study' & bmi_cat == 1
				qui replace overweight_se = r(mean) if study == `study' & overweight_se == .
			}
			else {
				qui sum unadjse if study == `study' & bmi_cat == 1
				qui replace overweight_se = r(mean) if study == `study' & overweight_se == .
			}
			
			qui sum `var' if study == `study' & bmi_cat == 2
			qui replace obese = r(mean) if study == `study' & obese == .
			
			if "`var'" == "adjor" | "`var'"== "adjhr" {
				qui sum adjse if study == `study' & bmi_cat == 2
				qui replace obese_se = r(mean) if study == `study' & obese_se == .
			}
			else {
				qui sum unadjse if study == `study' & bmi_cat == 2
				qui replace obese_se = r(mean) if study == `study' & obese_se == .
			}
		}
	}
}

keep if referent == 1
qui gen md = lnor_a
qui replace md = lnor_u if md == .
qui replace md = lnhr_a if md == .
qui replace md = lnhr_u if md == .

qui gen md_se = lnor_a_se
qui replace md_se = lnor_u_se if md_se == .
qui replace md_se = lnhr_a_se if md_se == .
qui replace md_se = lnhr_u_se if md_se == .

*SMD5
qui gen md5 = md*5
qui gen md5_se = md_se*5

*Keep only the useful variables
keep endnote albatross author year studyname studylocation overallrisk midyear ethnicity timeof meanyears meanage whichvariables* group outcome ncases* ncontrols* ///
bmimean* bmisdcases bmisdcontrol adjpvalue md md_se study_bmi ntotal sd_pooled md5 md5_se type or_hr *ncases *ncontrols *ntotal *_bmi overweight* obese*

drop ncases ncontrols

qui replace sd_pooled = bmisdcases if study_bmi<10
qui replace study_bmi = bmimeancases if study_bmi<10

drop if endnote == 24299 & outcome == "Metastatic"
drop if strpos(outcome, "Gleason")>0 | strpos(outcome, "grade")>0 | strpos(outcome,"ocali")>0 | outcome == "Non-aggressive"
replace outcome = "Advanced" if outcome == "Aggressive" | outcome == "Extraprostatic" | outcome == "Metastatic"

*Replace some studynames
qui replace studyname = "northern sweden health and disease cohort (nshdc)" if endnote == 41806
replace studyname = "multiethnic cohort" if endnote == 209014

*Get rid of duplicate studies (with the highest SEs)
qui replace studyname = lower(studyname)
sort studyname group outcome
qui egen study2 = group(studyname outcome group)
bysort study2: egen count = count(study2)
qui replace study2 = . if count == 1
qui replace count = . if count == 1 | count == 0
sort count studyname outcome group 
qui gen include = .
qui replace include = 1 if count == .
bysort study2: egen rank = rank(md5_se) if study2 !=., unique
qui replace include = 1 if count !=. & rank == 1

*Dal Maso is better than Crispo (due to adjustment)
replace include = 1 if endnote == 42842 // Dal Maso
replace include = 0 if endnote == 42687 // Crispo

*Get rid of studies with the longest f/u time (e.g. Moller)
qui replace include = 0 if endnote == 4352 & group != "All - 60 years"
qui replace include = 0 if endnote == 24299 & group != "All"
qui replace include = 0 if endnote == 24392 & group != "All - baseline"
qui replace include = 0 if endnote == 39525 & group != "All"
qui replace include = 0 if endnote == 42842 & group != "All" 

*Create categorical analysis dataset
preserve
keep if overweight !=. | obese != .
keep endnote author year outcome ntotal studyname studylocation midyear ethnicity timeof meanyears meanage whichvariables overallrisk ncases_total ncontrols_total ///
overweight* obese* normal* or_hr
order endnote studyname normal* overweight* obese*
rename ncases_total ncases
rename ncontrols_total ncontrols
*Drop duplicate studies (higher SEs in midspan & MEC)
drop if (endnote == 7831 & outcome == "All") | (endnote == 39525 & outcome == "All") 
save "PCa categorical analysis complete.dta", replace
restore

rename bmimeancases bmi_cases
rename bmimeancontrols bmi_controls
rename bmisdcases bmi_cases_sd
rename bmisdcontrols bmi_controls_sd
rename study_bmi_mean bmi_all
rename sd_pooled bmi_all_sd

drop bmimean normal_* overweight_* obese_* ncases_study ncontrols_study overweight obese study2 count rank
keep if include == 1

rename adjpvalue pvalue
replace pvalue = 2*normal(-abs(md5/md5_se)) if pvalue == .

order endnote albatross outcome md5 md5_se pvalue ncases ncontrols ntotal bmi_cases bmi_cases_sd bmi_controls bmi_controls_sd bmi_all bmi_all_sd studyname
rename ncases_total ncases
rename ncontrols_total ncontrols
save "PCa categories complete.dta", replace
}

*Continuous
{

********************************************************************************
**Continuous data**
********************************************************************************

set more off
clear
use "PCa continuous.dta", clear

gen or_hr = ""
replace or_hr = "OR" if unadjor != . | adjor != .
replace or_hr = "HR" if unadjhr != . | adjhr != .

*BMI subgroup has all the info I need
gen unit = 1 if strpos(bmisubgroup, "1")>0
replace unit = 2 if strpos(bmisubgroup, "2")>0
replace unit = 5 if strpos(bmisubgroup, "5")>0

*Get the SE
qui egen lciadj = ends(adj95ci), punct("-") trim head
destring lciadj, replace
qui egen uciadj = ends(adj95ci), punct("-") trim tail
destring uciadj, replace

qui egen lciunadj = ends(unadj95ci), punct("-") trim head
destring lciunadj, replace
qui egen uciunadj = ends(unadj95ci), punct("-") trim tail
destring uciunadj, replace

qui replace adjse = (ln(uciadj)-ln(lciadj))/3.92
qui replace unadjse = (ln(uciunadj)-ln(lciunadj))/3.92

*27155 has 99% CIs, not 95% CIs
qui replace adjse = (ln(uciadj)-ln(lciadj))/(2*2.575) if endnote == 27155
qui replace unadjse = (ln(uciunadj)-ln(lciunadj))/(2*2.575) if endnote == 27155

*Convert everything to 5 unit increases
qui gen md5 = ln(adjor)
qui replace md5 = ln(adjhr) if adjor == .

qui gen md5_se = adjse

qui replace md5 = md5*5 if unit == 1
qui replace md5_se = md5_se*5 if unit == 1
qui replace md5 = md5*2.5 if unit == 2
qui replace md5_se = md5_se*2.5 if unit == 2

qui replace ncontrols = nocancern if ncontrols == .
qui gen ntotal = ncases+ncontrols

*Choose which to include
qui gen include = 1
qui replace include = 0 if endnote == 4352 & group != "All - 60 years"
qui replace include = 0 if endnote == 7260 & group != "All"
qui replace include = 0 if endnote == 36541 & group != "All"
qui replace include = 0 if endnote == 44358 & group != "All"
qui replace include = 0 if endnote == 100004 
qui replace include = 0 if endnote == 20795 

drop if outcome == "Nonaggressive" | strpos(outcome,"Gleason")>0 | strpos(outcome,"grade")>0 | strpos(outcome,"ocali")>0
drop if outcome == "Metastatic" & endnote == 28720
replace outcome = "Advanced" if outcome == "Aggressive"
drop if include == 0

rename adjpvalue pvalue
replace pvalue = 2*normal(-abs(md5/md5_se)) if pvalue == .

keep endnote albatross author year studyname studylocation overallrisk midyear ethnicity timeof meanyears meanage whichvariables* group outcome ncases ncontrols ///
pvalue md5 md5_se ntotal include type or_hr

order endnote albatross outcome md5 md5_se pvalue ncases ncontrols ntotal studyname

save "PCa continuous complete.dta", replace

}

*MD from categorical
{
*MD from categorical results
use "PCa MD from Cat.dta", clear

*Keep only the useful variables
keep endnote albatross author year studyname studylocation overallrisk midyear ethnicity timeof meanyears meanage whichvariables* group outcome ncases ncontrols ///
meanmedian* sdcases sdcontrol pvalue smd smd_se study_bmi ntotal sd_pooled md5 md5_se type

sort endnote

rename meanmediancases bmi_cases
rename meanmediancontrols bmi_controls
rename sdcases bmi_cases_sd
rename sdcontrols bmi_controls_sd
rename study_bmi_mean bmi_all
rename sd_pooled bmi_all_sd

gen or_hr = "OR"

order endnote albatross outcome md5 md5_se pvalue ncases ncontrols ntotal bmi_cases bmi_cases_sd bmi_controls bmi_controls_sd bmi_all bmi_all_sd studyname

save "PCa MD from Cat complete.dta", replace

}

*Albatross
{
*Albatross
use "Albatross - PCa.dta",clear

gen ntotal = ncases+ncontrols, a(ncontrols)
gen e = 1 if adjhr >1 | adjor > 1
replace e = -1 if adjhr <1 | adjor < 1
rename adjpvalue pvalue

keep endnote albatross author year studyname studylocation overallrisk midyear ethnicity timeof meanyears meanage whichvariables* group outcome ncases* ncontrols* ntotal pvalue e
drop if outcome == "Localized" | outcome == "Nonagressive"
replace outcome = "Advanced" if outcome == "Aggressive"

save "Albatross - PCa complete.dta", replace
}

*Merge
{
********************************************************************************
**Merge**
********************************************************************************

set more off
clear

*First, get the continuous data up to scratch
use "PCa categories complete.dta", clear

*Append data
append using "PCa MD complete.dta"
append using "PCa continuous complete.dta"
append using "PCa MD from Cat complete.dta"
append using "IPD PCa - complete.dta"
replace author = studyname if type == "IPD"
replace timeofbmic = "Before" if author == "Krimpen"
replace timeofbmic = "Same time" if author == "PLCO" | author == "ProtecT"
replace or_hr = "OR" if type == "IPD"

*Retitle Lundqvist for each cohort
qui replace author = "Lundqvist (Younger cohort)" if author == "Lundqvist" & group == "Younger"
qui replace author = "Lundqvist (Older cohort)" if author == "Lundqvist" & group == "Older"

drop include

*Get rid of duplicate studies (with the highest SEs)
replace studyname = lower(studyname)
sort outcome studyname

*CAPS - keep the continuous estimates, but also keep the mean BMI etc.
foreach var of varlist bmi_cases-bmi_all_sd {
	qui sum `var' if endnote == 4352 & outcome == "All" & type == "MD"
	qui replace `var' = r(mean) if endnote == 4352 & outcome == "All" & type == "continuous"
	
	qui sum `var' if endnote == 4352 & outcome == "Advanced" & type == "categorical"
	qui replace `var' = r(mean) if endnote == 4352 & outcome == "Advanced" & type == "continuous"
}
drop if endnote == 4352 & type != "continuous"

*HPFS - keep the BMI closest to the outcome, and categorical over MD
drop if endnote == 209225 | endnote == 3101

*SEER - closest BMI to outcome
drop if endnote == 19122 & group != "All"

*The NCS - keep continuous, pass down BMI etc.
foreach var of varlist bmi_cases-bmi_all_sd {
	qui sum `var' if endnote == 208122 & outcome == "Advanced" 
	qui replace `var' = r(mean) if endnote == 44358 & outcome == "Advanced"
	
	qui sum `var' if endnote == 208122 & outcome == "All" 
	qui replace `var' = r(mean) if endnote == 44358 & outcome == "All" 
}
drop if endnote == 208122 & outcome == "All"

*Danish diet, cancer and health study - keep the categorical one, since the continuous uses *far* fewer people (n=688 versus 28690)
drop if (endnote == 5189 | endnote == 209072) & outcome == "All"

*Italian study - keep the categorical one
drop if endnote == 19607 & outcome == "All"

*JPHC - drop the albatross (Karahashi already has JPHC results)
drop if endnote == 11050 & outcome == "All"

*Malmo diet and cancer study - keep the categorical one
drop if (endnote == 6846 | endnote == 7584) & outcome == "All"

*NSHDC - keep the categorical one
drop if endnote == 19732 & outcome == "All"

*PHS - keep the one with the smaller SE
drop if endnote == 42921 & outcome == "All"

*Proteus - keep the categorical one
drop if endnote == 2582 & outcome == "All"

append using "Albatross - PCa complete.dta"

order e, a(md5_se)
replace e = md5 if e == .

save "PCa - pre-graph.dta", replace

********************************************************************************
}

*Graphs
use "PCa - pre-graph.dta", clear
*Giles had two results...
drop if author == "Giles" & group == "All - 21 years"
replace year = 2014 if endnote == 210062
sort outcome year
rename ncases cases
rename ncontrols controls
format cases controls %9.0f

rename year year_num
tostring year_num, gen(year)
replace year = "IPD" if type == "IPD"

*Boehm wasn't before, it was same time
replace timeof = "Same time" if author == "Boehm"

*GRAPH
metan md5 md5_se if outcome == "All" & or_hr == "HR", random second(fixed) effect("HR") eform label(namevar=author, yearvar=year) sortby(year author) xlabel(0.7, 0.8, 0.9, 1.1, 1.2, 1.3, 1.4) force textsize(110) xtitle("HR for prostate cancer per 5 kg/m{superscript:2} increase BMI", size(vsmall)) rcols(cases controls)
graph export "$cd_graphs\Final graphs\For processing\PCa All (HR).tif", replace as(tif) width(1200)
metan md5 md5_se if outcome == "All" & or_hr == "OR", random second(fixed) effect("OR") by(timeof) eform label(namevar=author, yearvar=year) sortby(year author) xlabel(0.4, 0.5, 0.66, 0.8, 1.25, 1.5, 2, 2.5) force textsize(105) xtitle("OR for prostate cancer per 5 kg/m{superscript:2} increase BMI", size(vsmall)) rcols(cases controls)
graph export "$cd_graphs\Final graphs\For processing\PCa All (OR).tif", replace as(tif) width(1200)

metafunnel md5 md5_se if outcome == "All" & or_hr == "HR" & md5_se<1, xtitle("HR for prostate cancer per 5 kg/m{superscript:2} increase in BMI") eform ytitle("SE of log(HR)")
graph export "$cd_graphs\Final graphs\For processing\PCa All - Funnel (HR).tif", replace as(tif) width(1200)
metafunnel md5 md5_se if outcome == "All" & or_hr == "OR" & md5_se<1, xtitle("OR for prostate cancer per 5 kg/m{superscript:2} increase in BMI") eform ytitle("SE of log(OR)")
graph export "$cd_graphs\Final graphs\For processing\PCa All - Funnel (OR).tif", replace as(tif) width(1200)

metan md5 md5_se if outcome == "Advanced" & or_hr == "HR", random second(fixed) effect("HR") eform label(namevar=author, yearvar=year) sortby(year author) xlabel(0.66, 0.8, 1.25, 1.5) force textsize(125) xtitle("HR for advanced prostate cancer per 5 kg/m{superscript:2} increase BMI", size(vsmall)) rcols(cases controls)
graph export "$cd_graphs\Final graphs\For processing\PCa Advanced (HR).tif", replace as(tif) width(1200)
metan md5 md5_se if outcome == "Advanced" & or_hr == "OR", random second(fixed) effect("OR") by(timeof) eform label(namevar=author, yearvar=year) sortby(year author) xlabel(0.66, 0.8, 1.25, 1.5) force textsize(125) xtitle("OR for advanced prostate cancer per 5 kg/m{superscript:2} increase BMI", size(vsmall)) rcols(cases controls)
graph export "$cd_graphs\Final graphs\For processing\PCa Advanced (OR).tif", replace as(tif) width(1200)

metafunnel md5 md5_se if outcome == "Advanced" & or_hr == "HR" & md5_se<1, xtitle("HR for advanced prostate cancer per 5 kg/m{superscript:2} increase in BMI") eform ytitle("SE of log(HR)")
graph export "$cd_graphs\Final graphs\For processing\PCa Advanced - Funnel (HR).tif", replace as(tif) width(1200)
metafunnel md5 md5_se if outcome == "Advanced" & or_hr == "OR" & md5_se<1, xtitle("OR for advanced prostate cancer per 5 kg/m{superscript:2} increase in BMI") eform ytitle("SE of log(OR)")
graph export "$cd_graphs\Final graphs\For processing\PCa Advanced - Funnel (OR).tif", replace as(tif) width(1200)

*Albatross

qui gen included = "Included - OR" if albatross != 1
qui replace included = "Included - HR" if albatross != 1 & or_hr == "HR"
qui replace outcome = "Advanced" if outcome == "Aggressive"
qui replace included = "Excluded" if albatross == 1

qui gen n = ntotal
qui replace n = cases+controls if n == .
qui gen p = pvalue
qui replace p = pvalue if p == .
qui replace p = 0.9999 if p == 1
qui replace e = 0.01 if e == 0
qui replace e = 0.1 if endnote == 11050
qui replace e = 0.1 if endnote == 8063

qui gen r = cases/controls

save "PCa complete.dta", replace

albatross n p e if outcome == "All", type(smd) contours(0.02 0.05 0.1) r(r) by(included) adjust color
graph export "$cd_graphs\Final graphs\For processing\PCa All - Albatross Adjusted.tif", replace
albatross n p e if outcome == "Advanced", type(smd) contours(0.02 0.05 0.1) r(r) by(included) adjust color
graph export "$cd_graphs\Final graphs\For processing\PCa Advanced - Albatross Adjusted.tif", replace

*Meta-regression - PCa
use "PCa complete.dta", clear

replace meanyears = 5 if author == "Krimpen"
replace meanyears = 0 if author == "ProtecT" | author == "PLCO"
replace meanage = 61.8 if author == "ProtecT"
replace meanage = 68.8 if author == "PLCO"
replace meanage = 61.1 if author == "Krimpen"
replace midyear = 2004 if author == "ProtecT"
replace midyear = 1998 if author == "PLCO"
replace midyear = 1995 if author == "Krimpen"

*Mark case-control studies
gen case_control = 0
{
qui replace case_control = 1 if endnoteid == 19122
qui replace case_control = 1 if endnoteid == 28720
qui replace case_control = 1 if endnoteid == 4352
qui replace case_control = 1 if endnoteid == 24897
qui replace case_control = 1 if endnoteid == 24897
qui replace case_control = 1 if endnoteid == 24897
qui replace case_control = 1 if endnoteid == 24897
qui replace case_control = 1 if endnoteid == 24055
qui replace case_control = 1 if endnoteid == 44677
qui replace case_control = 1 if endnoteid == 23481
qui replace case_control = 1 if endnoteid == 22642
qui replace case_control = 1 if endnoteid == 44036
qui replace case_control = 1 if endnoteid == 43538
qui replace case_control = 1 if endnoteid == 20401
qui replace case_control = 1 if endnoteid == 42842
qui replace case_control = 1 if endnoteid == 42930
qui replace case_control = 1 if endnoteid == 208335
qui replace case_control = 1 if endnoteid == 42467
qui replace case_control = 1 if endnoteid == 18681
qui replace case_control = 1 if endnoteid == 16233
qui replace case_control = 1 if endnoteid == 28720
qui replace case_control = 1 if endnoteid == 36021
qui replace case_control = 1 if endnoteid == 8029
qui replace case_control = 1 if endnoteid == 8774
qui replace case_control = 1 if endnoteid == 207948
qui replace case_control = 1 if endnoteid == 8698
qui replace case_control = 1 if endnoteid == 4634
qui replace case_control = 1 if endnoteid == 4352
qui replace case_control = 1 if endnoteid == 210051
qui replace case_control = 1 if endnoteid == 209177

*Cross-sectional
qui replace case_control = 1 if endnoteid == 18982 | author == "PLCO" | author == "ProtecT"

label define case_control 0 "Cohort/Nested case-control" 1 "Case-control/Cross-sectional"
label values case_control case_control
}

metan md5 md5_se if outcome == "All" & or_hr == "OR", random second(fixed) effect("OR") by(case_control) eform label(namevar=author, yearvar=year) sortby(year author) xlabel(0.4, 0.5, 0.66, 0.8, 1.25, 1.5, 2, 2.5) force textsize(105) xtitle("OR for prostate cancer per 5 kg/m{superscript:2} increase BMI", size(vsmall)) rcols(cases controls)
graph export "$cd_graphs\Final graphs\For processing\PCa All (OR cohort sensitivity).tif", replace as(tif) width(1200)

*Ethnicity
qui gen ethnic_bin = 1
qui replace ethnic_bin = 2 if ethnicity == "Chinese" | ethnicity == "Japan" | ethnicity == "Japanese" | ethnicity == "Japanese (100)" | ethnicity == "Korean" | ///
ethnicity == "Black" | ethnicity=="Multiethnic" | author == "Whittemore (Blacks)" | ethnicity == "Asian" | ethnicity == "Black (100%)" | ethnicity == "Chinese-American" | ///
ethnicity == "Hawaii" | ethnicity == "Iran" | ethnicity == "Japanese (100%)" | ethnicity == "Japanese-American"
capture label define ethnic_bin 1 "Caucasian" 2 "Non-Caucasian" 
label values ethnic_bin ethnic_bin
encode overallrisk, gen(rob)
order rob, a(overallrisk)
replace rob = 0 if rob == 2

*USA/Europe
gen region = 0 if strpos(studylocation,"USA") > 0 | strpos(studylocation,"Hawaii") > 0
replace region = 1 if strpos(studylocation,"Greece") > 0 | strpos(studylocation,"Czech") > 0 | strpos(studylocation,"Denmark") > 0 | strpos(studylocation,"Europe") > 0 | strpos(studylocation,"Finland") > 0 ///
 | strpos(studylocation,"Italy") > 0 | strpos(studylocation,"Netherlands") > 0 | strpos(studylocation,"Norway") > 0 | strpos(studylocation,"Scandanavia") > 0 | strpos(studylocation,"Sweden") > 0 | strpos(studylocation,"UK") > 0
label define region 0 "USA" 1 "Europe"
label values region region 
*Mean centre mid-year of recruitment, bmi_all, mean years between BMI calculation, mean age at diagnosis
preserve
foreach var of varlist midyear bmi_all meanyears meanage {
	 qui su `var'
	 replace `var' = `var' - r(mean)
}
metareg md5 ethnic_bin rob midyear bmi_all meanyears meanage region if outcome == "All" & albatross != 1 & or_hr == "HR", wsse(md5_se) eform
metareg md5 ethnic_bin rob midyear bmi_all meanyears meanage region if outcome == "All" & albatross != 1 & or_hr == "OR", wsse(md5_se) eform
*metareg md5 ethnic_bin rob midyear meanyears meanage region if outcome == "Advanced" & albatross != 1 & or_hr == "HR", wsse(md5_se) eform //not enough studies
metareg md5 rob midyear meanyears meanage region if outcome == "Advanced" & albatross != 1 & or_hr == "OR", wsse(md5_se) eform
restore

*Flip usa/europe to estimate constant for USA
preserve
replace region = region - 1
replace region = 1 if region == -1
foreach var of varlist midyear bmi_all meanyears meanage {
	 qui su `var'
	 replace `var' = `var' - r(mean)
}
metareg md5 ethnic_bin rob midyear bmi_all meanyears meanage region if outcome == "All" & albatross != 1 & or_hr == "HR", wsse(md5_se) eform
metareg md5 ethnic_bin rob midyear bmi_all meanyears meanage region if outcome == "All" & albatross != 1 & or_hr == "OR", wsse(md5_se) eform
*metareg md5 ethnic_bin rob midyear meanyears meanage region if outcome == "Advanced" & albatross != 1 & or_hr == "HR", wsse(md5_se) eform //not enough studies
metareg md5 rob midyear meanyears meanage region if outcome == "Advanced" & albatross != 1 & or_hr == "OR", wsse(md5_se) eform
restore

*Pre-PSA era 
*Note: publication year is probably the best metric to use here, as anything published 1995 and before likely had little PSA use, and anything published 2005+ had much PSA use)
gen psa_era = 0 if year_num < 2000
replace psa_era = 1 if psa_era == .
*Mean centre variables
preserve
foreach var of varlist midyear bmi_all meanyears meanage {
	 qui su `var'
	 replace `var' = `var' - r(mean)
}
metareg md5 rob psa_era bmi_all meanyears meanage if outcome == "All" & albatross != 1 & or_hr == "HR", wsse(md5_se) eform
metareg md5 rob psa_era bmi_all meanyears meanage if outcome == "All" & albatross != 1 & or_hr == "OR", wsse(md5_se) eform
*metareg md5 ethnic_bin rob psa_era meanyears meanage if outcome == "Advanced" & albatross != 1 & or_hr == "HR", wsse(md5_se) eform //not enough studies
*metareg md5 rob psa_era meanyears meanage if outcome == "Advanced" & albatross != 1 & or_hr == "OR", wsse(md5_se) eform //not enough studies
restore
*Flip psa_era to estimate constant for post-PSA era
replace psa_era = psa_era - 1
replace psa_era = 1 if psa_era == -1
preserve
foreach var of varlist midyear bmi_all meanyears meanage {
	 qui su `var'
	 replace `var' = `var' - r(mean)
}
metareg md5 rob psa_era bmi_all meanyears meanage if outcome == "All" & albatross != 1 & or_hr == "HR", wsse(md5_se) eform
metareg md5 rob psa_era bmi_all meanyears meanage if outcome == "All" & albatross != 1 & or_hr == "OR", wsse(md5_se) eform
*metareg md5 ethnic_bin rob psa_era meanyears meanage if outcome == "Advanced" & albatross != 1 & or_hr == "HR", wsse(md5_se) eform //not enough studies
*metareg md5 rob psa_era meanyears meanage if outcome == "Advanced" & albatross != 1 & or_hr == "OR", wsse(md5_se) eform //not enough studies
restore



*All PCa
metareg md5 ethnic_bin rob midyear bmi_all meanyears meanage if outcome == "All" & albatross != 1 & or_hr == "HR", wsse(md5_se) eform
metareg md5 ethnic_bin rob midyear bmi_all meanyears meanage if outcome == "All" & albatross != 1 & or_hr == "OR", wsse(md5_se) eform
*Advanced PCa
metareg md5 ethnic_bin rob midyear meanyears meanage if outcome == "Advanced" & albatross != 1 & or_hr == "HR", wsse(md5_se) eform
metareg md5 ethnic_bin rob midyear bmi_all meanyears if outcome == "Advanced" & albatross != 1 & or_hr == "HR", wsse(md5_se) eform //BMI and age can't be in model at same time
metareg md5 rob midyear meanyears meanage if outcome == "Advanced" & albatross != 1 & or_hr == "OR", wsse(md5_se) eform
metareg md5 midyear bmi_all meanyears meanage if outcome == "Advanced" & albatross != 1 & or_hr == "OR", wsse(md5_se) eform //BMI and ROB can't in the same model
*Metabias
metabias md5 md5_se if albatross != 1 & outcome == "All" & or_hr == "HR", egger
metabias md5 md5_se if albatross != 1 & outcome == "All" & or_hr == "OR", egger
metabias md5 md5_se if albatross != 1 & outcome == "Advanced" & or_hr == "HR", egger
metabias md5 md5_se if albatross != 1 & outcome == "Advanced" & or_hr == "OR", egger

*Categorical graphs
use "PCa categorical analysis complete", clear
rename obese_ncases obese_cases
rename obese_ncontrols obese_controls
rename overweight_ncases overweight_cases
rename overweight_ncontrols overweight_controls
rename normal_ncases normal_cases
rename normal_ncontrols normal_controls
replace overweight = ln(overweight)
replace obese = ln(obese)
append using "IPD PCa - complete (categorical)"
drop *logpsa*
replace type = "AD" if type == ""
replace normal_ntotal = normal_cases+normal_controls if normal_ntotal == .
replace overweight_ntotal = overweight_cases+normal_controls if overweight_ntotal == .
replace obese_ntotal = obese_cases+normal_controls if obese_ntotal == .
replace or_hr = "OR" if or_hr == ""
replace author = studyname if author == ""

*IPD isn't right for whatever reason, manual entry (it's fine for all PCa!)
replace overweight_cases = 6 if studyname == "Krimpen" & outcome == "Advanced"
replace overweight_cases = 580 if studyname == "PLCO" & outcome == "Advanced"
replace overweight_cases = 258 if studyname == "ProtecT" & outcome == "Advanced"
replace overweight_ntotal = 889 if studyname == "Krimpen" & outcome == "Advanced"
replace overweight_ntotal = 16640 if studyname == "PLCO" & outcome == "Advanced"
replace overweight_ntotal = 21103 if studyname == "ProtecT" & outcome == "Advanced"

replace obese_cases = 3 if studyname == "Krimpen" & outcome == "Advanced"
replace obese_cases = 254 if studyname == "PLCO" & outcome == "Advanced"
replace obese_cases = 108 if studyname == "ProtecT" & outcome == "Advanced"
replace obese_ntotal = 160 if studyname == "Krimpen" & outcome == "Advanced"
replace obese_ntotal = 7902 if studyname == "PLCO" & outcome == "Advanced"
replace obese_ntotal = 9424 if studyname == "ProtecT" & outcome == "Advanced"

gen n1 = normal_ntotal
gen n2 = overweight_ntotal
gen n3 = obese_ntotal

replace n1 = normal_cases+normal_controls

*Wright doesn't have numbers for N
replace n1 = . if endnote == 41001
replace n2 = . if endnote == 41001
replace n3 = . if endnote == 41001
gen n4 = n1+n2+n3
replace n4 = n1+n3 if n2 == .
gen n5 = normal_cases+overweight_cases+obese_cases if endnote != 41001
replace n5 = normal_cases + obese_cases if overweight_cases == . & endnote != 41001

format %9.0f n1-n3

replace timeof = "Before" if author == "Krimpen"
replace timeof = "Same time" if author == "ProtecT" | author == "PLCO"

gen or_hr2 = or_hr
replace or_hr2 = "OR - before" if or_hr == "OR"
replace or_hr2 = "OR - same time" if timeof == "Same time"

rename year year_num
tostring year_num, gen(year)
replace year = "IPD" if type == "IPD"

*Boehm wasn't before, it was same time
replace timeof = "Same time" if author == "Boehm"
replace or_hr2 = "OR - same time" if author == "Boehm"

metan overweight overweight_se if outcome == "All", random second(fixed) effect("HR|OR") by(or_hr2) eform label(namevar=author, yearvar=year) sortby(year author) xlabel(0.5, 0.67, 1.5, 2) force textsize(110) nooverall xtitle("HR and OR for prostate cancer, overweight versus normal weight", size(vsmall)) rcols(n1 n2)
graph export "$cd_graphs\Final graphs\For processing\PCa All - Overweight.tif", replace as(tif) width(1200)
metan obese obese_se if outcome == "All", random second(fixed) by(or_hr2) effect("HR|OR") eform label(namevar=author, yearvar=year) sortby(year author) xlabel(0.5, 0.67, 1.5, 2) force textsize(110) nooverall xtitle("HR and OR for prostate cancer, obese versus normal weight", size(vsmall)) rcols(n1 n3)
graph export "$cd_graphs\Final graphs\For processing\PCa All - Obese.tif", replace as(tif) width(1200)

metan overweight overweight_se if outcome == "Advanced", random second(fixed) effect("HR|OR") by(or_hr) eform label(namevar=author, yearvar=year) sortby(year author) xlabel(0.5, 0.67, 1.5, 2) force textsize(150) nooverall xtitle("HR and OR for advanced prostate cancer, overweight versus normal weight", size(vsmall)) rcols(n1 n2)
graph export "$cd_graphs\Final graphs\For processing\PCa Advanced - Overweight.tif", replace as(tif) width(1200)
metan obese obese_se if outcome == "Advanced", random second(fixed) effect("HR|OR") by(or_hr) eform label(namevar=author, yearvar=year) sortby(year author) xlabel(0.5, 0.67, 1.5, 2) force textsize(150) nooverall xtitle("HR and OR for advanced prostate cancer, obese versus normal weight", size(vsmall)) rcols(n1 n3)
graph export "$cd_graphs\Final graphs\For processing\PCa Advanced - Obese.tif", replace as(tif) width(1200)

*Combine OR categories
metan overweight overweight_se if outcome == "All", random by(or_hr) eform nograph
metan obese obese_se if outcome == "All", random by(or_hr) eform nograph

metan overweight overweight_se if outcome == "Advanced", random effect("HR|OR") by(or_hr) eform nograph
metan obese obese_se if outcome == "Advanced", random seffect("HR|OR") by(or_hr) eform nograph

*Extra analysis - OR for obesity (PCa all), looking at difference between IPD and non-IPD
metan obese obese_se if outcome == "All" & or_hr == "OR", random second(fixed) by(type) eform label(namevar=author, yearvar=year) sortby(year author) xlabel(0.5, 0.67, 1.5, 2) force textsize(110) nooverall xtitle("OR for prostate cancer, obese versus normal weight (IPD vs AD)", size(vsmall)) rcols(n1 n2)

*Tables
use "PCa complete.dta", clear
keep endnote author year studyname studylocation ethnicity midyear which* md5 md5_se ntotal cases p type timeof included outcome albatross e
order endnote author year studyname studylocation ethnicity midyear which* md5 ntotal p type
replace albatross = 0 if albatross == .
replace endnote = -1 if author == "ProtecT"
replace endnote = -2 if author == "PLCO"
replace endnote = -3 if author == "Krimpen"
replace endnote = 165701 if author == "Lundqvist (Older cohort)"
replace endnote = 165702 if author == "Lundqvist (Younger cohort)"
sort outcome albatross included time type
gen or = round(exp(md5),0.01)
gen l = round(exp(md5-1.96*md5_se),0.01)
gen u = round(exp(md5+1.96*md5_se),0.01)
gen x = "("
gen y = ")"
gen z = " "
gen a = " to "
egen ci = concat(or z x l a u y) if md5 != .
order ci, a(md5)
forvalues i = 1/99 {
	replace ci = subinstr(ci,"(.`i'","(0.`i'",.)
	replace ci = subinstr(ci," .`i'"," 0.`i'",.)
}

forvalues i = 1/9 {
dis "(0.`var' "
	replace ci = subinstr(ci,".`i' ",".`i'0 ",.)
	replace ci = subinstr(ci,".`i')",".`i'0)",.)
}

replace ci = subinstr(ci,"1","1.00",1) if or == 1
replace ci = subinstr(ci,"(1 ","(1.00 ",.)
replace ci = subinstr(ci," 1)"," 1.00)",.)

gen b = "0"
egen ci2 = concat(b ci) if or < 1
replace ci = ci2 if ci2 != ""

egen ntotal2 = concat(ntotal z x cases y)
order ntotal2, a(ntotal)

qui replace ci = "Positive" if e > 0 & albatross == 1
qui replace ci = "Negative" if e == -1 & albatross == 1

drop md5 md5_se l-a ci2 b or ntotal cases e

qui replace ethnicity = "" if ethnicity == studylocation | ethnicity == "Australian" | ethnicity == "Canadian" | ethnicity == "Chinese" | ethnicity == "Czech" | ethnicity == "Finnish" | ///
ethnicity == "Greece" | ethnicity == "Greek" | ethnicity == "Italian" | ethnicity == "Japanese" | ethnicity == "Korean" | ethnicity == "Norwegian" | ethnicity == "Swedish" 
replace endnote = 248971 if author == "Whittemore (Whites)"
replace endnote = 248972 if author == "Whittemore (Blacks)"
replace endnote = 248973 if author == "Whittemore (Chinese-American)"
replace endnote = 248974 if author == "Whittemore (Japanese-American)"

gsort -outcome albatross included timeof year author

save "Tables\PCa - Supplementary tables.dta", replace 

********************************************************************************
*Catgeorical table
use "PCa categorical analysis complete", clear
rename obese_ncases obese_cases
rename obese_ncontrols obese_controls
rename overweight_ncases overweight_cases
rename overweight_ncontrols overweight_controls
rename normal_ncases normal_cases
rename normal_ncontrols normal_controls
replace overweight = ln(overweight)
replace obese = ln(obese)
append using "IPD PCa - complete (categorical)"
drop *logpsa*
replace type = "AD" if type == ""
replace normal_ntotal = normal_cases+normal_controls if normal_ntotal == .
replace overweight_ntotal = overweight_cases+normal_controls if overweight_ntotal == .
replace obese_ntotal = obese_cases+normal_controls if obese_ntotal == .
replace or_hr = "OR" if or_hr == ""
replace author = studyname if author == ""

gen n1 = normal_ntotal
gen n2 = overweight_ntotal
gen n3 = obese_ntotal

replace n1 = normal_cases+normal_controls

*Wright doesn't have numbers for N
replace n1 = . if endnote == 41001
replace n2 = . if endnote == 41001
replace n3 = . if endnote == 41001

*IPD isn't right for whatever reason, manual entry (it's fine for all PCa!)
replace overweight_cases = 6 if studyname == "Krimpen" & outcome == "Advanced"
replace overweight_cases = 580 if studyname == "PLCO" & outcome == "Advanced"
replace overweight_cases = 258 if studyname == "ProtecT" & outcome == "Advanced"
replace overweight_ntotal = 889 if studyname == "Krimpen" & outcome == "Advanced"
replace overweight_ntotal = 16640 if studyname == "PLCO" & outcome == "Advanced"
replace overweight_ntotal = 21103 if studyname == "ProtecT" & outcome == "Advanced"

replace obese_cases = 3 if studyname == "Krimpen" & outcome == "Advanced"
replace obese_cases = 254 if studyname == "PLCO" & outcome == "Advanced"
replace obese_cases = 108 if studyname == "ProtecT" & outcome == "Advanced"
replace obese_ntotal = 160 if studyname == "Krimpen" & outcome == "Advanced"
replace obese_ntotal = 7902 if studyname == "PLCO" & outcome == "Advanced"
replace obese_ntotal = 9424 if studyname == "ProtecT" & outcome == "Advanced"

format %9.0f n1-n3

replace timeof = "Before" if author == "Krimpen"
replace timeof = "Same time" if author == "ProtecT" | author == "PLCO"

gen or_hr2 = or_hr
replace or_hr2 = "OR - before" if or_hr == "OR"
replace or_hr2 = "OR - same time" if timeof == "Same time"

keep endnote studyname normal* overweight* obese* author year n1-n3 or_hr2 outcome endnote
sort outcome or_hr2 year author

foreach var in overweight obese {
	gen `var'_or = exp(`var')
	gen `var'_l95 = exp(`var'-1.96*`var'_se)
	gen `var'_u95 = exp(`var'+1.96*`var'_se)
}
format %9.2f overweight_or-obese_u95

order outcome or_hr2 author normal_bmi normal_cases normal_controls overweight_bmi overweight_cases overweight_controls obese_bmi obese_cases obese_control overweight_or-obese_u95

save "Tables\PCa - Supplementary table (categorical).dta", replace 




********************************************************************************
********************************************************************************

*End of PCa analyses

********************************************************************************
********************************************************************************


}

********************************************************************************
*End of analyses
