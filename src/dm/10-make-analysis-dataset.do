capture log close
set more off
clear all

log using "~/13beaches/src/dm/10-make-analysis-dataset.log", text replace

*----------------------------------------
* 10-make-analysis-dataset.do
* ben arnold
*
* merge the appended epidemiology
* data to the appended water quality data
*
* calculate daily measures of symptom
* incidence, prevalence, and cumulative
* incidence
*
* calculate days of missed of daily
* activities due to GI illness
*
* calculate days of missed work
* due to GI illness
*
* calculate medical visits due to GI illness
*
* ...and due to all syndromes (not used in the analysis, but just for posterity)
*   Syndrome 1 (GI): stomachache, diarrhea, nausea, vomiting
*   Syndrome 2 (eyes): watery eyes, eye infection
*   Syndrome 3 (resp): sore throat, cough/runny nose, cold
*   Syndrome 4 (ears): earache
*   Syndrome 5 (uti): urinary tract infection
*   Syndrome 6 (skin): infected cuts, skin rash
*----------------------------------------

*----------------------------------------
* input files:
*	13beaches-epi.dta
*   13beaches-wq.dta
*
* output files:
*	13beaches-analysis.dta
*
*----------------------------------------


*----------------------------------------
* Load the combined epidemiology dataset
*----------------------------------------

use "~/dropbox/13beaches/data/final/13beaches-epi.dta", clear


*---------------------------------------------
* Create daily prevalence, incidence, and 
* at risk indicators for 0 to 14 days for
* each symptom
*---------------------------------------------

* calculate length of follow-up
* (clean one data entry error on a telephone date)
replace teledate = mdy(7,21,2009) if indid=="BB1001675701"

gen fudays = teledate - intdate
	label var fudays "Days between beach visit and CATI interview"

local symps "diarrhea nausea stomach vomiting rash eyeinfection earache fever urinarytractinfection cough sorethroat runnynose cold cut"

foreach symp of local symps {

	* first day of symptoms
	gen `symp'st = (`symp'stdt - intdate)
		label var `symp'st "First day of `symp'"

	* recode duration for periods < 1 day
	replace `symp'days = 1 if `symp'days == 0
	
	* last day of symptoms
	gen `symp'end = `symp'st + `symp'days - 1
		* have to comment this out b/c don't have teledate for many indivs
		* replace `symp'end = (teledate - intdate) if `symp'still==1
		replace `symp'days = (`symp'end - `symp'st + 1) if (`symp'days==. & `symp'st!=. & `symp'end!=.)
		label var `symp'end "Last day of `symp'"
	* prevalent illness
	forvalues i = 0/12 {
		gen byte `symp'p`i' = 0
		label var `symp'p`i' "Prevalent `symp', day `i'"
		replace `symp'p`i' = 1 if `i'>=`symp'st & `i'<= `symp'end
	}
	* incident illness
	forvalues i = 0/12 {
		gen byte `symp'i`i' = 0
		label var `symp'i`i' "Incident `symp', day `i'"
		replace `symp'i`i' = 1 if `i'==`symp'st
	}
	* at risk
	forvalues i = 0/12 {	
		gen byte `symp'risk`i' = 1
		label var `symp'risk`i' "At risk for `symp', day `i'"
		replace `symp'risk`i' = 0 if `i'>`symp'st
	}
	* cumulative incidence indicators
	forvalues i = 0/12 {
		gen byte `symp'ci`i' = 0
		replace `symp'ci`i' = 1 if (`i'>=`symp'st)
		label var `symp'ci`i' "Incident `symp' by day `i'"
	}
		

}


*---------------------------------------------
* Zero-out Mission Bay prevalence indicators
* and symptoms that were not measured in that
* cohort
*---------------------------------------------
local symps "diarrhea nausea stomach vomiting rash eyeinfection earache fever urinarytractinfection cough sorethroat runnynose cold cut"
foreach symp of local symps {
	* Mission Bay has no information about which days an indiv was ill,
	* so we cannot extract information about daily prevalence
	forvalues i = 0/12 {
		replace `symp'p`i' = . if beach=="Mission Bay"
	}
}

local symps "cold urinarytractinfection"
foreach symp of local symps {
	forvalues i = 0/12 {
		replace `symp'i`i' = . if beach=="Mission Bay"
		replace `symp'ci`i' = . if beach=="Mission Bay"
		replace `symp'risk`i' = . if beach=="Mission Bay"
	}
}

*---------------------------------------------
* Composite Outcome: GI Illness
* As defined in Wade 2010, Colford 2012, Arnold 2013, Yau 2014
* note that for mission bay we do not have
* information about illness duration and so 
* we cannot calculate daily prevalence.  Assume
* that GI illness started on the first date of
* any one of the symptoms
*---------------------------------------------

* Identify prevalent days
forvalues i = 0/12 {
	gen byte gip`i' = (diarrheap`i'==1) | (vomitingp`i'==1) | (nauseap`i'==1 & stomachp`i'==1) | (nauseap`i'==1 & stopdaily_gas==1) | (stomachp`i'==1 & stopdaily_gas==1)
		label var gip`i' "Prevalent GI illness, day `i'"
		replace gip`i' = . if beach=="Mission Bay"
}
* Identify the start and end date
gen gist= .
	label var gist "First day of GI illness"
forvalues i = 12(-1)0 {
	replace gist = `i' if gip`i' == 1
}
gen giend = .
	label var giend "Last day of GI illness"
forvalues i = 0/12 {
	replace giend = `i' if gip`i' == 1
}
* backfill any prevalent days
forvalues i = 0/12 {
	replace gip`i' = 1 if (`i' >= gist) & (`i' <= giend)
}
* calculate duration
gen gidays = giend - gist + 1
	label var gidays "Days of GI illness"

* Figure out the GI illness start date for
* Mission Bay -- assume it is the minimum
* start date for the GI symptoms
gen mbmin = .
forvalues i = 12(-1)0 {
	replace mbmin = `i' if (diarrheaci`i'==1 | vomitingci`i'==1 | nauseaci`i'==1 | stomachci`i'==1) & (beach=="Mission Bay")
}
replace gist = mbmin if beach=="Mission Bay"
drop mbmin

* incident illness
forvalues i = 0/12 {
	gen byte gii`i' = 0
	label var gii`i' "Incident GI illness, day `i'"
	replace gii`i' = 1 if `i'==gist
}
* days at risk
forvalues i = 0/12 {	
	gen byte girisk`i' = 1
	label var girisk`i' "At risk for GI illness, day `i'"
	replace girisk`i' = 0 if `i'> gist
}
* cumulative incidence indicators
forvalues i = 0/12 {
	gen byte gici`i' = 0
	replace gici`i' = 1 if (`i'>=gist)
	label var gici`i' "Incident GI illness by day `i'"	
}	

*---------------------------------------------
* URI (defined consistent with Wade et al. 2010)
* any 2 of sore throat, cough, runny nose, code, fever
*---------------------------------------------

* Identify prevalent days
forvalues i = 0/12 {
	egen temp = rowtotal(sorethroatp`i' coughp`i' runnynosep`i' coldp`i' feverp`i')
	gen byte urip`i' = temp > 1
		label var urip`i' "Prevalent URI, day `i'"		
	drop temp
}
* Identify the start and end date
gen urist= .
	label var urist "First day of URI"
forvalues i = 12(-1)0 {
	replace urist = `i' if urip`i' == 1
}
gen uriend = .
	label var uriend "Last day of URI"
forvalues i = 0/12 {
	replace uriend = `i' if urip`i' == 1
}
* backfill any prevalent days
forvalues i = 0/12 {
	replace urip`i' = 1 if (`i' >= urist) & (`i' <= uriend)
}
* calculate duration
gen uridays = uriend - urist + 1
	label var uridays "Days of URI"
	
* estimate the start date of URI for Mission Bay
* since we do not have daily prevalence information,
* we have to assume the earliest start date of 2+ URI symptoms
gen mbmin=.
gen counter = 0
forvalues i=0/12 {
	replace counter = counter+1 if (sorethroati`i'==1) & (beach=="Mission Bay")
	replace counter = counter+1 if (coughi`i'==1) & (beach=="Mission Bay")
	replace counter = counter+1 if (runnynosei`i'==1) & (beach=="Mission Bay")
	replace counter = counter+1 if (feveri`i'==1) & (beach=="Mission Bay") 
	replace mbmin = `i' if (counter>=2) & (mbmin==.) & (beach=="Mission Bay")
}
replace urist = mbmin if beach=="Mission Bay"
drop mbmin counter



* incident illness
forvalues i = 0/12 {
	gen byte urii`i' = 0
	label var urii`i' "Incident URI, day `i'"
	replace urii`i' = 1 if `i'==urist
}
* days at risk
forvalues i = 0/12 {	
	gen byte uririsk`i' = 1
	label var uririsk`i' "At risk for URI, day `i'"
	replace uririsk`i' = 0 if `i'> urist
}
* cumulative incidence indicators
forvalues i = 0/12 {
	gen byte urici`i' = 0
	replace urici`i' = 1 if (`i'>=urist)
	label var urici`i' "Incident URI by day `i'"
}	



*---------------------------------------------
* Count of days missed of daily activities due to
* GI illness
* Survey Question:
* "Did the [symptoms] prevent you from performing 
*  daily activities such as school, recreation, 
*  or vacation activities, or work around the home?"
*---------------------------------------------

gen byte dailygi = 0
	replace dailygi = stopdailydays_gas if (stopdailydays_gas!=.)
	* recode single value of 98 to 0
	replace dailygi = 0 if dailygi==98
	* recode values >10 to 10
	replace dailygi = 10 if dailygi>10
	label var dailygi "Days of missed activities due to GI illness"
	

*---------------------------------------------
* Count of days missed of paid work due to
* GI illness, including direct and indirect impacts
* Survey Questions:
* "When your [symptoms] began, were you working 
*  for pay either inside or outside the home? 
*  Please include jobs for which you were self employed."
* and
* "How many TOTAL days did other household members lose 
*  time from work because of [your/NAME]'s symptoms?"
*---------------------------------------------

gen byte workgi = 0
	replace workgi = stayhomedays_gas if (working==1 & stayhome_gas!=.)
	replace workgi = 0 if workgi==.
	* add secondary days missed
	replace workgi = workgi + othersmissdays_gas if (othersmissdays_gas>0 & othersmissdays_gas!=.)
	label var workgi "Days of paid work missed due to GI illness"

*---------------------------------------------
* Count of medical visits due to GI illness
* and overall
* Includes healthcare provider 
* in-person visits, phone consultations, and
* emergency room visits
* For people who reported separately for
* healthcare providers and ER visits, take the
* higher number reported
*
* For the overall medical visit count, take
* a conservative approach due to the limitations
* of the instrument, by just taking the highest
* number of visits reported for each individual
* across the different symptom categories 
* (GI, eye, respiratory, ear, uti, skin) due
* to the possibility for co-moribidity and 
* shared visits due to multiple symptoms
*---------------------------------------------
gen byte medgi = 0
	replace medgi = visitdoctimes_gas if visitdoc_gas==1
	replace medgi = visitertimes_gas if (emergencyroom_gas==1 & visitertimes_gas > medgi)
	replace medgi = 1 if inlist(medgi,0,.) & (phonedoc_gas==1)
	replace medgi = 0 if medgi ==.
	label var medgi "Num. medical visits due to GI illness, incl. phone consults and ER visits"


gen byte medvisits = 0
	replace medvisits = visitdoctimes_gas if visitdoc_gas==1
	replace medvisits = visitertimes_gas if (emergencyroom_gas==1 & visitertimes_gas > medvisits)
	
	replace medvisits = visitdoctimes_eye if (visitdoc_eye==1 & visitdoctimes_eye > medvisits)
	replace medvisits = visitertimes_eye if (emergencyroom_eye==1 & visitertimes_eye > medvisits)
	
	replace medvisits = visitdoctimes_res if (visitdoc_res==1 & visitdoctimes_res > medvisits)
	replace medvisits = visitertimes_res if (emergencyroom_res==1 & visitertimes_res > medvisits)
	
	replace medvisits = visitdoctimes_ear if (visitdoc_ear==1 & visitdoctimes_ear > medvisits)
	replace medvisits = visitertimes_ear if (emergencyroom_ear==1 & visitertimes_ear > medvisits)
	
	replace medvisits = visitdoctimes_uti if (visitdoc_uti==1 & visitdoctimes_uti > medvisits)
	replace medvisits = visitertimes_uti if (emergencyroom_uti==1 & visitertimes_uti > medvisits)
	
	replace medvisits = visitdoctimes_skn if (visitdoc_skn==1 & visitdoctimes_skn > medvisits)
	replace medvisits = visitertimes_skn if (emergencyroom_skn==1 & visitertimes_skn > medvisits)
	
	replace medvisits = 1 if inlist(medvisits,0,.) & (phonedoc_gas==1 | phonedoc_eye==1 | phonedoc_res==1 | phonedoc_ear==1 | phonedoc_uti==1 | phonedoc_skn==1)
	
	replace medvisits = 0 if medvisits==.
	label var medvisits "Num. medical visits (GI, Eye, Resp, Ear, UTI, Skin), incl. phone consults and ER"

	
*---------------------------------------------
* create final covariates
*---------------------------------------------


* create an age category variable for stratification
gen agestrat = .
replace agestrat = 1 if age<=4
replace agestrat = 2 if age>4 & age<=10
replace agestrat = 3 if age>10 & age!=.
	label define agestrat 1 "(0, 4]" 2 "(4, 10]" 3 ">10"
	label values agestrat agestrat
	label var agestrat "Age category used for stratification"
	
* create an age category variable
egen agecat = cut(age), at(0,5(10)75,200) icodes
	replace agecat=9 if age==.
	label define agecat 0 "0-4" 1 "5-14" 2 "15-24" 3 "25-34" 4 "35-44" 5 "45-54" 6 "55-64" 7 "65-74" 8 "75+" 9 "Missing"
	label values agecat agecat
	label var agecat "Age category"

* identify females
gen byte female = (sex==2)
	label var female "Female"
	
* create a GI contact at baseline analysis variable
gen byte gicontactbase = gicontact_base
	replace gicontactbase = 2 if gicontactbase==.
	label define ynm 0 "No" 1 "Yes" 2 "Missing"
	label values gicontactbase ynm
	label var gicontactbase "Contact with person with GI symptoms at enrollment"
	
	
* identify swimmers
gen swimmer = (bodycontact==1) | (headunder==1) | (swallwater==1)
	label values swimmer yesno
	label var swimmer "Individual is a swimmer"

	
*----------------------------------------
* Identify point source and non-point source
* conditions
*----------------------------------------

gen byte pointsource = 1
	replace pointsource = 0 if inlist(beach,"Doheny","Mission Bay","Malibu","Surfside")
	label values pointsource yesno
	label var pointsource "Point-source (vs. non-point source) conditions"
	order beach beachcode pointsource berm groundwater hhid indid

*----------------------------------------
* Identify marine vs. freshwater beaches
*----------------------------------------
gen byte marine = 1
	replace marine = 0 if inlist(beach,"Huntington","Silver","Washington Park","West")
	label var marine "Marine water"
	label values marine yesno

*----------------------------------------
* Merge the combined Epi and WQ datasets
*----------------------------------------
sort beach beachcode intdate
tempfile epi
save `epi'

use "~/dropbox/13beaches/data/final/13beaches-wq.dta", clear
rename coldate intdate
drop marine
sort beach beachcode intdate
tempfile wq
save `wq'

use `epi', clear
merge m:1 beach beachcode intdate using `wq'

tab beach if _merge==1
tab beach if _merge==2

* there are individuals in the Avalon, Doheny, and Mission Bay
* datasets that do not have matching wq data
* flag these individuals to exclude them
* from WQ - health analyses
gen nowq = _merge==1
	label var nowq "No water quality data available"

drop if _merge==2
drop _merge




*---------------------------------------------
* save an analysis dataset
*---------------------------------------------

compress
label data "13 Beaches NIH-R03 analysis dataset, created by 10-make-analysis.dataset.do"
saveold "~/dropbox/13beaches/data/final/13beaches-analysis.dta", replace version(12)
outsheet using "~/dropbox/13beaches/data/final/13beaches-analysis.csv", comma replace

desc, s
notes

* write a codebook to separate file
log close
log using "~/dropbox/13beaches/data/final/13beaches-analysis-codebook.txt", text replace
desc, s
notes
codebook
log close

exit




