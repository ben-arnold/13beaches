capture log close
set more off
clear all

log using "~/13beaches/src/dm/4-append-epi-data.log", text replace

*----------------------------------------
* 4-append-epi-data.do
* ben arnold
*
* append NEEAR, ADM, and Mission Bay
* epi data into a single file
*
*----------------------------------------

*----------------------------------------
* input files:
*	neear-epi.dta
*	adm-epi.dta
*	mb-epi.dta
*
* output files:
*	13beaches-epi.dta / .csv
*	13beaches-epi-varlist.dta / .csv
*----------------------------------------


*----------------------------------------
* read in NEEAR and append ADM
*----------------------------------------
use "~/dropbox/13beaches/data/final/neear-epi.dta", clear

desc, s

append using "~/dropbox/13beaches/data/final/adm-epi.dta"

desc, s

*----------------------------------------
* append Mission Bay
*----------------------------------------

append using "~/dropbox/13beaches/data/final/mb-epi.dta"

desc, s


* make beach labels more informative
replace beach = "Avalon" if beach=="AV"
replace beach = "Boqueron" if beach=="BB"
replace beach = "Doheny" if beach=="DO"
replace beach = "Edgewater" if beach=="EB"
replace beach = "Fairhope" if beach=="FB"
replace beach = "Goddard" if beach=="GB"
replace beach = "Huntington" if beach=="HB"
replace beach = "Malibu" if beach=="MA"
replace beach = "Surfside" if beach=="MB"
replace beach = "Silver" if beach=="SB"
replace beach = "West" if beach=="WB"
replace beach = "Washington Park" if beach=="WP"
replace beach = "Mission Bay" if beach=="MissionBay"

* format dates
format intdate teledate %d


* variable label that wasn't coming through for some reason
label define sanddry 1 "All dry" 2 "Mostly dry" 3 "Mostly wet" 4 "All wet" 8 "DK" 9 "MD"
label values sanddry sanddry

* add the Yes/No label to a few straggling variables that remain unlabeled
local yesnovars "eatfood cough gicontact_int fevertemptaken working othersmiss_gas othersmiss_eye othersmiss_res othersmiss_ear othersmiss_uti othersmiss_skn"
foreach var of varlist `yesnovars' {
	label values `var' yesno
}

* assume that 39 people with no GI illness at enrollment are "No"
replace gibase=0 if gibase==.

* assume that 381 people with no information about water contact had none
replace anycontact=0 if anycontact==.
* assume that 501 people with no information about body immersion contact had none
replace bodycontact=0 if bodycontact==.
* assume that 550 people with no information about head immersion contact had none
replace headunder=0 if headunder==.
* assume that 1648 people with no information about water ingestion had none
replace swallwater=0 if swallwater==.

*---------------------------------------------
* add variable notes
*---------------------------------------------

* drop any existing notes
notes drop _dta
notes drop *

* general notes
note: TS 13-Beaches Epidemiology Data, created by 4-append-epi-data.do
note: Please contact Ben Arnold (benarnold@berkeley.edu) and Tim Wade (wade.tim@epa.gov) if you have specific questions about the dataset
note: Some variables were not collected at all beaches. See 13beaches-epi-varlist.dta for a summary table of variable inclusion by beach group (NEEAR, Avalon/Doheny/Malibu, Mision Bay).
note: Doheny and Malibu have a -berm- variable that indicates when a sand berm was open to allow freshwater to flow into the ocean. See Colford et al (Water Res, 2012, 46, 2176-2186) and Arnold et al (Epidemiology, 2013, 24, 845-853) for details. 
note: Avalon has a -groundwater- variable that indicates that groundwater flow was above median. See Yau et al (Water Res, 2014, 59, 23-36) for details.


* Race variable for NEEAR
note race: "In the NEEAR beaches, hispanic individuals were all coded as non-white, which is different from Avalon/Doheney/Malibu/Mission Bay, where some people identified as white, hispanic"

* sunhr for MB
note sunhr: "sunhr and sunmin not calculable for Mission Bay -- only collected as categorical data, stored as sunhr_cat"
note sunhr_cat: "Mission Bay only collected categorical information about exposure to direct sunlight"

* stmatch
note stmatch: "Site- and Time-specific match code only created for Avalon/Doheny/Malibu beaches for merge to site- and time-specific water quality data"

* hhinc
note hhinc: "Household income categories collected at Avalon/Doheny/Malibu"
note hhinc_mb: "Household income categories collected at Mission Bay, different from Av/Do/Ma beaches"


*----------------------------------------
* save combined epidemiology dataset
*----------------------------------------
order beach beachcode hhid indid intdate teledate
compress
label data "13 beaches epi data, created by 4-append-epi-data.do"
saveold "~/dropbox/13beaches/data/final/13beaches-epi.dta", replace version(12)
outsheet using "~/dropbox/13beaches/data/final/13beaches-epi.csv", comma replace

* write a codebook to separate file
log close
log using "~/dropbox/13beaches/data/final/13beaches-epi-codebook.txt", text replace
desc, s
notes
codebook
log close


*---------------------------------------------
* variable reconciliation
*---------------------------------------------

log using "~/13beaches/src/dm/4-append-epi-data.log", text append


use "~/dropbox/13beaches/data/temp/neear-epi-vars.dta", clear
sort name
tempfile neear
save `neear'

use "~/dropbox/13beaches/data/temp/adm-epi-vars.dta", clear
sort name
merge 1:1 name using `neear'

gen byte inNEEAR = inlist(_merge,2,3)
	label var inNEEAR "In NEEAR data"
gen byte inADM = inlist(_merge,1,3)
	label var inADM "In ADM data"
	drop _merge

save `neear', replace

use "~/dropbox/13beaches/data/temp/mb-epi-vars.dta", clear
sort name
merge 1:1 name using `neear'
gen byte inMB = inlist(_merge,1,3)
	label var inMB "In Mission Bay data"
	drop _merge

* set unmatched values to 0
replace inNEEAR = 0 if inNEEAR==.
replace inADM = 0 if inADM==.


* tabulate discordance
tab inNEEAR inADM
tab inNEEAR inMB
tab inADM inMB

* list discordance
list name if inNEEAR & !inADM
list name if inADM & !inNEEAR

list name if (inNEEAR & inADM) & (!inMB)
list name if inMB & !(inNEEAR | inADM)


label data "13 beaches variable list, created by 4-append-epi-data.do"
saveold "~/dropbox/13beaches/data/final/13beaches-epi-varlist.dta", replace version(12)
outsheet using "~/dropbox/13beaches/data/final/13beaches-epi-varlist.csv", comma replace




log close
exit




