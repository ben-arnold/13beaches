capture log close
set more off
clear all

log using "~/13beaches/src/dm/2-format-adm-epi.log", text replace

*----------------------------------------
* 2-format-adm-epi.do
* ben arnold
*
* append and format the Avalon, Doheny, 
* & Malibu datasets so that they can be 
* appended to the other beaches datasets 
* (NEEAR, Mission Bay). 
*
* some variables need to be renamed to
* be consistent with the NEEAR data.
* use the NEEAR data as the name template
*
* version 3 (23 apr 2015)
* updated doheny beach codes to separate out site E
*
* version 2 (21 feb 2015)
* retained beachcode to merge to water quality data
* 
* version 1 (20 Feb 2014)
*
*----------------------------------------

*----------------------------------------
* input files:
*	avalon2.dta
*   doheny2.dta
*	malibu2.dta
*
* output files:
*	adm-epi.dta
*
*----------------------------------------


*----------------------------------------
* read in the dataset for each beach
* restrict to completed interviews
* label the beach
*----------------------------------------


use "~/dropbox/beaches/avalon/data/final/avalon2.dta", clear
keep if pout==1
capture drop beach
gen beach = "AV"
	label var beach "Study beach"
tempfile avalon
save `avalon'


use "~/dropbox/beaches/doheny/data/final/doheny2.dta", clear
keep if pout==1
capture drop beach
gen beach = "DO"
	label var beach "Study beach"
drop psid
rename personid psid
tempfile doheny
save `doheny'

use "~/dropbox/beaches/malibu/data/final/malibu2.dta", clear
keep if pout==1
capture drop beach
gen beach = "MA"
	label var beach "Study beach"
	
	drop recrdt

append using `avalon'
append using `doheny'


*----------------------------------------
* Create a standard beach code for Avalon
* to differentiate sites A/B/C from site D
* In Doheny and Malibu, differentiate site C
* (in the lagoon) from the other sites
* And in Doheny differentiate site E from
* the others because it was nearly 1 mile
* south on the beach
*----------------------------------------
gen str beachcode = ""
	replace beachcode = "Avalon-ABC" if beach=="AV" & siteid!=4
	replace beachcode = "Avalon-D" if beach=="AV" & siteid==4
	replace beachcode = "Doheny-ABD" if beach=="DO" & (siteid!=3 & siteid!=5)
	replace beachcode = "Doheny-C" if beach=="DO" & siteid==3
	replace beachcode = "Doheny-E" if beach=="DO" & siteid==5
	replace beachcode= "Malibu-ABDE" if beach=="MA" & siteid!=3
	replace beachcode= "Malibu-C" if beach=="MA" & siteid==3
	label var beachcode "Water quality sampling location beach code"
order beach beachcode


*----------------------------------------
* rename variables and subset to be
* consistent with the NEEAR data
*----------------------------------------

* make unique household and individual IDs
replace hhid = beach + hhid
	label var hhid "Household ID"

gen indid = beach + string(psid,"%12.0g")
	label var indid "Individual ID"

order beach hhid indid

rename coldate intdate
	label var intdate "Beach interview date"
rename catidt teledate
	label var teledate "Phone interview date"

drop gast
label var gichron "Chronic GI problems/Crohn's/IBS"

rename ageyrx age
label var age "age in years"

gen watertime = (total*60 + water)
	label var watertime "Time in water (mins)"
	replace watertime = 0 if anycontact==0
	
label var asth "chronic resp prob"

rename windy wsurf
rename cut cutbase
rename eggs eggs_base
rename fish fish_base
rename shellfishlist fish_int
rename raw rawmeat_base
rename rawmeatlist rawmeat_int
rename eggslist eggs_int

rename sick gicontact_base
rename sun sunbase
rename wearnoseplugs noseplugs
rename weareyegoggles mask
rename wearearplugs earplugs
drop pool
rename pooly pool
rename publicpooly publicpool
rename differentbeachy differentbeach
rename prob skinchron

gen swam = aftdays >=1
	label var swam "went bathing or swimming anywhere since the beach interview?"


rename racecat1 race

label var samebeach "Swim/wade: same beach"

replace sex = 2 if sex==7

rename animy anim_base
rename animalcontactlisty anim_int
gen anim_any = (anim_base==1)|(anim_int==1)
	label var anim_any "Any animal contact"

gen eggs_any = (eggs_base==1)|(eggs_int==1)
	label var eggs_any "Any undercooked eggs"
	
gen fish_any = (fish_base==1)|(fish_int==1)
	label var fish_any "Any raw fish"
	
gen rawmeat_any = (rawmeat_base==1)|(rawmeat_int==1)
	label var rawmeat_any "Any raw meat"


drop nausea cough diarrhea earache
rename anysunburnlist sunburnlist

rename siteid swimloc

* minor data cleaning
replace vomitingstill = 0 if inlist(vomitingstill,6,9)

* standardize direct sunlight in minutes
rename dsun sunmin
replace sunmin = sunhr / 60


* remove "x" recode suffixes on variables
local recodes "algaewash algaemouth wetsuit sandwash sunburndays sanddry sandmouth sunburnstill"
foreach var of local recodes {
	rename `var'x `var'
}


* recode variables from 2=0
local vlist "block bsand fevertemptaken food pool preg *still *list *allergy *_int stayhome_* stopdaily_* othermiss_* phonedoc_* visitdoc_* emergencyroom_* hospitalized_* anyprescrdrug_* anyotcmeds_* earplugs noseplugs mask wetsuit dig sandwash sandmouth algae algaemouth algaewash cutbase prothat shade repel eatfood drinksy anim_base gicontact_base fish_base rawmeat_base eggs_base working asth skinchron swam samebeach publicpool differentbeach"
label define yesno 0 "No" 1 "Yes"
foreach var of varlist `vlist' {
	recode `var' 2=0
	label values `var' yesno
}

* remove "list" suffixes from symptoms
local symps "cold cough cut diarrhea earache eyeinfection fever headache nausea rash runnynose sorethroat stomach sunburn urinarytractinfection vomiting wateryeyes"
foreach var of local symps {
	capture drop `var'
	rename `var'list `var'
}

*----------------------------------------
* process the medical information into
* the format used in NEEAR
*----------------------------------------
* create variable shells for each symptom
drop *_gas
local stubs "gas eye res ear uti skn"
local slabs "GI_illness EYE_infec Respiratory Ear_infec UTI Skin_rash"
local vlist "stayhome stayhomedays stopdaily stopdailydays othermiss othermissdays phonedoc visitdoc visitdoctimes emergencyroom visitertimes hospitalized hospitalizeddays anyprescrdrug ownmoneyfordrugs anyotcmeds ownmoneyforotc"
local i = 1
foreach stub of local stubs {
	local slab = word("`slabs'",`i')
	foreach var of local vlist {
		qui gen byte `var'_`stub' = .
			local vlab : var label `var'_1
			local vlab2 = substr("`vlab'",strpos("`vlab'",":")+1,.)
			local vlab3 = "`vlab2'" + " (`slab')"
			di as res "`var'_`stub': `vlab3'"
			label var `var'_`stub' "`vlab3'"
			local vallab : value label `var'_1
			label values `var'_`stub' `vallab'
	}
	local i = `i'+1
}

* fill in the variables using data in the 6 symptom bins
* loop over the 6 symptom bins
forvalues i = 1/6 {
	di as res "----------------"_n "syndrome bin `i'" _n "----------------" _n
	* loop over the different stubs
	forvalues j = 1/6 {
		local stub = word("`stubs'",`j')
		di as res "----------------" _n "syndrome type: `stub'" _n "----------------" _n
		* replace across all relevant variables
		foreach var of local vlist {
			replace `var'_`stub' = `var'_`i' if (synd_no_`i' == `j')
		}
	}
}

* rename hospitalizeddays hospitaldays for comparability to NEEAR
* rename othermiss and othermissdays for comparability to NEEAR
foreach stub of local stubs {
	rename hospitalizeddays_`stub' hospitaldays_`stub'
	rename othermiss_`stub' othersmiss_`stub'
	rename othermissdays_`stub' othersmissdays_`stub'
}


* minor data cleaning
replace stopdaily_eye = 0 if stopdaily_eye==9
replace stopdaily_res = 0 if stopdaily_res==9
replace stopdaily_ear = 0 if stopdaily_ear==9
replace stopdaily_uti = 0 if stopdaily_uti==9
replace stopdaily_skn = 0 if stopdaily_skn==9

* recode missing values for aggregate symptoms to 0 for filter No questions
rename anycuts anycut 
rename anystomachache anystomach
local symps "diarrhea nausea vomiting fever headache sorethroat cough cold runnynose earache wateryeyes eyeinfection cut rash sunburn stomach"
foreach symp of local symps {
	replace `symp'=0 if any`symp'==2
}
* no filter for urinary tract infection, so assume missing values are no
replace urinarytractinfection=0 if urinarytractinfection==.



*----------------------------------------
* subset the dataset to relevant variables
*----------------------------------------
# delimit ;
keep beach beachcode indid hhid intdate teledate hhmem age sex 
watertime 
earplugs noseplugs mask wetsuit
sunbase sunhr sunmin dig bsand sanddry sandmouth sandwash algae algaemouth algaewash
block prothat shade repel  
eatfood food drinksy come pool publicpool
anim_base anim_int anim_any 
gicontact_base 
fish_base fish_int fish_any
rawmeat_base rawmeat_int rawmeat_any
eggs_base eggs_int eggs_any
rawfood
asth gichron skinchron allergy preg 
anycontact bodycontact headunder mouthwater swallwater wave
diarrhea diarrheanumber nausea vomiting vomitingnumber urinarytractinfection fever feverstill feverdays fevertemp* headache sorethroat cough cold runnynose earache wateryeyes eyeinfection cut rash sunburn stomach
*still *days *allergy
working 
hisp race
hhinc
*_gas *_eye *_res *_ear *_uti *_skn
*stdt
stopdaily_gas
gibase cutbase utibase vomitbase sorebase earbase eyebase rashbase
watertime
siteid stmatch
swam samebeach differentbeach
groundwater berm
;
# delimit cr

* drop variables not contained in NEEAR
drop hcgi3* uridays befdays aftdays fudays siteid2



*---------------------------------------------
* output variable information for harmonization
*---------------------------------------------
preserve
desc, replace clear
label data "variables for adm-epi.dta"
save "~/dropbox/13beaches/data/temp/adm-epi-vars.dta", replace
restore



*----------------------------------------
* save the data for combining it with
* other beaches
*----------------------------------------
label data "Avalon, Doheny, Malibu epi data, created by 2-format-adm-epi.do"
save "~/dropbox/13beaches/data/final/adm-epi.dta", replace

codebook, c

log close


* write a codebook
log using "~/dropbox/13beaches/data/final/adm-epi-codebook.txt", text replace
desc, s
*aorder
codebook
log close




exit





