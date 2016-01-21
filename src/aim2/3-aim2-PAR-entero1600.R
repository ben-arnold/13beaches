



# --------------------------------------
# 3-aim2-PAR-entero1600.R
# ben arnold (benarnold@berkeley.edu)
#
# description:
# estimate the Diarrhea risk attributable
# to swimming in water that exceeds
# 35 CFU/100ml based on the Entero EPA 1600
# indicator.  Estimate the effect for 
# the overall population, E(Y|R,W) - E(Y|R=0,W),
# which is the population attributable risk (PAR).
# also divide the PAR by the baseline incidence to
# estimate the population attributable fraction (PAF)
#
# Analyses are conducted for EPA 1600
# Among Swimmers with Body Immmersion
#
# The exposure categories are <=35 and >35 CFU/100ml
#
# For diarrhea episodes and days of regular activities
# missed due to GI illness, we have stratified the 
# estimation by Age.  Days missed from paid work
# and medical visits associated with GI illness
# are too rare to stratify the data by age, unfortunately.
#
# --------------------------------------

# --------------------------------------
# input files:
#	13beaches-analysis.csv
# 
# output files:
#	aim2-PARentero1600-diar.RData
#	aim2-PARentero1600-gi.RData
#	aim2-PARentero1600-daily.RData
#	aim2-PARentero1600-workgi.RData
#	aim2-PARentero1600-medgi.RData
# --------------------------------------


# --------------------------------------
# preamble
# --------------------------------------

rm(list=ls())
library(sandwich)
library(lmtest)

# source the base functions
source("~/13beaches/src/aim1/0-aim1-base-functions.R")
source("~/13beaches/src/aim2/0-aim2-base-functions.R")

# set the number of iterations for the bootstraps
Nboot <- 1000

# --------------------------------------
# load the and pre-preprocess the 
# analysis dataset
# (refer to the base functions script
# for details on the pre-processing)
# --------------------------------------

ad <- preprocess.13beaches("~/dropbox/13beaches/data/final/13beaches-analysis.csv")

# drop individuals with no water quality information
table(ad$nowq)
ad <- subset(ad,nowq==0)
dim(ad)

# subset to non-missing exposure categories
table(ad$bodycontact)
ad <- subset(ad,ad$bodycontact=="Yes")
	dim(ad)
ad <- subset(ad,is.na(ad$entero35)==FALSE)
	dim(ad)



# --------------------------------------
# Summarize 
# the number of swimmers and number
# of cumulative incident episodes by
# Entero 1600 category
# --------------------------------------
calcNs <- function(Y,A) {
  		n <- tapply(Y,A,sum)
  		N <- tapply(Y,A,function(x) length(x))
  		cbind(n,N)
}


# --------------------------------------
# Diarrhea, Episodes
# --------------------------------------
# cumulative incident cases and individuals at risk
# by exposure category
N.all       <- calcNs(ad$diarrheaci10,ad$entero35)
N.age0to4   <- calcNs(ad$diarrheaci10[ad$agestrat=="(0, 4]"],ad$entero35[ad$agestrat=="(0, 4]"])
N.age5to10  <- calcNs(ad$diarrheaci10[ad$agestrat=="(4, 10]"],ad$entero35[ad$agestrat=="(4, 10]"])
N.age11plus <- calcNs(ad$diarrheaci10[ad$agestrat==">10"],ad$entero35[ad$agestrat==">10"])
N.ps        <- calcNs(ad$diarrheaci10[ad$pointsource=="Yes"],ad$entero35[ad$pointsource=="Yes"])
N.nps       <- calcNs(ad$diarrheaci10[ad$pointsource=="No"], ad$entero35[ad$pointsource=="No"])

# Estimate and bootstrap PAR & PAF
# see the aim2 base functions for details of the bootAR function
set.seed(346435)
fmla <- formula(diarrheaci10~entero35+marine+pointsource+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach)

# all ages and stratified by age category
AR35cfu.diar <- bootAR(fn="AR35cfu",fmla=fmla,data=ad,ID=ad$hhid,strata=ad$beach,iter=Nboot)

	# remove age category from the adjustment covariates
	fmla <- formula(diarrheaci10~entero35+marine+pointsource+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach)


AR35cfu.diar.0to4 <- bootAR(fn="AR35cfu",fmla=fmla,data=subset(ad,agestrat=="(0, 4]"),ID=ad$hhid[ad$agestrat=="(0, 4]"],strata=ad$beach[ad$agestrat=="(0, 4]"],iter=Nboot)

AR35cfu.diar.5to10 <- bootAR(fn="AR35cfu",fmla=fmla,data=subset(ad,agestrat=="(4, 10]"),ID=ad$hhid[ad$agestrat=="(4, 10]"],strata=ad$beach[ad$agestrat=="(4, 10]"],iter=Nboot)

AR35cfu.diar.11plus <- bootAR(fn="AR35cfu",fmla=fmla,data=subset(ad,agestrat==">10"),ID=ad$hhid[ad$agestrat==">10"],strata=ad$beach[ad$agestrat==">10"],iter=Nboot)

	# add age category to the adjustment covariates, remove pointsource and marine
	fmla <- formula(diarrheaci10~entero35+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach)

AR35cfu.diar.ps <- bootAR(fn="AR35cfu",fmla=fmla,data=subset(ad,pointsource=="Yes"),ID=ad$hhid[ad$pointsource=="Yes"],strata=ad$beach[ad$pointsource=="Yes"],iter=Nboot)

AR35cfu.diar.nps <- bootAR(fn="AR35cfu",fmla=fmla,data=subset(ad,pointsource=="No"),ID=ad$hhid[ad$pointsource=="No"],strata=ad$beach[ad$pointsource=="No"],iter=Nboot)

# save results down
save(N.all,N.age0to4,N.age5to10,N.age11plus,N.ps,N.nps,AR35cfu.diar,AR35cfu.diar.0to4,AR35cfu.diar.5to10,AR35cfu.diar.11plus,AR35cfu.diar.ps,AR35cfu.diar.nps,file="~/dropbox/13beaches/aim2-results/rawoutput/aim2-PARentero1600-diar.RData")


# --------------------------------------
# GI Illness, Episodes
# --------------------------------------
# cumulative incident cases and individuals at risk
# by exposure category
N.all       <- calcNs(ad$gici10,ad$entero35)
N.age0to4   <- calcNs(ad$gici10[ad$agestrat=="(0, 4]"],ad$entero35[ad$agestrat=="(0, 4]"])
N.age5to10  <- calcNs(ad$gici10[ad$agestrat=="(4, 10]"],ad$entero35[ad$agestrat=="(4, 10]"])
N.age11plus <- calcNs(ad$gici10[ad$agestrat==">10"],ad$entero35[ad$agestrat==">10"])
N.ps        <- calcNs(ad$gici10[ad$pointsource=="Yes"],ad$entero35[ad$pointsource=="Yes"])
N.nps       <- calcNs(ad$gici10[ad$pointsource=="No"], ad$entero35[ad$pointsource=="No"])

# Estimate and bootstrap PAR & PAF
# see the aim2 base functions for details of the bootAR function
set.seed(346435)
fmla <- formula(gici10~entero35+marine+pointsource+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach)

# all ages and stratified by age category
AR35cfu.gi <- bootAR(fn="AR35cfu",fmla=fmla,data=ad,ID=ad$hhid,strata=ad$beach,iter=Nboot)

	# remove age category from the adjustment covariates
	fmla <- formula(gici10~entero35+marine+pointsource+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach)


AR35cfu.gi.0to4 <- bootAR(fn="AR35cfu",fmla=fmla,data=subset(ad,agestrat=="(0, 4]"),ID=ad$hhid[ad$agestrat=="(0, 4]"],strata=ad$beach[ad$agestrat=="(0, 4]"],iter=Nboot)

AR35cfu.gi.5to10 <- bootAR(fn="AR35cfu",fmla=fmla,data=subset(ad,agestrat=="(4, 10]"),ID=ad$hhid[ad$agestrat=="(4, 10]"],strata=ad$beach[ad$agestrat=="(4, 10]"],iter=Nboot)

AR35cfu.gi.11plus <- bootAR(fn="AR35cfu",fmla=fmla,data=subset(ad,agestrat==">10"),ID=ad$hhid[ad$agestrat==">10"],strata=ad$beach[ad$agestrat==">10"],iter=Nboot)

	# add age category to the adjustment covariates, remove pointsource and marine
	fmla <- formula(gici10~entero35+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach)

AR35cfu.gi.ps <- bootAR(fn="AR35cfu",fmla=fmla,data=subset(ad,pointsource=="Yes"),ID=ad$hhid[ad$pointsource=="Yes"],strata=ad$beach[ad$pointsource=="Yes"],iter=Nboot)

AR35cfu.gi.nps <- bootAR(fn="AR35cfu",fmla=fmla,data=subset(ad,pointsource=="No"),ID=ad$hhid[ad$pointsource=="No"],strata=ad$beach[ad$pointsource=="No"],iter=Nboot)

# save results down
save(N.all,N.age0to4,N.age5to10,N.age11plus,N.ps,N.nps,AR35cfu.gi,AR35cfu.gi.0to4,AR35cfu.gi.5to10,AR35cfu.gi.11plus,AR35cfu.gi.ps,AR35cfu.gi.nps,file="~/dropbox/13beaches/aim2-results/rawoutput/aim2-PARentero1600-gi.RData")


# --------------------------------------
# Days missed of daily activities
# --------------------------------------
# cumulative incident cases and individuals at risk
# by exposure category
N.all       <- calcNs(ad$dailygi,ad$entero35)
N.age0to4   <- calcNs(ad$dailygi[ad$agestrat=="(0, 4]"],ad$entero35[ad$agestrat=="(0, 4]"])
N.age5to10  <- calcNs(ad$dailygi[ad$agestrat=="(4, 10]"],ad$entero35[ad$agestrat=="(4, 10]"])
N.age11plus <- calcNs(ad$dailygi[ad$agestrat==">10"],ad$entero35[ad$agestrat==">10"])
N.ps        <- calcNs(ad$dailygi[ad$pointsource=="Yes"],ad$entero35[ad$pointsource=="Yes"])
N.nps       <- calcNs(ad$dailygi[ad$pointsource=="No"], ad$entero35[ad$pointsource=="No"])


# Estimate and bootstrap PAR & PAF
# see the aim2 base functions for details of the bootAR function
set.seed(76857)
fmla <- formula(dailygi~entero35+marine+pointsource+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach)

# all ages and stratified by age category
AR35cfu.dailygi <- bootAR(fn="AR35cfu",fmla=fmla,data=ad,ID=ad$hhid,strata=ad$beach,iter=Nboot)

	# remove age category from the adjustment covariates
	fmla <- formula(dailygi~entero35+marine+pointsource+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach)

AR35cfu.dailygi.0to4 <- bootAR(fn="AR35cfu",fmla=fmla,data=subset(ad,agestrat=="(0, 4]"),ID=ad$hhid[ad$agestrat=="(0, 4]"],strata=ad$beach[ad$agestrat=="(0, 4]"],iter=Nboot)

AR35cfu.dailygi.5to10 <- bootAR(fn="AR35cfu",fmla=fmla,data=subset(ad,agestrat=="(4, 10]"),ID=ad$hhid[ad$agestrat=="(4, 10]"],strata=ad$beach[ad$agestrat=="(4, 10]"],iter=Nboot)

AR35cfu.dailygi.11plus <- bootAR(fn="AR35cfu",fmla=fmla,data=subset(ad,agestrat==">10"),ID=ad$hhid[ad$agestrat==">10"],strata=ad$beach[ad$agestrat==">10"],iter=Nboot)

	# add age category to the adjustment covariates, remove pointsource and marine
	fmla <- formula(dailygi~entero35+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach)

AR35cfu.dailygi.ps <- bootAR(fn="AR35cfu",fmla=fmla,data=subset(ad,pointsource=="Yes"),ID=ad$hhid[ad$pointsource=="Yes"],strata=ad$beach[ad$pointsource=="Yes"],iter=Nboot)

AR35cfu.dailygi.nps <- bootAR(fn="AR35cfu",fmla=fmla,data=subset(ad,pointsource=="No"),ID=ad$hhid[ad$pointsource=="No"],strata=ad$beach[ad$pointsource=="No"],iter=Nboot)

# save results down
save(N.all,N.age0to4,N.age5to10,N.age11plus,N.ps,N.nps,AR35cfu.dailygi,AR35cfu.dailygi.0to4,AR35cfu.dailygi.5to10,AR35cfu.dailygi.11plus,AR35cfu.dailygi.ps,AR35cfu.dailygi.nps,file="~/dropbox/13beaches/aim2-results/rawoutput/aim2-PARentero1600-dailygi.RData")


# --------------------------------------
# Days missed of paid work
# --------------------------------------
# cumulative incident cases and individuals at risk
# by exposure category
N.all       <- calcNs(ad$workgi,ad$entero35)
N.age0to4   <- calcNs(ad$workgi[ad$agestrat=="(0, 4]"],ad$entero35[ad$agestrat=="(0, 4]"])
N.age5to10  <- calcNs(ad$workgi[ad$agestrat=="(4, 10]"],ad$entero35[ad$agestrat=="(4, 10]"])
N.age11plus <- calcNs(ad$workgi[ad$agestrat==">10"],ad$entero35[ad$agestrat==">10"])
N.ps        <- calcNs(ad$workgi[ad$pointsource=="Yes"],ad$entero35[ad$pointsource=="Yes"])
N.nps       <- calcNs(ad$workgi[ad$pointsource=="No"], ad$entero35[ad$pointsource=="No"])


# Estimate and bootstrap PAR & PAF
# see the aim2 base functions for details of the bootAR function

# the outcome is too rare to stratify by age
set.seed(2342)
fmla <- formula(workgi~entero35+marine+pointsource+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach)
AR35cfu.workgi <- bootAR(fn="AR35cfu",fmla=fmla,data=ad,ID=ad$hhid,strata=ad$beach,iter=Nboot)

# save results down
save(N.all,N.age0to4,N.age5to10,N.age11plus,N.ps,N.nps,AR35cfu.workgi,file="~/dropbox/13beaches/aim2-results/rawoutput/aim2-PARentero1600-workgi.RData")

# --------------------------------------
# Medical visits
# --------------------------------------
# cumulative incident cases and individuals at risk
# by exposure category
N.all       <- calcNs(ad$medgi,ad$entero35)
N.age0to4   <- calcNs(ad$medgi[ad$agestrat=="(0, 4]"],ad$entero35[ad$agestrat=="(0, 4]"])
N.age5to10  <- calcNs(ad$medgi[ad$agestrat=="(4, 10]"],ad$entero35[ad$agestrat=="(4, 10]"])
N.age11plus <- calcNs(ad$medgi[ad$agestrat==">10"],ad$entero35[ad$agestrat==">10"])
N.ps        <- calcNs(ad$medgi[ad$pointsource=="Yes"],ad$entero35[ad$pointsource=="Yes"])
N.nps       <- calcNs(ad$medgi[ad$pointsource=="No"], ad$entero35[ad$pointsource=="No"])


# Estimate and bootstrap PAR & PAF
# see the aim2 base functions for details of the bootAR function

# the outcome is too rare to stratify by age
set.seed(746356)
fmla <- formula(medgi~entero35+marine+pointsource+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach)
AR35cfu.medgi <- bootAR(fn="AR35cfu",fmla=fmla,data=ad,ID=ad$hhid,strata=ad$beach,iter=Nboot)

# save results down
save(N.all,N.age0to4,N.age5to10,N.age11plus,N.ps,N.nps,AR35cfu.medgi,file="~/dropbox/13beaches/aim2-results/rawoutput/aim2-PARentero1600-medgi.RData")




