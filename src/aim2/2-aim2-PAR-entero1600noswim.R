



# --------------------------------------
# 2-aim2-PAR-entero1600noswim.R
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
#	aim2-PARentero1600noswim-diar.RData
#	aim2-PARentero1600noswim-gi.RData
#	aim2-PARentero1600noswim-daily.RData
#	aim2-PARentero1600noswim-workgi.RData
#	aim2-PARentero1600noswim-medgi.RData
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

# subset to non-missing exposure categories
table(ad$swim35)
table(is.na(ad$swim35))
ad <- subset(ad,is.na(ad$swim35)==FALSE)
	dim(ad)



# --------------------------------------
# Summarize 
# the number of swimmers and number
# of cumulative incident episodes by
# exposure category
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
N.all       <- calcNs(ad$diarrheaci10,ad$swim35)
N.age0to4   <- calcNs(ad$diarrheaci10[ad$agestrat=="(0, 4]"],ad$swim35[ad$agestrat=="(0, 4]"])
N.age5to10  <- calcNs(ad$diarrheaci10[ad$agestrat=="(4, 10]"],ad$swim35[ad$agestrat=="(4, 10]"])
N.age11plus <- calcNs(ad$diarrheaci10[ad$agestrat==">10"],ad$swim35[ad$agestrat==">10"])
N.ps        <- calcNs(ad$diarrheaci10[ad$pointsource=="Yes"],ad$swim35[ad$pointsource=="Yes"])
N.nps       <- calcNs(ad$diarrheaci10[ad$pointsource=="No"], ad$swim35[ad$pointsource=="No"])

# Estimate and bootstrap PAR & PAF
# see the aim2 base functions for details of the bootAR function
set.seed(346435)
fmla <- formula(diarrheaci10~swim35+marine+pointsource+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach)

# all ages and stratified by age category
AR35noswim.diar <- bootAR(fn="AR35noswim",fmla=fmla,data=ad,ID=ad$hhid,strata=ad$beach,iter=Nboot)

	# remove age category from the adjustment covariates
	fmla <- formula(diarrheaci10~swim35+marine+pointsource+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach)


AR35noswim.diar.0to4 <- bootAR(fn="AR35noswim",fmla=fmla,data=subset(ad,agestrat=="(0, 4]"),ID=ad$hhid[ad$agestrat=="(0, 4]"],strata=ad$beach[ad$agestrat=="(0, 4]"],iter=Nboot)

AR35noswim.diar.5to10 <- bootAR(fn="AR35noswim",fmla=fmla,data=subset(ad,agestrat=="(4, 10]"),ID=ad$hhid[ad$agestrat=="(4, 10]"],strata=ad$beach[ad$agestrat=="(4, 10]"],iter=Nboot)

AR35noswim.diar.11plus <- bootAR(fn="AR35noswim",fmla=fmla,data=subset(ad,agestrat==">10"),ID=ad$hhid[ad$agestrat==">10"],strata=ad$beach[ad$agestrat==">10"],iter=Nboot)

	# add age category to the adjustment covariates, remove pointsource and marine
	fmla <- formula(diarrheaci10~swim35+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach)

AR35noswim.diar.ps <- bootAR(fn="AR35noswim",fmla=fmla,data=subset(ad,pointsource=="Yes"),ID=ad$hhid[ad$pointsource=="Yes"],strata=ad$beach[ad$pointsource=="Yes"],iter=Nboot)

AR35noswim.diar.nps <- bootAR(fn="AR35noswim",fmla=fmla,data=subset(ad,pointsource=="No"),ID=ad$hhid[ad$pointsource=="No"],strata=ad$beach[ad$pointsource=="No"],iter=Nboot)

# save results down
save(N.all,N.age0to4,N.age5to10,N.age11plus,N.ps,N.nps,AR35noswim.diar,AR35noswim.diar.0to4,AR35noswim.diar.5to10,AR35noswim.diar.11plus,AR35noswim.diar.ps,AR35noswim.diar.nps,file="~/dropbox/13beaches/aim2-results/rawoutput/aim2-PARentero1600noswim-diar.RData")


# --------------------------------------
# GI Illness, Episodes
# --------------------------------------
# cumulative incident cases and individuals at risk
# by exposure category
N.all       <- calcNs(ad$gici10,ad$swim35)
N.age0to4   <- calcNs(ad$gici10[ad$agestrat=="(0, 4]"],ad$swim35[ad$agestrat=="(0, 4]"])
N.age5to10  <- calcNs(ad$gici10[ad$agestrat=="(4, 10]"],ad$swim35[ad$agestrat=="(4, 10]"])
N.age11plus <- calcNs(ad$gici10[ad$agestrat==">10"],ad$swim35[ad$agestrat==">10"])
N.ps        <- calcNs(ad$gici10[ad$pointsource=="Yes"],ad$swim35[ad$pointsource=="Yes"])
N.nps       <- calcNs(ad$gici10[ad$pointsource=="No"], ad$swim35[ad$pointsource=="No"])

# Estimate and bootstrap PAR & PAF
# see the aim2 base functions for details of the bootAR function
set.seed(346435)
fmla <- formula(gici10~swim35+marine+pointsource+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach)

# all ages and stratified by age category
AR35noswim.gi <- bootAR(fn="AR35noswim",fmla=fmla,data=ad,ID=ad$hhid,strata=ad$beach,iter=Nboot)

	# remove age category from the adjustment covariates
	fmla <- formula(gici10~swim35+marine+pointsource+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach)


AR35noswim.gi.0to4 <- bootAR(fn="AR35noswim",fmla=fmla,data=subset(ad,agestrat=="(0, 4]"),ID=ad$hhid[ad$agestrat=="(0, 4]"],strata=ad$beach[ad$agestrat=="(0, 4]"],iter=Nboot)

AR35noswim.gi.5to10 <- bootAR(fn="AR35noswim",fmla=fmla,data=subset(ad,agestrat=="(4, 10]"),ID=ad$hhid[ad$agestrat=="(4, 10]"],strata=ad$beach[ad$agestrat=="(4, 10]"],iter=Nboot)

AR35noswim.gi.11plus <- bootAR(fn="AR35noswim",fmla=fmla,data=subset(ad,agestrat==">10"),ID=ad$hhid[ad$agestrat==">10"],strata=ad$beach[ad$agestrat==">10"],iter=Nboot)

	# add age category to the adjustment covariates, remove pointsource and marine
	fmla <- formula(gici10~swim35+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach)

AR35noswim.gi.ps <- bootAR(fn="AR35noswim",fmla=fmla,data=subset(ad,pointsource=="Yes"),ID=ad$hhid[ad$pointsource=="Yes"],strata=ad$beach[ad$pointsource=="Yes"],iter=Nboot)

AR35noswim.gi.nps <- bootAR(fn="AR35noswim",fmla=fmla,data=subset(ad,pointsource=="No"),ID=ad$hhid[ad$pointsource=="No"],strata=ad$beach[ad$pointsource=="No"],iter=Nboot)

# save results down
save(N.all,N.age0to4,N.age5to10,N.age11plus,N.ps,N.nps,AR35noswim.gi,AR35noswim.gi.0to4,AR35noswim.gi.5to10,AR35noswim.gi.11plus,AR35noswim.gi.ps,AR35noswim.gi.nps,file="~/dropbox/13beaches/aim2-results/rawoutput/aim2-PARentero1600noswim-gi.RData")


# --------------------------------------
# Days missed of daily activities
# --------------------------------------
# cumulative incident cases and individuals at risk
# by exposure category
N.all       <- calcNs(ad$dailygi,ad$swim35)
N.age0to4   <- calcNs(ad$dailygi[ad$agestrat=="(0, 4]"],ad$swim35[ad$agestrat=="(0, 4]"])
N.age5to10  <- calcNs(ad$dailygi[ad$agestrat=="(4, 10]"],ad$swim35[ad$agestrat=="(4, 10]"])
N.age11plus <- calcNs(ad$dailygi[ad$agestrat==">10"],ad$swim35[ad$agestrat==">10"])
N.ps        <- calcNs(ad$dailygi[ad$pointsource=="Yes"],ad$swim35[ad$pointsource=="Yes"])
N.nps       <- calcNs(ad$dailygi[ad$pointsource=="No"], ad$swim35[ad$pointsource=="No"])


# Estimate and bootstrap PAR & PAF
# see the aim2 base functions for details of the bootAR function
set.seed(76857)
fmla <- formula(dailygi~swim35+marine+pointsource+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach)

# all ages and stratified by age category
AR35noswim.dailygi <- bootAR(fn="AR35noswim",fmla=fmla,data=ad,ID=ad$hhid,strata=ad$beach,iter=Nboot)

	# remove age category from the adjustment covariates
	fmla <- formula(dailygi~swim35+marine+pointsource+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach)

AR35noswim.dailygi.0to4 <- bootAR(fn="AR35noswim",fmla=fmla,data=subset(ad,agestrat=="(0, 4]"),ID=ad$hhid[ad$agestrat=="(0, 4]"],strata=ad$beach[ad$agestrat=="(0, 4]"],iter=Nboot)

AR35noswim.dailygi.5to10 <- bootAR(fn="AR35noswim",fmla=fmla,data=subset(ad,agestrat=="(4, 10]"),ID=ad$hhid[ad$agestrat=="(4, 10]"],strata=ad$beach[ad$agestrat=="(4, 10]"],iter=Nboot)

AR35noswim.dailygi.11plus <- bootAR(fn="AR35noswim",fmla=fmla,data=subset(ad,agestrat==">10"),ID=ad$hhid[ad$agestrat==">10"],strata=ad$beach[ad$agestrat==">10"],iter=Nboot)

	# add age category to the adjustment covariates, remove pointsource and marine
	fmla <- formula(dailygi~swim35+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach)

AR35noswim.dailygi.ps <- bootAR(fn="AR35noswim",fmla=fmla,data=subset(ad,pointsource=="Yes"),ID=ad$hhid[ad$pointsource=="Yes"],strata=ad$beach[ad$pointsource=="Yes"],iter=Nboot)

AR35noswim.dailygi.nps <- bootAR(fn="AR35noswim",fmla=fmla,data=subset(ad,pointsource=="No"),ID=ad$hhid[ad$pointsource=="No"],strata=ad$beach[ad$pointsource=="No"],iter=Nboot)

# save results down
save(N.all,N.age0to4,N.age5to10,N.age11plus,N.ps,N.nps,AR35noswim.dailygi,AR35noswim.dailygi.0to4,AR35noswim.dailygi.5to10,AR35noswim.dailygi.11plus,AR35noswim.dailygi.ps,AR35noswim.dailygi.nps,file="~/dropbox/13beaches/aim2-results/rawoutput/aim2-PARentero1600noswim-dailygi.RData")


# --------------------------------------
# Days missed of paid work
# --------------------------------------
# cumulative incident cases and individuals at risk
# by exposure category
N.all       <- calcNs(ad$workgi,ad$swim35)
N.age0to4   <- calcNs(ad$workgi[ad$agestrat=="(0, 4]"],ad$swim35[ad$agestrat=="(0, 4]"])
N.age5to10  <- calcNs(ad$workgi[ad$agestrat=="(4, 10]"],ad$swim35[ad$agestrat=="(4, 10]"])
N.age11plus <- calcNs(ad$workgi[ad$agestrat==">10"],ad$swim35[ad$agestrat==">10"])
N.ps        <- calcNs(ad$workgi[ad$pointsource=="Yes"],ad$swim35[ad$pointsource=="Yes"])
N.nps       <- calcNs(ad$workgi[ad$pointsource=="No"], ad$swim35[ad$pointsource=="No"])


# Estimate and bootstrap PAR & PAF
# see the aim2 base functions for details of the bootAR function

# the outcome is too rare to stratify by age
set.seed(2342)
fmla <- formula(workgi~swim35+marine+pointsource+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach)
AR35noswim.workgi <- bootAR(fn="AR35noswim",fmla=fmla,data=ad,ID=ad$hhid,strata=ad$beach,iter=Nboot)

# save results down
save(N.all,N.age0to4,N.age5to10,N.age11plus,N.ps,N.nps,AR35noswim.workgi,file="~/dropbox/13beaches/aim2-results/rawoutput/aim2-PARentero1600noswim-workgi.RData")

# --------------------------------------
# Medical visits
# --------------------------------------
# cumulative incident cases and individuals at risk
# by exposure category
N.all       <- calcNs(ad$medgi,ad$swim35)
N.age0to4   <- calcNs(ad$medgi[ad$agestrat=="(0, 4]"],ad$swim35[ad$agestrat=="(0, 4]"])
N.age5to10  <- calcNs(ad$medgi[ad$agestrat=="(4, 10]"],ad$swim35[ad$agestrat=="(4, 10]"])
N.age11plus <- calcNs(ad$medgi[ad$agestrat==">10"],ad$swim35[ad$agestrat==">10"])
N.ps        <- calcNs(ad$medgi[ad$pointsource=="Yes"],ad$swim35[ad$pointsource=="Yes"])
N.nps       <- calcNs(ad$medgi[ad$pointsource=="No"], ad$swim35[ad$pointsource=="No"])


# Estimate and bootstrap PAR & PAF
# see the aim2 base functions for details of the bootAR function

# the outcome is too rare to stratify by age
set.seed(746356)
fmla <- formula(medgi~swim35+marine+pointsource+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach)
AR35noswim.medgi <- bootAR(fn="AR35noswim",fmla=fmla,data=ad,ID=ad$hhid,strata=ad$beach,iter=Nboot)

# save results down
save(N.all,N.age0to4,N.age5to10,N.age11plus,N.ps,N.nps,AR35noswim.medgi,file="~/dropbox/13beaches/aim2-results/rawoutput/aim2-PARentero1600noswim-medgi.RData")




