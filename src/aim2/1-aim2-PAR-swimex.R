



# --------------------------------------
# 1-aim2-PAR-swimex.R
# ben arnold (benarnold@berkeley.edu)
#
# description:
# estimate the Diarrhea risk attributable
# to swimming in water using swimmers and
# non-swimmers.  Estimate the effect for 
# the overall population, E(Y|R,W) - E(Y|R=0,W),
# which is the population attributable risk (PAR).
# also divide the PAR by the baseline incidence to
# estimate the population attributable fraction (PAF)
#
#
# The exposure categories are Body Immersion
# vs. no water contact.
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
#	aim2-PARswimex-diar.RData
#	aim2-PARswimex-gi.RData
#	aim2-PARswimex-daily.RData
#	aim2-PARswimex-workgi.RData
#	aim2-PARswimex-medgi.RData
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

# set the number of iterations for the bootsraps
Nboot <- 1000

# --------------------------------------
# load the and pre-preprocess the 
# analysis dataset
# (refer to the base functions script
# for details on the pre-processing)
# --------------------------------------

ad <- preprocess.13beaches("~/dropbox/13beaches/data/final/13beaches-analysis.csv")

# --------------------------------------
# Calculate the actual Ns for each cell
# and store them for potential use in 
# plotting and tables
# --------------------------------------
# Total population, Non-swimmers and body immersion swimmers
# cumulative incident cases and individuals at risk
N.swimex <- function(Y,age) {
	# Y   : binary outcome variable to summarize
	# age : age category (factor)
	calcNs <- function(x) {
  		n <- sum(x)
  		N <- length(x)
  		cbind(n,N)
	}
	N.all       <- calcNs(Y)
	N.age0to4   <- calcNs(Y[age=="(0, 4]"])
	N.age5to10  <- calcNs(Y[age=="(4, 10]"])
	N.age11plus <- calcNs(Y[age==">10"])
	Ns <- rbind(N.all,N.age0to4,N.age5to10,N.age11plus)
	rownames(Ns) <- c("All Ages","Age 0 to 4","Age 5 to 10","Age >10")
	return(Ns)
}


# --------------------------------------
# Diarrhea Episodes
# --------------------------------------
# Total population, Non-swimmers and body immersion swimmers
# cumulative incident cases and individuals at risk
N.total  <- N.swimex(ad$diarrheaci10,ad$agestrat)
N.noswim <- N.swimex(ad$diarrheaci10[ad$anycontact=="No"],ad$agestrat[ad$anycontact=="No"])
N.body   <- N.swimex(ad$diarrheaci10[ad$bodycontact=="Yes"],ad$agestrat[ad$bodycontact=="Yes"])

# Estimate and bootstrap PAR & PAF
# see the aim2 base functions for details of the bootAR function
set.seed(74835)
fmla <- formula(diarrheaci10~anycontact+bodycontact+marine+pointsource+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach)

# all ages and stratified by age category
PARswimex.diar <- bootAR(fn="ARswimex",fmla=fmla,data=ad,ID=ad$hhid,strata=ad$beach,iter=Nboot)

	# remove age category from the adjustment covariates
	fmla <- formula(diarrheaci10~anycontact+bodycontact+marine+pointsource+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach)

PARswimex.diar.0to4 <- bootAR(fn="ARswimex",fmla=fmla,data=subset(ad,agestrat=="(0, 4]"),ID=ad$hhid[ad$agestrat=="(0, 4]"],strata=ad$beach[ad$agestrat=="(0, 4]"],iter=Nboot)

PARswimex.diar.5to10 <- bootAR(fn="ARswimex",fmla=fmla,data=subset(ad,agestrat=="(4, 10]"),ID=ad$hhid[ad$agestrat=="(4, 10]"],strata=ad$beach[ad$agestrat=="(4, 10]"],iter=Nboot)

PARswimex.diar.11plus <- bootAR(fn="ARswimex",fmla=fmla,data=subset(ad,agestrat==">10"),ID=ad$hhid[ad$agestrat==">10"],strata=ad$beach[ad$agestrat==">10"],iter=Nboot)

# save results down
save(N.total,N.noswim,N.body,PARswimex.diar,PARswimex.diar.0to4,PARswimex.diar.5to10,PARswimex.diar.11plus,file="~/dropbox/13beaches/aim2-results/rawoutput/aim2-PARswimex-diar.RData")

# --------------------------------------
# GI Illness Episodes
# --------------------------------------
# Total population, Non-swimmers and body immersion swimmers
# cumulative incident cases and individuals at risk
N.total  <- N.swimex(ad$gici10,ad$agestrat)
N.noswim <- N.swimex(ad$gici10[ad$anycontact=="No"],ad$agestrat[ad$anycontact=="No"])
N.body   <- N.swimex(ad$gici10[ad$bodycontact=="Yes"],ad$agestrat[ad$bodycontact=="Yes"])

# Estimate and bootstrap PAR & PAF
# see the aim2 base functions for details of the bootAR function
set.seed(74835)
fmla <- formula(gici10~anycontact+bodycontact+marine+pointsource+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach)

# all ages and stratified by age category
PARswimex.gi <- bootAR(fn="ARswimex",fmla=fmla,data=ad,ID=ad$hhid,strata=ad$beach,iter=Nboot)

	# remove age category from the adjustment covariates
	fmla <- formula(gici10~anycontact+bodycontact+marine+pointsource+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach)

PARswimex.gi.0to4 <- bootAR(fn="ARswimex",fmla=fmla,data=subset(ad,agestrat=="(0, 4]"),ID=ad$hhid[ad$agestrat=="(0, 4]"],strata=ad$beach[ad$agestrat=="(0, 4]"],iter=Nboot)

PARswimex.gi.5to10 <- bootAR(fn="ARswimex",fmla=fmla,data=subset(ad,agestrat=="(4, 10]"),ID=ad$hhid[ad$agestrat=="(4, 10]"],strata=ad$beach[ad$agestrat=="(4, 10]"],iter=Nboot)

PARswimex.gi.11plus <- bootAR(fn="ARswimex",fmla=fmla,data=subset(ad,agestrat==">10"),ID=ad$hhid[ad$agestrat==">10"],strata=ad$beach[ad$agestrat==">10"],iter=Nboot)

# save results down
save(N.total,N.noswim,N.body,PARswimex.gi,PARswimex.gi.0to4,PARswimex.gi.5to10,PARswimex.gi.11plus,file="~/dropbox/13beaches/aim2-results/rawoutput/aim2-PARswimex-gi.RData")


# --------------------------------------
# Days missed of daily activities
# Associated with GI Illness
# --------------------------------------
# Total population, Non-swimmers and body immersion swimmers
# cumulative incident cases and individuals at risk
N.total  <- N.swimex(ad$dailygi,ad$agestrat)
N.noswim <- N.swimex(ad$dailygi[ad$anycontact=="No"],ad$agestrat[ad$anycontact=="No"])
N.body   <- N.swimex(ad$dailygi[ad$bodycontact=="Yes"],ad$agestrat[ad$bodycontact=="Yes"])

# Estimate and bootstrap PAR & PAF
# see the aim2 base functions for details of the bootAR function
set.seed(35435)
fmla <- formula(dailygi~anycontact+bodycontact+marine+pointsource+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach)

# all ages and stratified by age category
PARswimex.dailygi <- bootAR(fn="ARswimex",fmla=fmla,data=ad,ID=ad$hhid,strata=ad$beach,iter=Nboot)

	# remove age category from the adjustment covariates
	fmla <- formula(dailygi~anycontact+bodycontact+marine+pointsource+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach)

PARswimex.dailygi.0to4 <- bootAR(fn="ARswimex",fmla=fmla,data=subset(ad,agestrat=="(0, 4]"),ID=ad$hhid[ad$agestrat=="(0, 4]"],strata=ad$beach[ad$agestrat=="(0, 4]"],iter=Nboot)

PARswimex.dailygi.5to10 <- bootAR(fn="ARswimex",fmla=fmla,data=subset(ad,agestrat=="(4, 10]"),ID=ad$hhid[ad$agestrat=="(4, 10]"],strata=ad$beach[ad$agestrat=="(4, 10]"],iter=Nboot)

PARswimex.dailygi.11plus <- bootAR(fn="ARswimex",fmla=fmla,data=subset(ad,agestrat==">10"),ID=ad$hhid[ad$agestrat==">10"],strata=ad$beach[ad$agestrat==">10"],iter=Nboot)


# save results down
save(N.total,N.noswim,N.body,PARswimex.dailygi,PARswimex.dailygi.0to4,PARswimex.dailygi.5to10,PARswimex.dailygi.11plus,file="~/dropbox/13beaches/aim2-results/rawoutput/aim2-PARswimex-dailygi.RData")

# --------------------------------------
# Days missed of paid work
# Associated with GI Illness
# --------------------------------------
# Total population, Non-swimmers and body immersion swimmers
# cumulative incident cases and individuals at risk
N.total  <- N.swimex(ad$workgi,ad$agestrat)
N.noswim <- N.swimex(ad$workgi[ad$anycontact=="No"],ad$agestrat[ad$anycontact=="No"])
N.body   <- N.swimex(ad$workgi[ad$bodycontact=="Yes"],ad$agestrat[ad$bodycontact=="Yes"])

# Estimate and bootstrap PAR & PAF
# see the aim2 base functions for details of the bootAR function
set.seed(5234)
fmla <- formula(workgi~anycontact+bodycontact+marine+pointsource+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach)

# the outcome is too rare to stratify by age
PARswimex.workgi <- bootAR(fn="ARswimex",fmla=fmla,data=ad,ID=ad$hhid,strata=ad$beach,iter=Nboot)

# save results down
save(N.total,N.noswim,N.body,PARswimex.workgi,file="~/dropbox/13beaches/aim2-results/rawoutput/aim2-PARswimex-workgi.RData")

# --------------------------------------
# Medical visits
# Associated with GI Illness
# --------------------------------------
# Total population, Non-swimmers and body immersion swimmers
# cumulative incident cases and individuals at risk
N.total  <- N.swimex(ad$medgi,ad$agestrat)
N.noswim <- N.swimex(ad$medgi[ad$anycontact=="No"],ad$agestrat[ad$anycontact=="No"])
N.body   <- N.swimex(ad$medgi[ad$bodycontact=="Yes"],ad$agestrat[ad$bodycontact=="Yes"])

# Estimate and bootstrap PAR & PAF
# see the aim2 base functions for details of the bootAR function
set.seed(234234)
fmla <- formula(medgi~anycontact+bodycontact+marine+pointsource+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach)

# the outcome is too rare to stratify by age
PARswimex.medgi <- bootAR(fn="ARswimex",fmla=fmla,data=ad,ID=ad$hhid,strata=ad$beach,iter=Nboot)

# save results down
save(N.total,N.noswim,N.body,PARswimex.medgi,file="~/dropbox/13beaches/aim2-results/rawoutput/aim2-PARswimex-medgi.RData")










