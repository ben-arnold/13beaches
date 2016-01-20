



# --------------------------------------
# 11-aim1-sens-enteroQPCR-470cce-regs-noswim-negcontrol.R
# ben arnold (benarnold@berkeley.edu)
#
# description:
# estimate the association between water
# quality indicator concentrations and
# the risk of Diarrhea among swimmers
# for the 13 beaches study
#
# Analyses are conducted for EPA 1600
# Among Non-Swimmers (a negative control analysis)
#
# The exposure categories are <=35 and >35 CFU/100ml
#
#
# NOTE: this script is a direct excerpt from 
# 5-aim1-enteroQPCR-470cce-regs-body.R
# for the overall pooled analyses.
#
# --------------------------------------

# --------------------------------------
# preamble
# --------------------------------------

rm(list=ls())
library(sandwich)
library(lmtest)

# source the base functions
source("~/13beaches/src/aim1/0-aim1-base-functions.R")

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
# to make the robust CI calcs work
table(ad$anycontact)
ad <- subset(ad,ad$anycontact=="No")
	dim(ad)
ad <- subset(ad,is.na(ad$entero470)==FALSE)
	dim(ad)


# --------------------------------------
# Calculate the actual Ns for each cell
# and store them for plotting and tables
# --------------------------------------

# --------------------------------------
# Summarize 
# the number of swimmers and number
# of cumulative incident episodes by
# Entero 1600 category
# --------------------------------------
calcNs <- function(outcome,exposurecat) {
  n <- tapply(outcome,exposurecat,sum)
  N <- tapply(outcome,exposurecat,function(x) length(x))
  cbind(n,N)
}

# overall
N.all <- calcNs(ad$diarrheaci10,ad$entero470)
# age stratified
N.age0to4   <- calcNs(ad$diarrheaci10[ad$agestrat=="(0, 4]"],ad$entero470[ad$agestrat=="(0, 4]"] )
N.age5to10  <- calcNs(ad$diarrheaci10[ad$agestrat=="(4, 10]"],ad$entero470[ad$agestrat=="(4, 10]"] )
N.age11plus <- calcNs(ad$diarrheaci10[ad$agestrat==">10"],ad$entero470[ad$agestrat==">10"] )



# --------------------------------------
# Calculate unadjusted 
# Cumulative Incidence Rates
# and robust 95% CIs using a saturated, 
# unadjusted model
# (just using the models to get the robust
# SEs that account for HH clustering)
# --------------------------------------


#### Overall Estimates (pooled over beach type)

# All Ages
all.cifit <- glm(diarrheaci10~entero470,family=poisson(link="log"),data=ad)
	allci.VC <- cl(ad,fm=all.cifit,cluster=ad$hhid)
	allci <- coeftest(all.cifit, allci.VC)
	allci

# Ages 0 - 4
age0to4.cifit <- glm(diarrheaci10~entero470,family=poisson(link="log"),data=ad[ad$agestrat=="(0, 4]",])
	age0to4ci.VC <- cl(ad[ad$agestrat=="(0, 4]",],fm=age0to4.cifit,cluster=ad$hhid[ad$agestrat=="(0, 4]"])
	age0to4ci <- coeftest(age0to4.cifit, age0to4ci.VC)
	age0to4ci

# Ages 5 - 10
age5to10.cifit <- glm(diarrheaci10~entero470,family=poisson(link="log"),data=ad[ad$agestrat=="(4, 10]",])
	age5to10ci.VC <- cl(ad[ad$agestrat=="(4, 10]",],fm=age5to10.cifit,cluster=ad$hhid[ad$agestrat=="(4, 10]"])
	age5to10ci <- coeftest(age5to10.cifit, age5to10ci.VC)
	age5to10ci

# Ages >10
age11plus.cifit <- glm(diarrheaci10~entero470,family=poisson(link="log"),data=ad[ad$agestrat==">10",])
	age11plusci.VC <- cl(ad[ad$agestrat==">10",],fm=age11plus.cifit,cluster=ad$hhid[ad$agestrat==">10"])
	age11plusci <- coeftest(age11plus.cifit, age11plusci.VC)
	age11plusci



# function to get Estimates and CIs from a linear combination of regression coefficients
lccalc <- function(lc,x,vcv) {
	# lc : linear combination of coefficients
	# x : log-linear model object returned from coeftest (class=coeftest)
	# vcv : variance-covariance matrix of coefficients for robust SEs
	est <- exp(t(lc)%*%x[,1])
	se  <- sqrt( t(lc)%*%vcv%*%lc )
	lb <- exp(log(est)-1.96*se)
	ub <- exp(log(est)+1.96*se)
	res <- c(est,lb,ub)
	return(res)
}
getCI <- function(fit,vcv) {
	# fit : an object returned from coeftest w/ 2 parameters corresponding to <=35/>35 CFU cumulative incidence
	# vcv : variance-covariance matrix for coefficients in the fit object
	lc1 <- lc2 <- rep(0,nrow(fit))
	lc1[c(1)  ] <- 1
	lc2[c(1,2)] <- 1
	lcs <- list(lc1,lc2)
	res <- t(sapply(lcs,lccalc,x=fit,vcv=vcv))
	colnames(res) <- c("CI","CIlb","CIub")
	rownames(res) <- c("Entero<=35cfu","Entero>35cfu")
	return(res)
}


ci.all    <- getCI(allci,allci.VC)
ci.0to4   <- getCI(age0to4ci,age0to4ci.VC)
ci.5to10  <- getCI(age5to10ci,age5to10ci.VC)
ci.11plus <- getCI(age11plusci,age11plusci.VC)


# --------------------------------------
# esimate risk of Diarrhea associated
# with exposure to EPA qPCR above/below 470 CCE/100ml
# All ages
# --------------------------------------


# --------------------------------------
# Pooled estimates


# all beaches
all.fit <- glm(diarrheaci10~entero470 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,family=poisson(link="log"),data=ad)
	all.VC <- cl(ad,fm=all.fit,cluster=ad$hhid)
	overall.fit <- coeftest(all.fit, all.VC)
	summary(all.fit)
	overall.fit

# --------------------------------------
# Age-stratified estimates and LR tests of
# interaction

# reduced models for LR tests of indicator x age interactions
agestrat.ref <- glm(diarrheaci10~entero470 +agestrat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,family=poisson(link="log"),data=ad)

# Pooled estimate (Age 0-4, 5-10, >10), Body Immersion
agestrat.fit <- glm(diarrheaci10~ entero470*agestrat +female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,family=poisson(link="log"),data=ad)
	# commented out robust SE calcs b/c not used
	# agestrat.VC <- cl(ad,fm=agestrat.fit,cluster=ad$hhid)
	# agestrat.body <- coeftest(agestrat.fit, agestrat.VC) 
	lrtest(agestrat.ref,agestrat.fit)



# --------------------------------------
# Stratified Models 
# based on tests of interaction (above)
# stratify the results by non-point and
# point source conditions
# --------------------------------------

# Age-stratified results
age0to4   <- mpreg(diarrheaci10~entero470 +female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ad[ad$agestrat=="(0, 4]",],vcv=T)
age5to10  <- mpreg(diarrheaci10~entero470 +female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ad[ad$agestrat=="(4, 10]",],vcv=T)
age11plus <- mpreg(diarrheaci10~entero470 +female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ad[ad$agestrat==">10",],vcv=T)

# --------------------------------------
# Estimate adjusted CIRs
# From the overall pooled model and
# From stratified models
# --------------------------------------

# function to get CIR Estimates and CIs from simple stratified models
getCIR <- function(x) {
	# x : log-linear model object returned from coeftest (class=coeftest)
	# NOTE: assumes exposure of interest is the first covariate and there are no interactions
	est <- exp(x[2,1])
	se  <- x[2,2]	
	lb <- exp(log(est)-1.96*se)
	ub <- exp(log(est)+1.96*se)
	res <- c(est,lb,ub)
	return(res)
}


# bind all of the CIRs together, organized by
# non-point source vs. point source and then age, then overall

cir.all       <- getCIR(overall.fit)
cir.age0to4   <- getCIR(age0to4$fit)
cir.age5to10  <- getCIR(age5to10$fit)
cir.age11plus <- getCIR(age11plus$fit)


# --------------------------------------
# Number of persons at risk and cases
# by strata
# Print results to the log
# --------------------------------------
N.all
N.age0to4
N.age5to10
N.age11plus

# --------------------------------------
# Cumulative Incidence Estimates
# Print results to the log
# --------------------------------------
ci.all
ci.0to4
ci.5to10
ci.11plus


# --------------------------------------
# Cumulative Incidence Ratios
# Print results to the log
# --------------------------------------
cir.all
cir.age0to4
cir.age5to10
cir.age11plus



# --------------------------------------
# save the objects
# --------------------------------------

save(
	overall.fit,all.VC,age0to4,age5to10,age11plus,
	N.all,N.age0to4,N.age5to10,N.age11plus,	
	ci.all,ci.0to4,ci.5to10,ci.11plus,
	cir.all,cir.age0to4,cir.age5to10,cir.age11plus,
	file="~/dropbox/13beaches/aim1-results/rawoutput/aim1-sens-enteroQPCR-470cce-regs-noswim-negcontrol.Rdata"
	)









