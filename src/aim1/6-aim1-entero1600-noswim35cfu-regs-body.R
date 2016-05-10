



# --------------------------------------
# 6-aim1-entero1600-noswim35cfu-regs-body
# ben arnold (benarnold@berkeley.edu)
#
# description:
# estimate the association between water
# quality indicator concentrations and
# the risk of Diarrhea among swimmers
# for the 13 beaches study
#
# Analyses are conducted for EPA 1600
# Among Non-swimmers and Swimmers with Body Immmersion
#
#
# The exposure categories are non-swimmers, <=35 and >35 CFU/100ml
#
# This differs from the analysis in 3-aim1-entero1600-35cfu-regs-body.R
# in that we are using a different reference group. Here, we are considering
# non-swimmers as the reference group. In the other program, swimmers <=35 CFU are the
# reference group.
#
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

# subset to non-missing exposure categories
# to make the robust CI calcs work
table(ad$anycontact,ad$bodycontact)
ad <- subset(ad,ad$anycontact=="No" | ad$bodycontact=="Yes")
	dim(ad)
ad <- subset(ad,is.na(ad$swim35)==FALSE)
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
N.all <- calcNs(ad$diarrheaci10,ad$swim35)
# age stratified
N.age0to4   <- calcNs(ad$diarrheaci10[ad$agestrat=="(0, 4]"],ad$swim35[ad$agestrat=="(0, 4]"] )
N.age5to10  <- calcNs(ad$diarrheaci10[ad$agestrat=="(4, 10]"],ad$swim35[ad$agestrat=="(4, 10]"] )
N.age11plus <- calcNs(ad$diarrheaci10[ad$agestrat==">10"],ad$swim35[ad$agestrat==">10"] )
# point source vs. non-point source stratified
N.nps <- calcNs(ad$diarrheaci10[ad$pointsource=="No"], ad$swim35[ad$pointsource=="No"])
N.ps  <- calcNs(ad$diarrheaci10[ad$pointsource=="Yes"],ad$swim35[ad$pointsource=="Yes"])
# age and point source stratified
N.nps0to4   <- calcNs(ad$diarrheaci10[ad$agestrat=="(0, 4]" & ad$pointsource=="No"],ad$swim35[ad$agestrat=="(0, 4]" & ad$pointsource=="No"] )
N.nps5to10  <- calcNs(ad$diarrheaci10[ad$agestrat=="(4, 10]" & ad$pointsource=="No"],ad$swim35[ad$agestrat=="(4, 10]" & ad$pointsource=="No"] )
N.nps11plus <- calcNs(ad$diarrheaci10[ad$agestrat==">10" & ad$pointsource=="No"],ad$swim35[ad$agestrat==">10" & ad$pointsource=="No"] )

N.ps0to4   <- calcNs(ad$diarrheaci10[ad$agestrat=="(0, 4]" & ad$pointsource=="Yes"],ad$swim35[ad$agestrat=="(0, 4]" & ad$pointsource=="Yes"] )
N.ps5to10  <- calcNs(ad$diarrheaci10[ad$agestrat=="(4, 10]" & ad$pointsource=="Yes"],ad$swim35[ad$agestrat=="(4, 10]" & ad$pointsource=="Yes"] )
N.ps11plus <- calcNs(ad$diarrheaci10[ad$agestrat==">10" & ad$pointsource=="Yes"],ad$swim35[ad$agestrat==">10" & ad$pointsource=="Yes"] )



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
all.cifit <- glm(diarrheaci10~swim35,family=poisson(link="log"),data=ad)
	allci.VC <- cl(ad,fm=all.cifit,cluster=ad$hhid)
	allci <- coeftest(all.cifit, allci.VC)
	allci

# Ages 0 - 4
age0to4.cifit <- glm(diarrheaci10~swim35,family=poisson(link="log"),data=ad[ad$agestrat=="(0, 4]",])
	age0to4ci.VC <- cl(ad[ad$agestrat=="(0, 4]",],fm=age0to4.cifit,cluster=ad$hhid[ad$agestrat=="(0, 4]"])
	age0to4ci <- coeftest(age0to4.cifit, age0to4ci.VC)
	age0to4ci

# Ages 5 - 10
age5to10.cifit <- glm(diarrheaci10~swim35,family=poisson(link="log"),data=ad[ad$agestrat=="(4, 10]",])
	age5to10ci.VC <- cl(ad[ad$agestrat=="(4, 10]",],fm=age5to10.cifit,cluster=ad$hhid[ad$agestrat=="(4, 10]"])
	age5to10ci <- coeftest(age5to10.cifit, age5to10ci.VC)
	age5to10ci

# Ages >10
age11plus.cifit <- glm(diarrheaci10~swim35,family=poisson(link="log"),data=ad[ad$agestrat==">10",])
	age11plusci.VC <- cl(ad[ad$agestrat==">10",],fm=age11plus.cifit,cluster=ad$hhid[ad$agestrat==">10"])
	age11plusci <- coeftest(age11plus.cifit, age11plusci.VC)
	age11plusci


#### Stratified Estimates by non-point source and point source beach conditions

# Non-point source - All Ages
nps.cifit <- glm(diarrheaci10~swim35,family=poisson(link="log"),data=ad[ad$pointsource=="No",])
	npsci.VC <- cl(ad[ad$pointsource=="No",],fm=nps.cifit,cluster=ad$hhid[ad$pointsource=="No"])
	npsci <- coeftest(nps.cifit, npsci.VC)
	npsci

# Non-point source - Ages 0 - 4
nps0to4.cifit <- glm(diarrheaci10~swim35,family=poisson(link="log"),data=ad[ad$pointsource=="No" & ad$agestrat=="(0, 4]",])
	nps0to4ci.VC <- cl(ad[ad$pointsource=="No" & ad$agestrat=="(0, 4]",],fm=nps0to4.cifit,cluster=ad$hhid[ad$pointsource=="No" & ad$agestrat=="(0, 4]"])
	nps0to4ci <- coeftest(nps0to4.cifit, nps0to4ci.VC)
	nps0to4ci

# Non-point source - Ages 5 - 10
nps5to10.cifit <- glm(diarrheaci10~swim35,family=poisson(link="log"),data=ad[ad$pointsource=="No" & ad$agestrat=="(4, 10]",])
	nps5to10ci.VC <- cl(ad[ad$pointsource=="No" & ad$agestrat=="(4, 10]",],fm=nps5to10.cifit,cluster=ad$hhid[ad$pointsource=="No" & ad$agestrat=="(4, 10]"])
	nps5to10ci <- coeftest(nps5to10.cifit, nps5to10ci.VC)
	nps5to10ci

# Non-point source - Ages >10
nps11plus.cifit <- glm(diarrheaci10~swim35,family=poisson(link="log"),data=ad[ad$pointsource=="No" & ad$agestrat==">10",])
	nps11plusci.VC <- cl(ad[ad$pointsource=="No" & ad$agestrat==">10",],fm=nps11plus.cifit,cluster=ad$hhid[ad$pointsource=="No" & ad$agestrat==">10"])
	nps11plusci <- coeftest(nps11plus.cifit, nps11plusci.VC)
	nps11plusci


# Point source	- All Ages
ps.cifit <- glm(diarrheaci10~swim35,family=poisson(link="log"),data=ad[ad$pointsource=="Yes",])
	psci.VC <- cl(ad[ad$pointsource=="Yes",],fm=ps.cifit,cluster=ad$hhid[ad$pointsource=="Yes"])
	psci <- coeftest(ps.cifit, psci.VC)
	psci

# Point source - Ages 0 - 4
ps0to4.cifit <- glm(diarrheaci10~swim35,family=poisson(link="log"),data=ad[ad$pointsource=="Yes" & ad$agestrat=="(0, 4]",])
	ps0to4ci.VC <- cl(ad[ad$pointsource=="Yes" & ad$agestrat=="(0, 4]",],fm=ps0to4.cifit,cluster=ad$hhid[ad$pointsource=="Yes" & ad$agestrat=="(0, 4]"])
	ps0to4ci <- coeftest(ps0to4.cifit, ps0to4ci.VC)
	ps0to4ci

# Point source - Ages 5 - 10
ps5to10.cifit <- glm(diarrheaci10~swim35,family=poisson(link="log"),data=ad[ad$pointsource=="Yes" & ad$agestrat=="(4, 10]",])
	ps5to10ci.VC <- cl(ad[ad$pointsource=="Yes" & ad$agestrat=="(4, 10]",],fm=ps5to10.cifit,cluster=ad$hhid[ad$pointsource=="Yes" & ad$agestrat=="(4, 10]"])
	ps5to10ci <- coeftest(ps5to10.cifit, ps5to10ci.VC)
	ps5to10ci

# Point source - Ages >10
ps11plus.cifit <- glm(diarrheaci10~swim35,family=poisson(link="log"),data=ad[ad$pointsource=="Yes" & ad$agestrat==">10",])
	ps11plusci.VC <- cl(ad[ad$pointsource=="Yes" & ad$agestrat==">10",],fm=ps11plus.cifit,cluster=ad$hhid[ad$pointsource=="Yes" & ad$agestrat==">10"])
	ps11plusci <- coeftest(ps11plus.cifit, ps11plusci.VC)
	ps11plusci


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
	# fit : an object returned from coeftest w/ 3 parameters corresponding to non-swimmers/<=35/>35 CFU cumulative incidence
	# vcv : variance-covariance matrix for coefficients in the fit object
	lc1 <- lc2 <- lc3 <- rep(0,nrow(fit))
	lc1[c(1)  ] <- 1
	lc2[c(1,2)] <- 1
	lc3[c(1,3)] <- 1
	lcs <- list(lc1,lc2,lc3)
	res <- t(sapply(lcs,lccalc,x=fit,vcv=vcv))
	colnames(res) <- c("CI","CIlb","CIub")
	rownames(res) <- c("Non-swimmers","Entero<=35cfu","Entero>35cfu")
	return(res)
}


ci.all    <- getCI(allci,allci.VC)
ci.0to4   <- getCI(age0to4ci,age0to4ci.VC)
ci.5to10  <- getCI(age5to10ci,age5to10ci.VC)
ci.11plus <- getCI(age11plusci,age11plusci.VC)

ci.nps       <- getCI(npsci,npsci.VC)
ci.nps0to4   <- getCI(nps0to4ci,nps0to4ci.VC)
ci.nps5to10  <- getCI(nps5to10ci,nps5to10ci.VC)
ci.nps11plus <- getCI(nps11plusci,nps11plusci.VC)

ci.ps       <- getCI(psci,psci.VC)
ci.ps0to4   <- getCI(ps0to4ci,ps0to4ci.VC)
ci.ps5to10  <- getCI(ps5to10ci,ps5to10ci.VC)
ci.ps11plus <- getCI(ps11plusci,nps11plusci.VC)



# --------------------------------------
# esimate risk of Diarrhea associated
# with exposure to EPA 1600 above/below 35 CFU/100ml
# Body Immersion
# All ages
# --------------------------------------



# tests of interaction by environmental conditions for Avalon and Doheny
av.h.noint <-glm(diarrheaci10~ swim35+groundwater +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,family=poisson(link="log"),data=ad[ad$beach=="Avalon",])
av.h <-glm(diarrheaci10~ swim35*groundwater +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,family=poisson(link="log"),data=ad[ad$beach=="Avalon",])
	lrtest(av.h.noint,av.h)

dh.h.noint <-glm(diarrheaci10~ swim35+berm +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,family=poisson(link="log"),data=ad[ad$beach=="Doheny",])
dh.h <-glm(diarrheaci10~ swim35*berm +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,family=poisson(link="log"),data=ad[ad$beach=="Doheny",])
	lrtest(dh.h.noint,dh.h)


# --------------------------------------
# Freshwater beaches

# Huntington
hufit <- mpreg(diarrheaci10~swim35 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ad[ad$beach=="Huntington",],vcv=T)

# Silver
sifit <- mpreg(diarrheaci10~swim35 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ad[ad$beach=="Silver",],vcv=T)

# Washington Park
wpfit <- mpreg(diarrheaci10~swim35 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ad[ad$beach=="Washington Park",],vcv=T)

# West
wefit <- mpreg(diarrheaci10~swim35 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ad[ad$beach=="West",],vcv=T)

# --------------------------------------
# Marine beaches

# Avalon
avfit <- mpreg(diarrheaci10~swim35 +groundwater +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ad[ad$beach=="Avalon",],vcv=T)

# Boqueron
bofit <- mpreg(diarrheaci10~swim35 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ad[ad$beach=="Boqueron",],vcv=T)

# Doheny
dhfit <- mpreg(diarrheaci10~swim35 +berm +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ad[ad$beach=="Doheny",],vcv=T)

# Edgewater
edfit <- mpreg(diarrheaci10~swim35 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ad[ad$beach=="Edgewater",],vcv=T)

# Fairhope
fafit <- mpreg(diarrheaci10~swim35 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ad[ad$beach=="Fairhope",],vcv=T)

# Goddard
gdfit <- mpreg(diarrheaci10~swim35 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ad[ad$beach=="Goddard",],vcv=T)

# Malibu
mafit <- mpreg(diarrheaci10~swim35 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ad[ad$beach=="Malibu",],vcv=T)

# Mission Bay
mbfit <- mpreg(diarrheaci10~swim35 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ad[ad$beach=="Mission Bay",],vcv=T)

# Surfside
sufit <- mpreg(diarrheaci10~swim35 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ad[ad$beach=="Surfside",],vcv=T)


# --------------------------------------
# Pooled estimates
# (can't use the mpreg fn because we 
# need the actual glm returned object 
# for the LR tests)

# all beaches
all.fit <- glm(diarrheaci10~swim35 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,family=poisson(link="log"),data=ad)
	all.VC <- cl(ad,fm=all.fit,cluster=ad$hhid)
	overall.fit <- coeftest(all.fit, all.VC)
	summary(all.fit)
	overall.fit

# Interaction model with fresh v. marine beaches
mf.ref <- glm(diarrheaci10~swim35 +marine +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,family=poisson(link="log"),data=ad)
mf.fit <- glm(diarrheaci10~swim35*marine +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,family=poisson(link="log"),data=ad)
	mf.VC <- cl(ad,fm=mf.fit,cluster=ad$hhid)
	mf.body <- coeftest(mf.fit, mf.VC) 
	summary(mf.fit)
	lrtest(mf.ref,mf.fit)
	
	# test for additive interaction
	amf.ref <- glm(diarrheaci10~swim35 +marine +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,family=gaussian,data=ad)
	amf.fit <- glm(diarrheaci10~swim35*marine +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,family=gaussian,data=ad)
	lrtest(amf.ref,amf.fit)
	
# Interaction model with point v. non-point source beaches
ps.ref <- glm(diarrheaci10~swim35 +pointsource +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,family=poisson(link="log"),data=ad)
ps.fit <- glm(diarrheaci10~swim35*pointsource +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,family=poisson(link="log"),data=ad)
	ps.VC <- cl(ad,fm=ps.fit,cluster=ad$hhid)
	ps.body <- coeftest(ps.fit, ps.VC) 
	summary(ps.fit)
	lrtest(ps.ref,ps.fit)
	
	# test for additive interaction
	aps.ref <- glm(diarrheaci10~swim35 +pointsource +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,family=gaussian,data=ad)
	aps.fit <- glm(diarrheaci10~swim35*pointsource +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,family=gaussian,data=ad)
	lrtest(aps.ref,aps.fit)
	
# --------------------------------------
# Age-stratified estimates and LR tests of
# interaction

# reduced models for LR tests of indicator x age interactions
agestrat.ref <- glm(diarrheaci10~swim35 +agestrat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,family=poisson(link="log"),data=ad)

# Pooled estimate (Age 0-4, 5-10, >10), Body Immersion
agestrat.fit <- glm(diarrheaci10~ swim35*agestrat +female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,family=poisson(link="log"),data=ad)
	# commented out robust SE calcs b/c not used
	# agestrat.VC <- cl(ad,fm=agestrat.fit,cluster=ad$hhid)
	# agestrat.body <- coeftest(agestrat.fit, agestrat.VC) 
	lrtest(agestrat.ref,agestrat.fit)
	
	# test for additive interaction
	aage.ref <- glm(diarrheaci10~swim35 +agestrat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,family=gaussian,data=ad)
	aage.fit <- glm(diarrheaci10~ swim35*agestrat +female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,family=gaussian,data=ad)
	lrtest(aage.ref,aage.fit)

# --------------------------------------
# Stratified Models 
# based on tests of interaction (above)
# stratify the results by non-point and
# point source conditions
# --------------------------------------

# Pooled Non-point source and point-source
nps <- mpreg(diarrheaci10~swim35 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ad[ad$pointsource=="No",],vcv=T)
ps <- mpreg(diarrheaci10~swim35 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ad[ad$pointsource=="Yes",],vcv=T)

# Age-stratified results
age0to4   <- mpreg(diarrheaci10~swim35 +female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ad[ad$agestrat=="(0, 4]",],vcv=T)
age5to10  <- mpreg(diarrheaci10~swim35 +female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ad[ad$agestrat=="(4, 10]",],vcv=T)
age11plus <- mpreg(diarrheaci10~swim35 +female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ad[ad$agestrat==">10",],vcv=T)

# Non-point source estimates by age group
nps0to4 <- mpreg(diarrheaci10~swim35 +female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ad[ad$pointsource=="No" & ad$agestrat=="(0, 4]",],vcv=T)

nps5to10 <- mpreg(diarrheaci10~swim35 +female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ad[ad$pointsource=="No" & ad$agestrat=="(4, 10]",],vcv=T)

nps11plus <- mpreg(diarrheaci10~swim35 +female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ad[ad$pointsource=="No" & ad$agestrat==">10",],vcv=T)

# Point-source estimates by age group
ps0to4 <- mpreg(diarrheaci10~swim35 +female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ad[ad$pointsource=="Yes" & ad$agestrat=="(0, 4]",],vcv=T)

ps5to10 <- mpreg(diarrheaci10~swim35 +female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ad[ad$pointsource=="Yes" & ad$agestrat=="(4, 10]",],vcv=T)

ps11plus <- mpreg(diarrheaci10~swim35 +female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ad[ad$pointsource=="Yes" & ad$agestrat==">10",],vcv=T)


# --------------------------------------
# Estimate adjusted CIRs
# From the overall pooled model and
# From stratified models
# --------------------------------------

# function to get CIR Estimates and CIs from simple stratified models
getCIR <- function(x) {
	# x : log-linear model object returned from coeftest (class=coeftest)
	# NOTE: assumes exposure of interest is the first covariate and there are no interactions
	est <- exp(x[2:3,1])
	se  <- x[2:3,2]	
	lb <- exp(log(est)-1.96*se)
	ub <- exp(log(est)+1.96*se)
	res <- c(est[1],lb[1],ub[1],est[2],lb[2],ub[2])
	names(res) <- c("CIRbelow35","CIRbelow35lb","CIRbelow35ub","CIRabove35","CIRabove35lb","CIRabove35ub")
	return(res)
}


# bind all of the CIRs together, organized by
# non-point source vs. point source and then age, then overall

cir.all <- rbind(getCIR(nps$fit),getCIR(ps$fit),getCIR(overall.fit) )

cir.age0to4 <- rbind(getCIR(nps0to4$fit),getCIR(ps0to4$fit),getCIR(age0to4$fit) )

cir.age5to10 <- rbind(getCIR(nps5to10$fit),getCIR(ps5to10$fit),getCIR(age5to10$fit) )

cir.age11plus <- rbind(getCIR(nps11plus$fit),getCIR(ps11plus$fit),getCIR(age11plus$fit) )

rownames(cir.all) <- rownames(cir.age11plus) <- rownames(cir.age5to10) <- rownames(cir.age0to4) <-  c("Non-point source","Point source","Overall")

# --------------------------------------
# Number of persons at risk and cases
# by strata
# Print results to the log
# --------------------------------------
N.all
N.age0to4
N.age5to10
N.age11plus
N.nps
N.ps
N.nps0to4
N.nps5to10
N.nps11plus
N.ps0to4
N.ps5to10
N.ps11plus

# --------------------------------------
# Cumulative Incidence Estimates
# Print results to the log
# --------------------------------------
ci.all
ci.0to4
ci.5to10
ci.11plus

ci.nps
ci.nps0to4
ci.nps5to10
ci.nps11plus

ci.ps
ci.ps0to4
ci.ps5to10
ci.ps11plus


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
	nps,nps0to4,nps5to10,nps11plus,
	ps,ps0to4,ps5to10,ps11plus,
	
	hufit,sifit,wpfit,wefit,avfit,bofit,dhfit,edfit,fafit,gdfit,mafit,mbfit,sufit,
	
	N.all,N.age0to4,N.age5to10,N.age11plus,
	N.nps,N.nps0to4,N.nps5to10,N.nps11plus,
	N.ps,N.ps0to4,N.ps5to10,N.ps11plus,
	
	ci.all,ci.0to4,ci.5to10,ci.11plus,
	ci.nps,ci.nps0to4,ci.nps5to10,ci.nps11plus,
	ci.ps,ci.ps0to4,ci.ps5to10,ci.ps11plus,
	
	cir.all,cir.age0to4,cir.age5to10,cir.age11plus,

	file="~/dropbox/13beaches/aim1-results/rawoutput/aim1-entero1600-noswim35cfu-regs-body.Rdata"
	)









