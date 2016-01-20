



# --------------------------------------
# 2-aim1-entero1600-Quartile-regs-swall.R
# ben arnold (benarnold@berkeley.edu)
#
# description:
# estimate the association between water
# quality indicator concentrations and
# the risk of Diarrhea among swimmers
# for the 13 beaches study
#
# Analyses are conducted for EPA 1600
# Among beachgoers who swallowed water
#
# NOTE: this script is identical to
# 2-aim1-entero1600-Quartile-regs-body.R
#
# It just subsets the data to individuals
# who swallowed water rather than swimmers
# with body immersion
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
table(ad$swallwater)
ad <- subset(ad,ad$swallwater=="Yes")
	dim(ad)
table(ad$qentero1600)
ad <- subset(ad,is.na(ad$qentero1600)==FALSE)
	dim(ad)


# --------------------------------------
# Calculate the actual Ns for each cell
# and store them for plotting and tables
# --------------------------------------

# --------------------------------------
# Summarize 
# the number of swimmers and number
# of cumulative incident episodes by
# quartile of Entero 1600
# --------------------------------------
calcNs <- function(outcome,exposurecat) {
  n <- tapply(outcome,exposurecat,sum)
  N <- tapply(outcome,exposurecat,function(x) length(x))
  cbind(n,N)
}

# overall
N.all <- calcNs(ad$diarrheaci10,ad$qentero1600)
# age stratified
N.age0to4   <- calcNs(ad$diarrheaci10[ad$agestrat=="(0, 4]"],ad$qentero1600[ad$agestrat=="(0, 4]"] )
N.age5to10  <- calcNs(ad$diarrheaci10[ad$agestrat=="(4, 10]"],ad$qentero1600[ad$agestrat=="(4, 10]"] )
N.age11plus <- calcNs(ad$diarrheaci10[ad$agestrat==">10"],ad$qentero1600[ad$agestrat==">10"] )
# point source vs. non-point source stratified
N.nps <- calcNs(ad$diarrheaci10[ad$pointsource=="No"], ad$qentero1600[ad$pointsource=="No"])
N.ps  <- calcNs(ad$diarrheaci10[ad$pointsource=="Yes"],ad$qentero1600[ad$pointsource=="Yes"])
# age and point source stratified
N.nps0to4   <- calcNs(ad$diarrheaci10[ad$agestrat=="(0, 4]" & ad$pointsource=="No"],ad$qentero1600[ad$agestrat=="(0, 4]" & ad$pointsource=="No"] )
N.nps5to10  <- calcNs(ad$diarrheaci10[ad$agestrat=="(4, 10]" & ad$pointsource=="No"],ad$qentero1600[ad$agestrat=="(4, 10]" & ad$pointsource=="No"] )
N.nps11plus <- calcNs(ad$diarrheaci10[ad$agestrat==">10" & ad$pointsource=="No"],ad$qentero1600[ad$agestrat==">10" & ad$pointsource=="No"] )

N.ps0to4   <- calcNs(ad$diarrheaci10[ad$agestrat=="(0, 4]" & ad$pointsource=="Yes"],ad$qentero1600[ad$agestrat=="(0, 4]" & ad$pointsource=="Yes"] )
N.ps5to10  <- calcNs(ad$diarrheaci10[ad$agestrat=="(4, 10]" & ad$pointsource=="Yes"],ad$qentero1600[ad$agestrat=="(4, 10]" & ad$pointsource=="Yes"] )
N.ps11plus <- calcNs(ad$diarrheaci10[ad$agestrat==">10" & ad$pointsource=="Yes"],ad$qentero1600[ad$agestrat==">10" & ad$pointsource=="Yes"] )


# --------------------------------------
# Calculate unadjusted 
# Cumulative Incidence Rates
# and robust 95% CIs using a saturated, 
# unadjusted model
# (just using the models to get the robust
# SEs that account for HH clustering)
# --------------------------------------

# fit the saturated models
# for now, only doing this for the overall results
# (not stratified by point/nonpoint because the
#  sample sizes are really small with the double stratification)

# All Ages
allci <- mpreg(diarrheaci10~qentero1600,dat=ad,vcv=T)

# Ages 0 - 4
age0to4ci <- mpreg(diarrheaci10~qentero1600,dat=ad[ad$agestrat=="(0, 4]",],vcv=T)

# Ages 5 - 10
age5to10ci <- mpreg(diarrheaci10~qentero1600,dat=ad[ad$agestrat=="(4, 10]",],vcv=T)

# Ages >10
age11plusci <- mpreg(diarrheaci10~qentero1600,dat=ad[ad$agestrat==">10",],vcv=T)

# Point source
psci <- mpreg(diarrheaci10~qentero1600,dat=ad[ad$pointsource=="Yes",],vcv=T)

# Non-point source
npsci <- mpreg(diarrheaci10~qentero1600,dat=ad[ad$pointsource=="No",],vcv=T)

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
	# fit : an object returned from coeftest w/ 4 parameters corresponding to Q1-Q4 cumulative incidence
	# vcv : variance-covariance matrix for coefficients in the fit object
	lc1 <- lc2 <- lc3 <- lc4 <- rep(0,nrow(fit))
	lc1[c(1)  ] <- 1
	lc2[c(1,2)] <- 1
	lc3[c(1,3)] <- 1
	lc4[c(1,4)] <- 1
	lcs <- list(lc1,lc2,lc3,lc4)
	res <- t(sapply(lcs,lccalc,x=fit,vcv=vcv))
	colnames(res) <- c("CI","CIlb","CIub")
	rownames(res) <- paste("Q",1:4,sep="")
	return(res)
}



ci.all    <- getCI(allci$fit,allci$vcovCL)
ci.0to4   <- getCI(age0to4ci$fit,age0to4ci$vcovCL)
ci.5to10  <- getCI(age5to10ci$fit,age5to10ci$vcovCL)
ci.11plus <- getCI(age11plusci$fit,age11plusci$vcovCL)

ci.ps  <- getCI(psci$fit,psci$vcovCL)
ci.nps <- getCI(npsci$fit,npsci$vcovCL)

# --------------------------------------
# Test for interaction (effect modification)
# by environmental conditions at
# Avalon and Doheny beaches
# --------------------------------------

# Avalon
# Above versus below median groundwater flow
# See Yau et al. 2014 Wat Res for details
av.h.noint <-glm(diarrheaci10~ qentero1600+groundwater +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,family=poisson(link="log"),data=ad[ad$beach=="Avalon",])
av.h <-glm(diarrheaci10~ qentero1600*groundwater +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,family=poisson(link="log"),data=ad[ad$beach=="Avalon",])
	lrtest(av.h.noint,av.h)

# Doheny
# Sand berm open vs. closed
# See Colford et al. 2012 Wat Res for details
dh.h.noint <-glm(diarrheaci10~qentero1600+berm +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,family=poisson(link="log"),data=ad[ad$beach=="Doheny",])
dh.h <-glm(diarrheaci10~ qentero1600*berm +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,family=poisson(link="log"),data=ad[ad$beach=="Doheny",])
	lrtest(dh.h.noint,dh.h)

# --------------------------------------
# Beach Specific Stratified Estimates
# --------------------------------------

# --------------------------------------
# Freshwater beaches

# Huntington
hufit <- mpreg(diarrheaci10~qentero1600 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ad[ad$beach=="Huntington",],vcv=T)

# Silver
sifit <- mpreg(diarrheaci10~qentero1600 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ad[ad$beach=="Silver",],vcv=T)

# Washington Park
wpfit <- mpreg(diarrheaci10~qentero1600 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ad[ad$beach=="Washington Park",],vcv=T)

# West
wefit <- mpreg(diarrheaci10~qentero1600 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ad[ad$beach=="West",],vcv=T)

# --------------------------------------
# Marine beaches

# Avalon
avfit <- mpreg(diarrheaci10~qentero1600*groundwater +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ad[ad$beach=="Avalon",],vcv=T)

# Boqueron
bofit <- mpreg(diarrheaci10~qentero1600 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ad[ad$beach=="Boqueron",],vcv=T)

# Doheny
dhfit <- mpreg(diarrheaci10~qentero1600*berm +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ad[ad$beach=="Doheny",],vcv=T)

# Edgewater
edfit <- mpreg(diarrheaci10~qentero1600 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ad[ad$beach=="Edgewater",],vcv=T)

# Fairhope
fafit <- mpreg(diarrheaci10~qentero1600 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ad[ad$beach=="Fairhope",],vcv=T)

# Goddard
gdfit <- mpreg(diarrheaci10~qentero1600 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ad[ad$beach=="Goddard",],vcv=T)

# Malibu
mafit <- mpreg(diarrheaci10~qentero1600 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ad[ad$beach=="Malibu",],vcv=T)

# Mission Bay
mbfit <- mpreg(diarrheaci10~qentero1600 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ad[ad$beach=="Mission Bay",],vcv=T)

# Surfside
sufit <- mpreg(diarrheaci10~qentero1600 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ad[ad$beach=="Surfside",],vcv=T)

# --------------------------------------
# Pooled estimates
# (can't use the mpreg fn because we 
# need the actual glm returned object 
# for the LR tests)

# all beaches
all.fit <- glm(diarrheaci10~qentero1600 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,family=poisson(link="log"),data=ad)
	all.VC <- cl(ad,fm=all.fit,cluster=ad$hhid)
	overall.fit <- coeftest(all.fit, all.VC)
	summary(all.fit)
	overall.fit

# Interaction model with fresh v. marine beaches
mf.ref <- glm(diarrheaci10~qentero1600+marine +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,family=poisson(link="log"),data=ad)
mf.fit <- glm(diarrheaci10~qentero1600*marine +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,family=poisson(link="log"),data=ad)
	# commented out robust SE calcs b/c not used
	# mf.VC <- cl(ad,fm=mf.fit,cluster=ad$hhid)
	# mf.head <- coeftest(mf.fit, mf.VC)
	summary(mf.fit) 
	lrtest(mf.ref,mf.fit)
	

# Interaction model with point v. non-point source beaches
ps.ref <- glm(diarrheaci10~qentero1600+pointsource +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,family=poisson(link="log"),data=ad)
ps.fit <- glm(diarrheaci10~qentero1600*pointsource+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,family=poisson(link="log"),data=ad)
	# commented out robust SE calcs b/c not used
	# ps.VC <- cl(ad,fm=ps.fit,cluster=ad$hhid)
	# ps.head <- coeftest(ps.fit, ps.VC)
	summary(ps.fit)
	lrtest(ps.ref,ps.fit)
	

# --------------------------------------
# Age-stratified estimates and LR tests of
# interaction

# reduced models for LR tests of indicator x age interactions
agestrat.ref <- glm(diarrheaci10~qentero1600 +agestrat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,family=poisson(link="log"),data=ad)

# Pooled estimate (Age 0-4, 5-10, >10)
agestrat.fit <- glm(diarrheaci10~ qentero1600*agestrat +female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,family=poisson(link="log"),data=ad)
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

# Non-point source and point-source
nps <- mpreg(diarrheaci10~qentero1600 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ad[ad$pointsource=="No",],vcv=T)
ps <- mpreg(diarrheaci10~qentero1600 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ad[ad$pointsource=="Yes",],vcv=T)

# Age-stratified results
age0to4   <- mpreg(diarrheaci10~qentero1600 +female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ad[ad$agestrat=="(0, 4]",],vcv=T)
age5to10  <- mpreg(diarrheaci10~qentero1600 +female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ad[ad$agestrat=="(4, 10]",],vcv=T)
age11plus <- mpreg(diarrheaci10~qentero1600 +female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ad[ad$agestrat==">10",],vcv=T)

# Non-point source estimates, stratified by age group
nps0to4 <- mpreg(diarrheaci10~qentero1600 +female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ad[ad$pointsource=="No" & ad$agestrat=="(0, 4]",],vcv=T)

nps5to10 <- mpreg(diarrheaci10~qentero1600 +female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ad[ad$pointsource=="No" & ad$agestrat=="(4, 10]",],vcv=T)

nps11plus <- mpreg(diarrheaci10~qentero1600 +female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ad[ad$pointsource=="No" & ad$agestrat==">10",],vcv=T)

# Point-source estimates by age group
ps0to4 <- mpreg(diarrheaci10~qentero1600 +female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ad[ad$pointsource=="Yes" & ad$agestrat=="(0, 4]",],vcv=T)

ps5to10 <- mpreg(diarrheaci10~qentero1600 +female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ad[ad$pointsource=="Yes" & ad$agestrat=="(4, 10]",],vcv=T)

ps11plus <- mpreg(diarrheaci10~qentero1600 +female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ad[ad$pointsource=="Yes" & ad$agestrat==">10",],vcv=T)


# --------------------------------------
# Calculate adjusted CIRs and 95% CIs
# From the regression coefficients
# --------------------------------------

# function to get Estimates and CIs for quartiles of exposure
# NOTE: assumes that the quartiles exposure variables are #s 2-4
CIRquartiles <- function(x,vcv) {
	# x : log-linear model object returned from coeftest (class=coeftest)
	# vcv : variance-covariance matrix of coefficients for robust SEs
	lc2 <- c(0,1,rep(0,nrow(x)-2))
	lc3 <- c(0,0,1,rep(0,nrow(x)-3))
	lc4 <- c(0,0,0,1,rep(0,nrow(x)-4))	
	estci <- function(lc,x,vcv) {
		est <- exp(t(lc)%*%x[,1])
		se  <- sqrt( t(lc)%*%vcv%*%lc )
		lb <- exp(log(est)-1.96*se)
		ub <- exp(log(est)+1.96*se)
		res <- c(est,lb,ub)
		return(res)
	}
	cirs <- t(sapply(list(lc2,lc3,lc4),estci,x=x,vcv=vcv))
	colnames(cirs) <- c("CIR","CIRlb","CIRub")
	rownames(cirs) <- paste("Entero1600-Q",2:4,sep="")
	return(cirs)
}

cir.all <- CIRquartiles(overall.fit,all.VC)

cir.nps <- CIRquartiles(nps$fit,nps$vcovCL)
cir.ps  <- CIRquartiles(ps$fit,  ps$vcovCL)

cir.age0to4   <- CIRquartiles(age0to4$fit,  age0to4$vcovCL)
cir.age5to10  <- CIRquartiles(age5to10$fit, age5to10$vcovCL)
cir.age11plus <- CIRquartiles(age11plus$fit,age11plus$vcovCL)

cir.nps0to4   <- CIRquartiles(nps0to4$fit,  nps0to4$vcovCL)
cir.nps5to10  <- CIRquartiles(nps5to10$fit, nps5to10$vcovCL)
cir.nps11plus <- CIRquartiles(nps11plus$fit,nps11plus$vcovCL)

cir.ps0to4   <- CIRquartiles(ps0to4$fit,  ps0to4$vcovCL)
cir.ps5to10  <- CIRquartiles(ps5to10$fit, ps5to10$vcovCL)
cir.ps11plus <- CIRquartiles(ps11plus$fit,ps11plus$vcovCL)


# --------------------------------------
# Print N counts to the log file
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
# Print Cumulative Incidence estimates
# to the log file
# --------------------------------------
ci.all
ci.0to4
ci.5to10
ci.11plus
ci.ps
ci.nps


# --------------------------------------
# Print CIR estimates to the log file
# --------------------------------------

cir.all
cir.nps
cir.ps
cir.age0to4
cir.age5to10
cir.age11plus
cir.nps0to4
cir.nps5to10
cir.nps11plus
cir.ps0to4
cir.ps5to10
cir.ps11plus


# --------------------------------------
# save the results
# exclude glm objects and data frames
# (they are really large)
# --------------------------------------
save(
	overall.fit,all.VC,age0to4,age5to10,age11plus,
	nps,nps0to4,nps5to10,nps11plus,
	ps,ps0to4,ps5to10,ps11plus,
	
	hufit,sifit,wpfit,wefit,avfit,bofit,dhfit,edfit,fafit,gdfit,mafit,mbfit,sufit,
	
	N.all,N.age0to4,N.age5to10,N.age11plus,
	N.nps,N.nps0to4,N.nps5to10,N.nps11plus,
	N.ps,N.ps0to4,N.ps5to10,N.ps11plus,
	
	ci.all,ci.0to4,ci.5to10,ci.11plus,ci.ps,ci.nps,
	
	cir.all,cir.age0to4,cir.age5to10,cir.age11plus,
	cir.nps,cir.nps0to4,cir.nps5to10,cir.nps11plus,
	cir.ps,cir.ps0to4,cir.ps5to10,cir.ps11plus,

	file="~/dropbox/13beaches/aim1-results/rawoutput/aim1-entero1600-Quartile-regs-head.Rdata"
	)









