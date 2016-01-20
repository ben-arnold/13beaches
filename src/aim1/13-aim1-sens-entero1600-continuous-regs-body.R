



# --------------------------------------
# 11-aim1-sens-entero1600-continuous-regs-body.R
# ben arnold (benarnold@berkeley.edu)
#
# description:
# estimate the association between water
# quality indicator concentrations and
# the risk of Diarrhea among swimmers
# for the 13 beaches study
#
# Analyses are conducted for EPA 1600
# Among Swimmers with Body Immmersion
#
# The exposure is treated as continuous (log10 CFU per 100ml)
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

# drop individuals with no water quality information
table(ad$nowq)
ad <- subset(ad,nowq==0)
dim(ad)

# subset to non-missing exposure categories
# to make the robust CI calcs work
table(ad$bodycontact)
ad <- subset(ad,ad$bodycontact=="Yes")
	dim(ad)
ad <- subset(ad,is.na(ad$entero1600)==FALSE)
	dim(ad)

# --------------------------------------
# Adjusted CIR estimates
# --------------------------------------

# --------------------------------------
# Pooled estimates
# (can't use the mpreg fn because we 
# need the actual glm returned object 
# for the LR tests)

# all beaches
all.fit <- glm(diarrheaci10~entero1600 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,family=poisson(link="log"),data=ad)
	all.VC <- cl(ad,fm=all.fit,cluster=ad$hhid)
	overall.fit <- coeftest(all.fit, all.VC)
	summary(all.fit)
	overall.fit

# Interaction model with fresh v. marine beaches
mf.ref <- glm(diarrheaci10~entero1600 +marine +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,family=poisson(link="log"),data=ad)
mf.fit <- glm(diarrheaci10~entero1600*marine +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,family=poisson(link="log"),data=ad)
	# commented out robust SE calcs b/c not used
	# mf.VC <- cl(ad,fm=mf.fit,cluster=ad$hhid)
	# mf.body <- coeftest(mf.fit, mf.VC) 
	summary(mf.fit)
	lrtest(mf.ref,mf.fit)
	

# Interaction model with point v. non-point source beaches
ps.ref <- glm(diarrheaci10~entero1600 +pointsource +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,family=poisson(link="log"),data=ad)
ps.fit <- glm(diarrheaci10~entero1600*pointsource +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,family=poisson(link="log"),data=ad)
	# commented out robust SE calcs b/c not used
	# ps.VC <- cl(ad,fm=ps.fit,cluster=ad$hhid)
	# ps.body <- coeftest(ps.fit, ps.VC) 
	summary(ps.fit)
	lrtest(ps.ref,ps.fit)

# --------------------------------------
# Age-stratified estimates and LR tests of
# interaction

# reduced models for LR tests of indicator x age interactions
agestrat.ref <- glm(diarrheaci10~entero1600 +agestrat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,family=poisson(link="log"),data=ad)

# Pooled estimate (Age 0-4, 5-10, >10), Body Immersion
agestrat.fit <- glm(diarrheaci10~ entero1600*agestrat +female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,family=poisson(link="log"),data=ad)
	# commented out robust SE calcs b/c not used
	# agestrat.VC <- cl(ad,fm=agestrat.fit,cluster=ad$hhid)
	# agestrat.body <- coeftest(agestrat.fit, agestrat.VC) 
	lrtest(agestrat.ref,agestrat.fit)


# --------------------------------------
# Stratified Models 
# based on tests of interaction (above)
# stratify by age
# --------------------------------------

# Pooled Non-point source and point-source
nps <- mpreg(diarrheaci10~entero1600 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ad[ad$pointsource=="No",],vcv=T)
ps <- mpreg(diarrheaci10~entero1600 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ad[ad$pointsource=="Yes",],vcv=T)

# Age-stratified results
age0to4   <- mpreg(diarrheaci10~entero1600 +female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ad[ad$agestrat=="(0, 4]",],vcv=T)
age5to10  <- mpreg(diarrheaci10~entero1600 +female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ad[ad$agestrat=="(4, 10]",],vcv=T)
age11plus <- mpreg(diarrheaci10~entero1600 +female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ad[ad$agestrat==">10",],vcv=T)

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
# non-point source vs. point source and then age

cir.ps <- rbind(getCIR(ps$fit),getCIR(nps$fit),getCIR(overall.fit) )

cir.age<- rbind(getCIR(age0to4$fit),getCIR(age5to10$fit),getCIR(age11plus$fit),getCIR(overall.fit) )


colnames(cir.ps) <- colnames(cir.age) <- c("CIR","CIRlb","CIRub")
rownames(cir.ps) <- c("Point source","Non-point source","Overall")
rownames(cir.age) <- c("Age 0 to 4","Age 5 to 10", "Age >10", "Overall")



# --------------------------------------
# Cumulative Incidence Ratios
# Print results to the log
# --------------------------------------
cir.age
cir.ps


# --------------------------------------
# calculate point-wise marginal predicted 
# probabilities and SEs for plotting
# the dose-response functions
#
# do this using a bootstrap approach
# since there is not yet an easy 
# implementation of delta-method SEs for
# R.  The method is described in these
# articles:
# Muller, C. J. & MacLehose, R. F. Estimating predicted probabilities from logistic regression: different methods correspond to different target populations. Int J Epidemiol, 2014, 43, 962-970
#
# Ahern, J.; Hubbard, A. & Galea, S. Estimating the effects of potential public health interventions on population disease burden: a step-by-step illustration of causal inference methods. Am J Epidemiol, 2009, 169, 1140-1147 
#
# refer to the aim 1 base functions for the 
# underlying bootstrap and estimation methods
#
# --------------------------------------

set.seed(4327234)

# Total Population
pYcurve.all <- boot.pY(diarrheaci10~entero1600 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,data=ad,nameX="entero1600",ID=ad$hhid,strata=ad$beach,iter=1000)

# Age stratified
pYcurve.age0to4 <- boot.pY(diarrheaci10~entero1600 +female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,data=ad[ad$agestrat=="(0, 4]",],nameX="entero1600",ID=ad$hhid[ad$agestrat=="(0, 4]"],strata=ad$beach[ad$agestrat=="(0, 4]"],iter=1000)

pYcurve.age5to10 <- boot.pY(diarrheaci10~entero1600 +female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,data=ad[ad$agestrat=="(4, 10]",],nameX="entero1600",ID=ad$hhid[ad$agestrat=="(4, 10]"],strata=ad$beach[ad$agestrat=="(4, 10]"],iter=1000)

pYcurve.age11plus <- boot.pY(diarrheaci10~entero1600 +female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,data=ad[ad$agestrat==">10",],nameX="entero1600",ID=ad$hhid[ad$agestrat==">10"],strata=ad$beach[ad$agestrat==">10"],iter=1000)

# Point source v. non-point source stratified

pYcurve.ps <- boot.pY(diarrheaci10~entero1600 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,data=ad[ad$pointsource=="Yes",],nameX="entero1600",ID=ad$hhid[ad$pointsource=="Yes"],strata=ad$beach[ad$pointsource=="Yes"],iter=1000)

pYcurve.nps <- boot.pY(diarrheaci10~entero1600 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,data=ad[ad$pointsource=="No",],nameX="entero1600",ID=ad$hhid[ad$pointsource=="No"],strata=ad$beach[ad$pointsource=="No"],iter=1000)




# --------------------------------------
# save the objects
# --------------------------------------

save(
	overall.fit,all.VC,
	age0to4,age5to10,age11plus,
	ps, nps,
	
	cir.age,cir.ps, 
	
	pYcurve.all, pYcurve.age0to4, pYcurve.age5to10, pYcurve.age11plus, pYcurve.ps, pYcurve.nps,
	
	file="~/dropbox/13beaches/aim1-results/rawoutput/aim1-entero1600-continuous-regs-body.Rdata"
	)









