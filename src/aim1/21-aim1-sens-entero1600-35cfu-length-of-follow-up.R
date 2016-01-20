



# --------------------------------------
# 11-aim1-sens-entero1600-35cfu-length-of-follow-up.R
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
# The exposure categories are <=35 and >35 CFU/100ml
#
# Here, we conduct a sensitivity analysis
# to look at the effect of follow-up period
# choice, by including between 1 to 10 days
# of follow-up
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
# Convenience functions used below
# --------------------------------------
# function to calculate number of cases by a categorical exposure
calcNs <- function(outcome,exposurecat) {
  n <- tapply(outcome,exposurecat,sum)
  N <- tapply(outcome,exposurecat,function(x) length(x))
  cbind(n,N)
}


# function to get Estimates and CIs from a linear combination of regression coefficients
# for >=35 CFU vs < 35 CFU
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


# --------------------------------------
# Sensitivity analysis function --
# facilitates repeating the exact analysis
# with just a different outcome, which
# corresponds to different lengths of
# follow-up
# --------------------------------------

sens.entero35 <- function(ad,Y="diarrheaci3") {
	# ad : dataset for analysis
	# Y  : outcome of interest
	
	# create generic outcome variable
	ad$Y <- ad[,Y]

	# --------------------------------------
	# Summarize 
	# the number of swimmers and number
	# of cumulative incident episodes by
	# Entero 1600 category
	# --------------------------------------
	cat("\n---------------------------------\n")
	cat("Sensitivity analysis for Entero >35\n")
	cat("Outcome: ",Y)
	cat("\n---------------------------------\n")
	
	# overall
	N.all <- calcNs(ad$Y,ad$entero35)
	# age stratified
	N.age0to4   <- calcNs(ad$Y[ad$agestrat=="(0, 4]"],ad$entero35[ad$agestrat=="(0, 4]"] )
	N.age5to10  <- calcNs(ad$Y[ad$agestrat=="(4, 10]"],ad$entero35[ad$agestrat=="(4, 10]"] )
	N.age11plus <- calcNs(ad$Y[ad$agestrat==">10"],ad$entero35[ad$agestrat==">10"] )
	# point source vs. non-point source stratified
	N.nps <- calcNs(ad$Y[ad$pointsource=="No"], ad$entero35[ad$pointsource=="No"])
	N.ps  <- calcNs(ad$Y[ad$pointsource=="Yes"],ad$entero35[ad$pointsource=="Yes"])
	
	# --------------------------------------
	# Calculate unadjusted 
	# Cumulative Incidence Rates
	# and robust 95% CIs using a saturated, 
	# unadjusted model
	# (just using the models to get the robust
	# SEs that account for HH clustering)
	# --------------------------------------
	
	cat("\n----------------------------------------\n")
	cat("Cumlative Indicence estimates\n")
	cat("from a saturated model (to get robust SEs)")
	cat("\n----------------------------------------\n")

	## Overall Estimate and Age-stratified estimates
	
	# All Ages
	allci <- mpreg(Y~entero35,dat=ad,vcv=T)
	
	# Ages 0 - 4
	age0to4ci <- mpreg(Y~entero35,dat=ad[ad$agestrat=="(0, 4]",],vcv=T)
	
	# Ages 5 - 10
	age5to10ci <- mpreg(Y~entero35,dat=ad[ad$agestrat=="(4, 10]",],vcv=T)
	
	# Ages >10
	age11plusci <- mpreg(Y~entero35,dat=ad[ad$agestrat==">10",],vcv=T)
	
	## Stratified Estimates by non-point source and point source beach conditions
	
	# Non-point source - All Ages
	npsci <- mpreg(Y~entero35,dat=ad[ad$pointsource=="No",],vcv=T)
	
	# Point source	- All Ages
	psci <- mpreg(Y~entero35,dat=ad[ad$pointsource=="Yes",],vcv=T)
	
	ci.all       <- getCI(allci$fit,allci$vcovCL)
	ci.age0to4   <- getCI(age0to4ci$fit,age0to4ci$vcovCL)
	ci.age5to10  <- getCI(age5to10ci$fit,age5to10ci$vcovCL)
	ci.age11plus <- getCI(age11plusci$fit,age11plusci$vcovCL)
	ci.nps       <- getCI(npsci$fit,npsci$vcovCL)
	ci.ps        <- getCI(psci$fit,psci$vcovCL)
	
	# --------------------------------------
	# esimate risk of Diarrhea associated
	# with exposure to EPA 1600 above/below 35 CFU/100ml
	# Body Immersion
	# All ages
	# --------------------------------------
	cat("\n---------------------------------\n")
	cat("Adjusted Regression Estimates")
	cat("\n---------------------------------\n")
	
	# --------------------------------------
	# Pooled estimate
	
	cat("\n---------------------------------\n")
	cat("Overall")
	cat("\n---------------------------------\n")
	# all beaches
	overall.fit <- mpreg(Y~entero35 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ad,vcv=T)
	
	
	# --------------------------------------
	# Stratified Models 
	# based on tests of interaction (above)
	# stratify the results by non-point and
	# point source conditions
	# --------------------------------------
	
	cat("\n---------------------------------\n")
	cat("All ages, non-point source")
	cat("\n---------------------------------\n")
	# Pooled Non-point source and point-source
	nps <- mpreg(Y~entero35 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ad[ad$pointsource=="No",],vcv=T)
	cat("\n---------------------------------\n")
	cat("All ages, point source")
	cat("\n---------------------------------\n")
	ps <- mpreg(Y~entero35 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ad[ad$pointsource=="Yes",],vcv=T)
	
	# Age-stratified results
	cat("\n---------------------------------\n")
	cat("Ages 0 to 4")
	cat("\n---------------------------------\n")
	age0to4   <- mpreg(Y~entero35 +female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ad[ad$agestrat=="(0, 4]",],vcv=T)
	cat("\n---------------------------------\n")
	cat("Ages 5 to 10")
	cat("\n---------------------------------\n")
	age5to10  <- mpreg(Y~entero35 +female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ad[ad$agestrat=="(4, 10]",],vcv=T)
	cat("\n---------------------------------\n")
	cat("Ages >10")
	cat("\n---------------------------------\n")
	age11plus <- mpreg(Y~entero35 +female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ad[ad$agestrat==">10",],vcv=T)
	
	
	# --------------------------------------
	# Estimate adjusted CIRs
	# From the overall pooled model and
	# From stratified models
	# --------------------------------------
	
	cir.all <- getCIR(overall.fit$fit) 
	
	cir.age0to4   <- getCIR(age0to4$fit)
	cir.age5to10  <- getCIR(age5to10$fit)
	cir.age11plus <- getCIR(age11plus$fit)
	
	cir.nps <- getCIR(nps$fit)
	cir.ps  <- getCIR(ps$fit)
	
	
	# --------------------------------------
	# Number of persons at risk and cases
	# by strata
	# Print results to the log
	# --------------------------------------
	cat("\n---------------------------------\n")
	cat("Number at Risk and Number of Cases")
	cat("\n---------------------------------\n")
	cat("\n Overall\n")
	print(N.all)
	cat("\n Ages 0 to 4\n")
	print(N.age0to4)
	cat("\n Ages 5 to 10\n")
	print(N.age5to10)
	cat("\n Ages >10\n")
	print(N.age11plus)
	cat("\n All Ages, non-point source\n")
	print(N.nps)
	cat("\n All Ages, point-source\n")
	print(N.ps)
	
	# --------------------------------------
	# Cumulative Incidence Estimates
	# Print results to the log
	# --------------------------------------
	cat("\n---------------------------------\n")
	cat("Cumlative Indicence")
	cat("\n---------------------------------\n")
	cat("\n Overall\n")
	print(ci.all)
	cat("\n Ages 0 to 4\n")
	print(ci.age0to4)
	cat("\n Ages 5 to 10\n")
	print(ci.age5to10)
	cat("\n Ages >10\n")
	print(ci.age11plus)
	cat("\n All Ages, non-point source\n")
	print(ci.nps)
	cat("\n All Ages, point-source\n")
	print(ci.ps)
	
	# --------------------------------------
	# Cumulative Incidence Ratios
	# Print results to the log
	# --------------------------------------
	cat("\n---------------------------------\n")
	cat("Cumlative Indicence Ratios (CIR)")
	cat("\n---------------------------------\n")
	cat("\n Overall\n")
	print(cir.all)
	cat("\n Ages 0 to 4\n")
	print(cir.age0to4)
	cat("\n Ages 5 to 10\n")
	print(cir.age5to10)
	cat("\n Ages >10\n")
	print(cir.age11plus)
	cat("\n All Ages, non-point source\n")
	print(cir.nps)
	cat("\n All Ages, point-source\n")
	print(cir.ps)
	
	# --------------------------------------
	# Save results
	# --------------------------------------
	list(overall.fit=overall.fit,age0to4=age0to4,age5to10=age5to10,age11plus=age11plus,
		nps=nps,ps=ps,
		
		N.all=N.all,N.age0to4=N.age0to4,N.age5to10=N.age5to10,N.age11plus=N.age11plus,
		N.nps=N.nps,N.ps=N.ps,
		
		ci.all=ci.all,ci.age0to4=ci.age0to4,ci.age5to10=ci.age5to10, ci.age11plus=ci.age11plus,
		ci.nps=ci.nps,ci.ps=ci.ps,
		
		cir.all=cir.all,cir.age0to4=cir.age0to4,cir.age5to10=cir.age5to10, cir.age11plus=cir.age11plus,
		cir.nps=cir.nps,cir.ps=cir.ps)

}

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
ad <- subset(ad,is.na(ad$entero35)==FALSE)
	dim(ad)



# --------------------------------------
# run the analyses for different 
# lengths of follow-up
# --------------------------------------

d1  <- sens.entero35(ad=ad,Y="diarrheaci1")
d2  <- sens.entero35(ad=ad,Y="diarrheaci2")
d3  <- sens.entero35(ad=ad,Y="diarrheaci3")
d4  <- sens.entero35(ad=ad,Y="diarrheaci4")
d5  <- sens.entero35(ad=ad,Y="diarrheaci5")
d6  <- sens.entero35(ad=ad,Y="diarrheaci6")
d7  <- sens.entero35(ad=ad,Y="diarrheaci7")
d8  <- sens.entero35(ad=ad,Y="diarrheaci8")
d9  <- sens.entero35(ad=ad,Y="diarrheaci9")
d10 <- sens.entero35(ad=ad,Y="diarrheaci10")


# --------------------------------------
# collate the objects across follow-up
# periods (just for some of the objects)
# --------------------------------------
ds <- list(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10)

N.all <- sapply(ds,function(x) x$N.all)
N.ps  <- sapply(ds,function(x) x$N.ps)
N.nps <- sapply(ds,function(x) x$N.nps)

cir.all <- sapply(ds,function(x) x$cir.all)
cir.ps  <- sapply(ds,function(x) x$cir.ps)
cir.nps <- sapply(ds,function(x) x$cir.nps)

# --------------------------------------
# save the objects
# --------------------------------------
rm(ad)
save.image(
	file="~/dropbox/13beaches/aim1-results/rawoutput/aim1-sens-entero1600-35cfu-length-of-follow-up.Rdata"
	)









