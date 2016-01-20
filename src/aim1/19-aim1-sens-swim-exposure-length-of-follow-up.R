



# --------------------------------------
# 11-aim1-sensitivity-swim-exposure-length-of-follow-up.R
# ben arnold (benarnold@berkeley.edu)
#
# description:
# estimate the adsociation between 
# water exposure and the risk of diarrhea
# for the 13 beaches study
#
# Here, we conduct a sensitivity analysis
# to look at the effect of follow-up period
# choice, by including between 1 to 10 days
# of follow-up
#
#
# --------------------------------------

# --------------------------------------
# preamble
# --------------------------------------

rm(list=ls())
library(sandwich)
library(lmtest)
library(RColorBrewer)

# source the base functions
source("~/13beaches/src/aim1/0-aim1-base-functions.R")


# --------------------------------------
# Additional Functions specific to 
# Water Exposure Analyses
# --------------------------------------

# --------------------------------------
# Convenience Function to calculate Ns
# For Swim Exposure Analyses
# --------------------------------------

N.swimex <- function(data) {
	# data: data frame that includes the variables referenced below
	#       (designed to be subset to a specific swim exposure category:
	#        non-swimmers, body immersion, head immersion, swallowed water)
	calcNs <- function(x) {
  		n <- sum(x)
  		N <- length(x)
  		cbind(n,N)
		}
	Ns <- calcNs(data$Y)
	return(Ns)
}

# --------------------------------------
# Convenience Function to calculate 
# Cumulative Incidence and robust SEs
# for Swim Exposure Analyses
#
# Note: relies on the mpreg function 
# in the aim 1 base functions
# --------------------------------------

CI.swimex <- function(data) {
	# data: data frame that includes the variables referenced below
	#       (designed to be subset to a specific swim exposure category:
	#        non-swimmers, body immersion, head immersion, swallowed water)
	allci       <- mpreg(Y~1,dat=data)
	# function to get estimates of cumulative incidence and CIs from model objects
	getCI <- function(fit) {
		# fit : an object returned from coeftest w/ single parameters corresponding to log(cumulative incidence)
		est <- exp(fit[1,1])
		se <- fit[1,2]
		lb <- exp(log(est)-1.96*se)
		ub <- exp(log(est)+1.96*se)
		res <- c(est,lb,ub)
		return(res)
	}
	CIoverall <- getCI(allci)
	return(CIoverall)
}

# --------------------------------------
# convenience functions to obtain
# point estimates and SEs from
# model fits for plotting
# --------------------------------------

# swim exposure regs w/o interactions
# (need to sum 2 coefficients: anycontact + higher level of contact)
CIR.swimex <- function(fo,vcv) {
	# fo : fit object returned from lbreg or mpreg
	# vcv: variance covariance matrix from lbreg or mpreg
	nr <- nrow(vcv)
	lc <- c(0,1,1,rep(0,nr-3)) # linear combination of betas for the estimate
	est <- exp(t(lc)%*%fo[,1])
	se <- sqrt( t(lc)%*%vcv%*%lc )
	lb <- exp(log(est)-1.96*se)
	ub <- exp(log(est)+1.96*se)
	res <- c(est,lb,ub)
	return(res)
}

sens.swimex <- function(ad,Y="diarrheaci3",swim.exposure="bodycontact") {
	#-----------------------------------
	# sens.swimex
	# A wrapper function to run all of the analyses for a swim exposure level in the
	# 13 beaches data for the sensitivity analysis that varies length of follow-up from 1-10 days
	# ad    : analysis data frame with exposure, outcome, and covariates
	# Y     : string: variable name for the outcome
	# swim.exposure :string, variable name for the exposure level
	# file  : directory and file name of the output file to store all of the results (.RData)
	#-----------------------------------
	
	# create a general exposure variable
	# from the argument passed to the wrapper function
	# (used below)
	if(swim.exposure=="bodycontact") {
		cat("\n-----------------------------------")
		cat("\nRESULTS FOR BODY IMMERSION EXPOSURE")
		cat(paste("\nOutcome: ",Y,sep=""))
		cat("\n-----------------------------------\n")
	} else if(swim.exposure=="headunder") {
		cat("\n-----------------------------------")
		cat("\nRESULTS FOR HEAD IMMERSION EXPOSURE")
		cat(paste("\nOutcome: ",Y,sep=""))
		cat("\n-----------------------------------\n")
	} else if(swim.exposure=="swallwater") {
		cat("\n-------------------------------------")
		cat("\nRESULTS FOR SWALLOWED WATER EXPOSURE")
		cat(paste("\nOutcome: ",Y,sep=""))
		cat("\n-------------------------------------\n")
	} else {
		simpleError("Error: must specify one of -bodycontact- -headunder- or -swallwater- for swim.exposure")
	}
	ad$Y <- ad[,Y]
	ad$exposure <- ad[,swim.exposure]

		
	# --------------------------------------
	# Non-Swimmers
	# Summarize the number of people at risk 
	# and number of cumulative incident episodes
	#
	# Calculate unadjusted Cumulative Incidence Rates
	# and robust 95% CIs using an intercept model
	#
	# see the N.swimex function in the base functions
	# see the CI.swimex function in the base functions
	# --------------------------------------
	cat("\n-------------------------------------------")
	cat("\nCalculating Non-swimmer cases and exposed")
	cat("\nas well as cumulative incidence and robust")
	cat("\nSEs using intercept models (output suppressed)")
	cat("\nSee summary output later in the log file")
	cat("\n-------------------------------------------\n")
	Ns.noswim  <- N.swimex(subset(ad,anycontact=="No"))
	templog <- capture.output({ CIs.noswim <- CI.swimex(subset(ad,anycontact=="No"))})

	
	# --------------------------------------
	# Higher level exposure
	# (body immersion, head immersion, swallowed water)
	# Summarize the number of people at risk 
	# and number of cumulative incident episodes
	#
	# Calculate unadjusted Cumulative Incidence Rates
	# and robust 95% CIs using an intercept model
	#
	# see the N.swimex function in the base functions
	# see the CI.swimex function in the base functions
	# --------------------------------------
	cat("\n-------------------------------------------")
	cat("\nCalculating Exposed Swimmer cases and exposed")
	cat("\nas well as cumulative incidence and robust")
	cat("\nSEs using intercept models (output suppressed)")
	cat("\nSee summary output later in the log file")
	cat("\n-------------------------------------------\n")
	Ns  <- N.swimex(subset(ad,exposure=="Yes"))
	templog <- capture.output({CIs <- CI.swimex(subset(ad,exposure=="Yes")) })
		rm(templog)
		
	# --------------------------------------
	# regressions to estimate adjusted CIRs
	# compared to non-swimmers
	# --------------------------------------
	
	# adjustment covariates:
	# +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood
	
	# --------------------------------------
	# Pooled estimates
	# --------------------------------------
	### Pooled, overall estimates
	cat("\n-------------------------------------------")
	cat("\nOverall Pooled Results (all ages, all cond)")
	cat("\n-------------------------------------------\n")
	allfit <- mpreg(Y~anycontact+exposure +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ad,vcv=T)
	
	# --------------------------------------
	# Calculate Adjusted CIRs and 95% CIs
	# from the pooled and stratified models
	# group together age-stratified estimates
	# organized by overall/marine/freshwater
	# --------------------------------------
	CIRs <- CIR.swimex(allfit$fit,allfit$vcovCL)
		
	# --------------------------------------
	# Print the main results to the log file
	# --------------------------------------
	cat("\n-------------------------------------------")
	cat("\nSUMMARY OUTPUT")
	cat("\n-------------------------------------------\n")
	cat("\n-------------------------------------------")
	cat("\nNON-SWIMMERS -- Cases (n) and Exposed (N)")
	cat("\nOverall and stratified by marine/freshwater")
	cat("\n-------------------------------------------\n")
		print(Ns.noswim)
	cat("\n----------------------------------------------")
	cat("\nNON-SWIMMERS -- Cumulative Incidence (95% CIs)")
	cat("\nOverall and stratified by marine/freshwater")
	cat("\n----------------------------------------------\n")
		print(CIs.noswim)
	cat("\n-------------------------------------------")
	cat("\nSWIMMERS -- Cases (n) and Exposed (N)")
	cat("\nOverall and stratified by marine/freshwater")
	cat("\n-------------------------------------------\n")
		print(Ns)
	cat("\n-------------------------------------------")
	cat("\nSWIMMERS -- Cumulative Incidence (95% CIs)")
	cat("\nOverall and stratified by marine/freshwater")
	cat("\n-------------------------------------------\n")
		print(CIs)
	cat("\n----------------------------------------------")
	cat("\nAdjusted Cumulative Incidence Ratios (95% CIs)")
	cat("\nOverall and stratified by marine/freshwater")
	cat("\n----------------------------------------------\n")
		print(CIRs)
		
	# --------------------------------------
	# Save all of the results down to disk (.RData)
	# (exclude the data frames and glm objects
	#  to save space -- they are big objects)
	# --------------------------------------
	list(allfit=allfit,N.noswim=Ns.noswim,CI.noswim=CIs.noswim,N=Ns,CI=CIs,CIR=CIRs)
	
}


# --------------------------------------
# load the and pre-preprocess the 
# analysis dataset
# (refer to the base functions script
# for details on the pre-processing)
# --------------------------------------

ad <- preprocess.13beaches("~/dropbox/13beaches/data/final/13beaches-analysis.csv")

# --------------------------------------
# Run the separate analyses for each
# level of swimmer water exposure
# --------------------------------------

# --------------------------------------
# repeat the analysis for 1-10 days
# --------------------------------------
d1  <- sens.swimex(ad=ad,Y="diarrheaci1")
d2  <- sens.swimex(ad=ad,Y="diarrheaci2")
d3  <- sens.swimex(ad=ad,Y="diarrheaci3")
d4  <- sens.swimex(ad=ad,Y="diarrheaci4")
d5  <- sens.swimex(ad=ad,Y="diarrheaci5")
d6  <- sens.swimex(ad=ad,Y="diarrheaci6")
d7  <- sens.swimex(ad=ad,Y="diarrheaci7")
d8  <- sens.swimex(ad=ad,Y="diarrheaci8")
d9  <- sens.swimex(ad=ad,Y="diarrheaci9")
d10 <- sens.swimex(ad=ad,Y="diarrheaci10")


# --------------------------------------
# collate the objects across follow-up
# periods
# --------------------------------------
ds <- list(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10)
N.noswim <- sapply(ds,function(x) x$N.noswim)
N.swim <- sapply(ds,function(x) x$N)
CIRs <- sapply(ds,function(x) x$CIR)



# --------------------------------------
# Save the results
# --------------------------------------
rm(ad)
save.image("~/dropbox/13beaches/aim1-results/rawoutput/aim1-sensitivity-swim-exposure-length-of-follow-up.RData")



















