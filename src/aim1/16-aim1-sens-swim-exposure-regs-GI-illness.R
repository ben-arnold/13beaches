



# --------------------------------------
# 11-aim1-sens-swim-exposure-regs-GI-illness
# ben arnold (benarnold@berkeley.edu)
#
# description:
# estimate the adsociation between 
# water exposure and the risk of GI illness
# for the 13 beaches study
#
#  Note: this script is identical to
#  1-aim1-swim-exposure-regs.R, with a 
#  global replace for the outcome:
#  diarrheaci10 -> gici10
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

	N.all       <- calcNs(data$gici10)
	N.age0to4   <- calcNs(data$gici10[data$agestrat=="(0, 4]"])
	N.age5to10  <- calcNs(data$gici10[data$agestrat=="(4, 10]"])
	N.age11plus <- calcNs(data$gici10[data$agestrat==">10"])
	Noverall <- rbind(N.all,N.age0to4,N.age5to10,N.age11plus)
	
	N.marine       <- calcNs(data$gici10[data$marine=="Yes"])
	N.marine0to4   <- calcNs(data$gici10[data$marine=="Yes" & data$agestrat=="(0, 4]"])
	N.marine5to10  <- calcNs(data$gici10[data$marine=="Yes" & data$agestrat=="(4, 10]"])
	N.marine11plus <- calcNs(data$gici10[data$marine=="Yes" & data$agestrat==">10"])
	Nmarine <- rbind(N.marine,N.marine0to4,N.marine5to10,N.marine11plus)
	
	N.fresh       <- calcNs(data$gici10[data$marine=="No"])
	N.fresh0to4   <- calcNs(data$gici10[data$marine=="No" & data$agestrat=="(0, 4]"])
	N.fresh5to10  <- calcNs(data$gici10[data$marine=="No" & data$agestrat=="(4, 10]"])
	N.fresh11plus <- calcNs(data$gici10[data$marine=="No" & data$agestrat==">10"])
	Nfresh <- rbind(N.fresh,N.fresh0to4,N.fresh5to10,N.fresh11plus)
	rownames(Noverall) <- rownames(Nmarine) <- rownames(Nfresh) <- c("All Ages","Age 0 to 4","Age 5 to 10","Age >10")

	list(Noverall=Noverall,Nmarine=Nmarine,Nfresh=Nfresh)
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
	allci       <- mpreg(gici10~1,dat=data)
	age0to4ci   <- mpreg(gici10~1,dat=subset(data,agestrat=="(0, 4]"))
	age5to10ci  <- mpreg(gici10~1,dat=subset(data,agestrat=="(4, 10]"))
	age11plusci <- mpreg(gici10~1,dat=subset(data,agestrat==">10"))
	
	marineci       <- mpreg(gici10~1,dat=subset(data,marine=="Yes"))
	marine0to4ci   <- mpreg(gici10~1,dat=subset(data,marine=="Yes" & agestrat=="(0, 4]"))
	marine5to10ci  <- mpreg(gici10~1,dat=subset(data,marine=="Yes" & agestrat=="(4, 10]"))
	marine11plusci <- mpreg(gici10~1,dat=subset(data,marine=="Yes" & agestrat==">10"))
	
	freshci       <- mpreg(gici10~1,dat=subset(data,marine=="No"))
	fresh0to4ci   <- mpreg(gici10~1,dat=subset(data,marine=="No" & agestrat=="(0, 4]"))
	fresh5to10ci  <- mpreg(gici10~1,dat=subset(data,marine=="No" & agestrat=="(4, 10]"))
	fresh11plusci <- mpreg(gici10~1,dat=subset(data,marine=="No" & agestrat==">10"))
	
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
	
	ci.overall <- list(allci,age0to4ci,age5to10ci,age11plusci)
	CIoverall <- t(sapply(ci.overall,getCI))
	
	ci.marine <- list(marineci,marine0to4ci,marine5to10ci,marine11plusci)
	CImarine <- t(sapply(ci.marine,getCI))
	
	ci.fresh <- list(freshci,fresh0to4ci,fresh5to10ci,fresh11plusci)
	CIfresh <- t(sapply(ci.fresh,getCI))
	
	rownames(CIoverall) <- rownames(CImarine) <- rownames(CIfresh) <- c("All Ages","Age 0 to 4","Age 5 to 10","Age >10")
	colnames(CIoverall) <- colnames(CImarine) <- colnames(CIfresh) <- c("CI","CIlb","CIub")
	list(CIoverall=CIoverall,CImarine=CImarine,CIfresh=CIfresh)
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


est.swimex <- function(ad,swim.exposure="bodycontact",file="bodycontact-reg.RData") {
	#-----------------------------------
	# est.swimex
	# A wrapper function to run all of the analyses for a swim exposure level in the
	# 13 beaches data
	# ad       : analysis data frame with exposure, outcome, and covariates
	# exposure : string: exposure level variable (bodycontact, headunder, swallwater)
	# file     : directory and file name of the output file to store all of the results (.RData)
	#-----------------------------------
	
	
	# create a general exposure variable
	# from the argument passed to the wrapper function
	# (used below)
	if(swim.exposure=="bodycontact") {
		cat("\n-----------------------------------")
		cat("\nRESULTS FOR BODY IMMERSION EXPOSURE")
		cat("\n-----------------------------------\n")
	} else if(swim.exposure=="headunder") {
		cat("\n-------------------------------\nRESULTS FOR HEAD IMMERSION EXPOSURE\n-------------------------------\n")
	} else if(swim.exposure=="swallwater") {
		cat("\n-------------------------------\nRESULTS FOR SWALLOWED WATER EXPOSURE\n-------------------------------\n")
	} else {
		simpleError("Error: must specify one of -bodycontact- -headunder- or -swallwater- for swim.exposure")
	}
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
	# Note: marine and pointsource excluded from beach-stratified estimates 
	# (no variation since defined at the beach level)
	

	# tests of interaction by environmental conditions for Avalon and Doheny
	# (no evidence for interaction)
	cat("\n-------------------------------------------")
	cat("\nAvalon beach")
	cat("\nTest of interaction by environmental conditions")
	cat("\nAbove versus below median groundwater flow")
	cat("\nSee Yau et al. 2014 Wat Res for details")
	cat("\n-------------------------------------------\n")
	av.i.noint <-glm(gici10~anycontact+exposure +groundwater +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,family=poisson(link="log"),data=ad[ad$beach=="Avalon",])
	av.i <-glm(gici10~(anycontact+exposure)*groundwater +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,family=poisson(link="log"),data=ad[ad$beach=="Avalon",])
		print(lrtest(av.i.noint,av.i))
	
	cat("\n-------------------------------------------")
	cat("\nDoheny beach")
	cat("\nTest of interaction by environmental conditions")
	cat("\nSand berm open vs. closed")
	cat("\nSee Colford et al. 2012 Wat Res for details")
	cat("\n-------------------------------------------\n")
	dh.i.noint <-glm(gici10~anycontact+exposure +berm +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,family=poisson(link="log"),data=ad[ad$beach=="Doheny",])
	dh.i <-glm(gici10~(anycontact+exposure)*berm +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,family=poisson(link="log"),data=ad[ad$beach=="Doheny",])
		print(lrtest(dh.i.noint,dh.i))

	# --------------------------------------
	# Beach-specific stratified estimates
	# --------------------------------------
	# Freshwater beaches
	# --------------------------------------
	
	
	cat("\n-------------------------------------------")
	cat("\nBEACH STRATIFIED ESTIMATES")
	cat("\n-------------------------------------------\n")
	
	# Huntington
	cat("\n-------------------------------------------")
	cat("\nHuntington")
	cat("\n-------------------------------------------\n")
	hu.fit <- mpreg(gici10~anycontact+exposure +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ad[ad$beach=="Huntington",],vcv=T)
	
	# Silver
	cat("\n-------------------------------------------")
	cat("\nSilver")
	cat("\n-------------------------------------------\n")
	si.fit <- mpreg(gici10~anycontact+exposure +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ad[ad$beach=="Silver",],vcv=T)
	
	# Washington Park
	cat("\n-------------------------------------------")
	cat("\nWashington Park")
	cat("\n-------------------------------------------\n")
	wp.fit <- mpreg(gici10~anycontact+exposure +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ad[ad$beach=="Washington Park",],vcv=T)
	
	# West
	cat("\n-------------------------------------------")
	cat("\nWest")
	cat("\n-------------------------------------------\n")
	we.fit <- mpreg(gici10~anycontact+exposure +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ad[ad$beach=="West",],vcv=T)
	
	# --------------------------------------
	# Marine beaches
	# --------------------------------------
	
	# Avalon
	cat("\n-------------------------------------------")
	cat("\nAvalon")
	cat("\n-------------------------------------------\n")
	av.fit <- mpreg(gici10~anycontact+exposure +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ad[ad$beach=="Avalon",],vcv=T)
	
	# Boqueron
	cat("\n-------------------------------------------")
	cat("\nBoqueron")
	cat("\n-------------------------------------------\n")
	bo.fit <- mpreg(gici10~anycontact+exposure +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ad[ad$beach=="Boqueron",],vcv=T)
	
	# Doheny
	cat("\n-------------------------------------------")
	cat("\nDoheny")
	cat("\n-------------------------------------------\n")
	dh.fit <- mpreg(gici10~anycontact+exposure +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ad[ad$beach=="Doheny",],vcv=T)
	
	# Edgewater
	cat("\n-------------------------------------------")
	cat("\nEdgewater")
	cat("\n-------------------------------------------\n")
	ed.fit <- mpreg(gici10~anycontact+exposure +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ad[ad$beach=="Edgewater",],vcv=T)
	
	# Fairhope
	cat("\n-------------------------------------------")
	cat("\nFairhope")
	cat("\n-------------------------------------------\n")
	fa.fit <- mpreg(gici10~anycontact+exposure +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ad[ad$beach=="Fairhope",],vcv=T)
	
	# Goddard
	cat("\n-------------------------------------------")
	cat("\nGoddard")
	cat("\n-------------------------------------------\n")
	gd.fit <- mpreg(gici10~anycontact+exposure +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ad[ad$beach=="Goddard",],vcv=T)
	
	# Malibu
	cat("\n-------------------------------------------")
	cat("\nMalibu")
	cat("\n-------------------------------------------\n")
	ma.fit <- mpreg(gici10~anycontact+exposure +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ad[ad$beach=="Malibu",],vcv=T)
	
	# Mission Bay
	cat("\n-------------------------------------------")
	cat("\nMission Bay")
	cat("\n-------------------------------------------\n")
	mb.fit <- mpreg(gici10~anycontact+exposure +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ad[ad$beach=="Mission Bay",],vcv=T)
	
	# Surfside
	cat("\n-------------------------------------------")
	cat("\nSurfside")
	cat("\n-------------------------------------------\n")
	su.fit <- mpreg(gici10~anycontact+exposure +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ad[ad$beach=="Surfside",],vcv=T)

	# --------------------------------------
	# Pooled estimates
	# --------------------------------------
	### Pooled, overall estimates
	cat("\n-------------------------------------------")
	cat("\nOverall Pooled Results (all ages, all cond)")
	cat("\n-------------------------------------------\n")
	allfit <- glm(gici10~anycontact+exposure +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,family=poisson(link="log"),data=ad)
		all.fit.VC <- cl(ad,fm=allfit,cluster=ad$hhid)
		all.fit <- coeftest(allfit, all.fit.VC)
		print(all.fit)
	
	### Pooled estimate (point vs. non-point source conditions)
	cat("\n-------------------------------------------")
	cat("\nTest of interaction (effect modification)")
	cat("\nfor point vs. non-point beach conditions")
	cat("\n-------------------------------------------\n")
	psref <- glm(gici10~anycontact+exposure+pointsource +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,family=poisson(link="log"),data=ad)
	psfit <- glm(gici10~(anycontact+exposure)*pointsource +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,family=poisson(link="log"),data=ad)
		ps.fit.VC <- cl(ad,fm=psfit,cluster=ad$hhid)
		ps.fit <- coeftest(psfit, ps.fit.VC)
		print(ps.fit)
		print(lrtest(psref,psfit))
	
	
	### Pooled estimate (fresh vs. marine)
	cat("\n-------------------------------------------")
	cat("\nTest of interaction (effect modification)")
	cat("\nfor fresh vs. marine beach conditions")
	cat("\n-------------------------------------------\n")
	fmref <- glm(gici10~anycontact+exposure +marine +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,family=poisson(link="log"),data=ad)
	fmfit <- glm(gici10~(anycontact+exposure)*marine +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,family=poisson(link="log"),data=ad)
		fm.fit.VC <- cl(ad,fm=fmfit,cluster=ad$hhid)
		fm.fit <- coeftest(fmfit, fm.fit.VC)
		print(fm.fit)
		print(lrtest(fmref,fmfit))
	
	# --------------------------------------
	# Age-stratified estimates and LR tests of
	# interaction (uses a different age adjustement 
	# covariate compared to the other pooled models:
	# "agestrat" vs. "agecat"
	
	### reduced models for LR tests of exposure x age interactions
	agestratref <- glm(gici10~anycontact+exposure +agestrat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,family=poisson(link="log"),data=ad)
	
	### Pooled estimate (Age 0-4, 5-10, >10)
	cat("\n-------------------------------------------")
	cat("\nTest of interaction (effect modification)")
	cat("\nby age category (0-4, 5-10, >10)")
	cat("\n-------------------------------------------\n")
	agestratfit <- glm(gici10~(anycontact+exposure)*agestrat +pointsource+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,family=poisson(link="log"),data=ad)
		agestrat.fit.VC <- cl(ad,fm=agestratfit,cluster=ad$hhid)
		agestrat.fit <- coeftest(agestratfit, agestrat.fit.VC) 
		print(agestrat.fit)
		print(lrtest(agestratref,agestratfit))
	
	# --------------------------------------
	# Stratified models based on
	# above tests of interaction
	# --------------------------------------
	
	cat("\n-------------------------------------------")
	cat("\nSTRATIFIED, POOLED MODELS")
	cat("\n-------------------------------------------\n")
	
	# all conditions - All Ages
	# see all.fit (above)
	cat("\n-------------------------------------------")
	cat("\nAll Ages")
	cat("\n-------------------------------------------\n")
	print(all.fit)

	# all conditions - Age 0 to 4
	cat("\n-------------------------------------------")
	cat("\nAges 0 to 4")
	cat("\n-------------------------------------------\n")
	age0to4.fit <- mpreg(gici10~anycontact+exposure +female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ad[ad$agestrat=="(0, 4]",],vcv=T)
	
	# all conditions - Age 5 to 10
	cat("\n-------------------------------------------")
	cat("\nAges 5 to 10")
	cat("\n-------------------------------------------\n")
	age5to10.fit <- mpreg(gici10~anycontact+exposure +female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ad[ad$agestrat=="(4, 10]",],vcv=T)
	
	# all conditions - Age >10
	cat("\n-------------------------------------------")
	cat("\nAges >10")
	cat("\n-------------------------------------------\n")
	age11plus.fit <- mpreg(gici10~anycontact+exposure +female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ad[ad$agestrat==">10",],vcv=T)
	
	# marine beaches - All Ages
	cat("\n-------------------------------------------")
	cat("\nAll Ages -- Marine Beaches")
	cat("\n-------------------------------------------\n")
	marine.fit <- mpreg(gici10~anycontact+exposure +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ad[ad$marine=="Yes",],vcv=T)
	
	# marine beaches - Age 0 to 4
	cat("\n-------------------------------------------")
	cat("\nAges 0 to 4 -- Marine Beaches")
	cat("\n-------------------------------------------\n")
	marine0to4.fit <- mpreg(gici10~anycontact+exposure +female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ad[ad$marine=="Yes" & ad$agestrat=="(0, 4]",],vcv=T)
	
	# marine beaches - Age 5 to 10
	cat("\n-------------------------------------------")
	cat("\nAges 5 to 10 -- Marine Beaches")
	cat("\n-------------------------------------------\n")
	marine5to10.fit <- mpreg(gici10~anycontact+exposure +female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ad[ad$marine=="Yes" & ad$agestrat=="(4, 10]",],vcv=T)
	
	# marine beaches - Age >10
	cat("\n-------------------------------------------")
	cat("\nAges >10 -- Marine Beaches")
	cat("\n-------------------------------------------\n")
	marine11plus.fit <- mpreg(gici10~anycontact+exposure +female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ad[ad$marine=="Yes" & ad$agestrat==">10",],vcv=T)
	
	# freshwater beaches - All Ages
	cat("\n-------------------------------------------")
	cat("\nAll Ages -- Freshwater Beaches")
	cat("\n-------------------------------------------\n")
	fresh.fit <- mpreg(gici10~anycontact+exposure +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ad[ad$marine=="No",],vcv=T)
	
	# freshwater beaches - Age 0 to 4
	cat("\n-------------------------------------------")
	cat("\nAges 0 to 4 -- Freshwater Beaches")
	cat("\n-------------------------------------------\n")
	fresh0to4.fit <- mpreg(gici10~anycontact+exposure +female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ad[ad$marine=="No" & ad$agestrat=="(0, 4]",],vcv=T)
	
	# freshwater beaches - Age 5 to 10
	cat("\n-------------------------------------------")
	cat("\nAges 5 to 10 -- Freshwater Beaches")
	cat("\n-------------------------------------------\n")
	fresh5to10.fit <- mpreg(gici10~anycontact+exposure +female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ad[ad$marine=="No" & ad$agestrat=="(4, 10]",],vcv=T)
	
	# freshwater beaches - Age >10
	cat("\n-------------------------------------------")
	cat("\nAges >10 -- Freshwater Beaches")
	cat("\n-------------------------------------------\n")
	fresh11plus.fit <- mpreg(gici10~anycontact+exposure +female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ad[ad$marine=="No" & ad$agestrat==">10",],vcv=T)
	

	# --------------------------------------
	# Calculate Adjusted CIRs and 95% CIs
	# from the pooled and stratified models
	# group together age-stratified estimates
	# organized by overall/marine/freshwater
	# --------------------------------------
	
	CIRoverall <- rbind(
		CIR.swimex(all.fit,all.fit.VC),
		CIR.swimex(age0to4.fit$fit,age0to4.fit$vcovCL),
		CIR.swimex(age5to10.fit$fit,age5to10.fit$vcovCL),
		CIR.swimex(age11plus.fit$fit,age11plus.fit$vcovCL)
	)
	CIRfresh <- rbind(
		CIR.swimex(fresh.fit$fit,fresh.fit$vcovCL),
		CIR.swimex(fresh0to4.fit$fit,fresh0to4.fit$vcovCL),
		CIR.swimex(fresh5to10.fit$fit,fresh5to10.fit$vcovCL),
		CIR.swimex(fresh11plus.fit$fit,fresh11plus.fit$vcovCL)
	)
	CIRmarine <- rbind(
		CIR.swimex(marine.fit$fit,marine.fit$vcovCL),
		CIR.swimex(marine0to4.fit$fit,marine0to4.fit$vcovCL),
		CIR.swimex(marine5to10.fit$fit,marine5to10.fit$vcovCL),
		CIR.swimex(marine11plus.fit$fit,marine11plus.fit$vcovCL)
	)
	rownames(CIRoverall) <- rownames(CIRfresh) <- rownames(CIRmarine) <- c("All Ages","Age 0 to 4","Age 5 to 10","Age >10")
	colnames(CIRoverall) <- colnames(CIRfresh) <- colnames(CIRmarine) <- c("CIR","CIRlb","CIRub")
	CIRs <- list(CIRoverall=CIRoverall,CIRmarine=CIRmarine,CIRfresh=CIRfresh)
	
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
	save(
		hu.fit,si.fit,wp.fit,we.fit,av.fit,bo.fit,dh.fit,ed.fit,fa.fit,gd.fit,ma.fit,mb.fit,su.fit,
		all.fit,age0to4.fit,age5to10.fit,age11plus.fit,
		marine.fit,marine0to4.fit,marine5to10.fit,marine11plus.fit,
		fresh.fit,fresh0to4.fit,fresh5to10.fit,fresh11plus.fit,
		Ns.noswim,CIs.noswim,
		Ns,CIs,CIRs,
	file=file)
	
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
# BODY IMMERSION
# --------------------------------------
est.swimex(ad=ad,swim.exposure="bodycontact",file="~/dropbox/13beaches/aim1-results/rawoutput/aim1-sens-swim-exposure-regs-body-GI-illness.RData")

# --------------------------------------
# HEAD IMMERSION
# --------------------------------------
est.swimex(ad=ad,swim.exposure="headunder",file="~/dropbox/13beaches/aim1-results/rawoutput/aim1-sens-swim-exposure-regs-head-GI-illness.RData")

# --------------------------------------
# SWALLOWED WATER
# --------------------------------------
est.swimex(ad=ad,swim.exposure="swallwater",file="~/dropbox/13beaches/aim1-results/rawoutput/aim1-sens-swim-exposure-regs-swall-GI-illness.RData")




