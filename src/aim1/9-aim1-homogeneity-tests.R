
# ------------------------------------
# 9-aim1-homogeneity-tests.R
#
# Calculate Cochran's Q statistic and
# the related Higgens' I squared for
# swim exposure analyses and for
# the Enterococcus analyses
#
# ------------------------------------

# ------------------------------------
# preamble
# ------------------------------------
rm(list=ls())


# ------------------------------------
# Function to calculate Q, P(Q), and I2
# ------------------------------------

I2calc <- function(studyCIRs,overallCIR) {
  # calculate Cochran's Q, P(Q), and Higgens I2
	# given study-specific CIRs and an overall CIR
	#
	# See Hedges & Vivea 1998, equation (7)
	# and Higgins et al. BMJ 2003;327:557–60
  #
  # Note: this function was validated against Figure 1 of Higgins 2003
	#
	# studyCIRs : a 3 column matrix with CIR, CIRlb, CIRub as columns
	# overallCIR : the overall CIR across studies (scalar)
	
	# calculate the weights for each study
	# equal to the inverse variance of each estimate
	wi = 1/ ( ( (log(studyCIRs[,3]) - log(studyCIRs[,1]))/1.96 ) ^2 )
	
	# calculate squared deviations for each estimate 
	# from the overall estimate on the log scale
	bdiff2 <- ( log(studyCIRs[,1]) - log(overallCIR)  )^2
	
	# weight the squared differences by the inverse variance
	bwdiff2 <- bdiff2*wi
	
	# calculate the weighted mean squared deviations
	Q <- sum(bwdiff2)
	
	# calculate a P-value using the Chi-squared distribution with k-1 df
	k <- nrow(studyCIRs)
	P <- 1-pchisq(Q,df=k-1)
	
	# calculate I2
	I2 <- 100*( (Q - (k-1))/Q)
	
	cat("\n Cochran's Q = ",Q,"\n Pr(Q,df=",k-1,")   = ",P,"\n I squared   = ",I2,"\n\n",sep="")
	
	return(rbind(Q,P,I2))

}


# ------------------------------------
# Swim exposure heterogeneity tests
# ------------------------------------

# ------------------------------------
# Calculate CIRs from beach regression output
# ------------------------------------
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

# beach name labels that correspond to the beach extraction order used below
fCIRlab <- c("Huntington","Silver","Washington Park","West","Avalon","Boqueron","Edgewater","Doheny","Fairhope","Goddard","Malibu","Mission Bay","Surfside")


# body immersion results
load("~/dropbox/13beaches/aim1-results/rawoutput/aim1-swim-exposure-regs-body.Rdata")
bCIRs <- rbind(
	CIR.swimex(hu.fit$fit,hu.fit$vcovCL),
	CIR.swimex(si.fit$fit,si.fit$vcovCL),
	CIR.swimex(wp.fit$fit,wp.fit$vcovCL),
	CIR.swimex(we.fit$fit,we.fit$vcovCL),
	CIR.swimex(av.fit$fit,av.fit$vcovCL),
	CIR.swimex(bo.fit$fit,bo.fit$vcovCL),
	CIR.swimex(ed.fit$fit,ed.fit$vcovCL),
	CIR.swimex(dh.fit$fit,dh.fit$vcovCL),
	CIR.swimex(fa.fit$fit,fa.fit$vcovCL),
	CIR.swimex(gd.fit$fit,gd.fit$vcovCL),
	CIR.swimex(ma.fit$fit,ma.fit$vcovCL),
	CIR.swimex(mb.fit$fit,mb.fit$vcovCL),
	CIR.swimex(su.fit$fit,su.fit$vcovCL)
	)
cir.body <- CIRs$CIRoverall[1,1]
	
# head immersion results
load("~/dropbox/13beaches/aim1-results/rawoutput/aim1-swim-exposure-regs-head.Rdata")
hCIRs <- rbind(
	CIR.swimex(hu.fit$fit,hu.fit$vcovCL),
	CIR.swimex(si.fit$fit,si.fit$vcovCL),
	CIR.swimex(wp.fit$fit,wp.fit$vcovCL),
	CIR.swimex(we.fit$fit,we.fit$vcovCL),
	CIR.swimex(av.fit$fit,av.fit$vcovCL),
	CIR.swimex(bo.fit$fit,bo.fit$vcovCL),
	CIR.swimex(ed.fit$fit,ed.fit$vcovCL),
	CIR.swimex(dh.fit$fit,dh.fit$vcovCL),
	CIR.swimex(fa.fit$fit,fa.fit$vcovCL),
	CIR.swimex(gd.fit$fit,gd.fit$vcovCL),
	CIR.swimex(ma.fit$fit,ma.fit$vcovCL),
	CIR.swimex(mb.fit$fit,mb.fit$vcovCL),
	CIR.swimex(su.fit$fit,su.fit$vcovCL)
	)
cir.head <- CIRs$CIRoverall[1,1]
		
# swallowed water results
load("~/dropbox/13beaches/aim1-results/rawoutput/aim1-swim-exposure-regs-swall.Rdata")
sCIRs <- rbind(
	CIR.swimex(hu.fit$fit,hu.fit$vcovCL),
	CIR.swimex(si.fit$fit,si.fit$vcovCL),
	CIR.swimex(wp.fit$fit,wp.fit$vcovCL),
	CIR.swimex(we.fit$fit,we.fit$vcovCL),
	CIR.swimex(av.fit$fit,av.fit$vcovCL),
	CIR.swimex(bo.fit$fit,bo.fit$vcovCL),
	CIR.swimex(ed.fit$fit,ed.fit$vcovCL),
	CIR.swimex(dh.fit$fit,dh.fit$vcovCL),
	CIR.swimex(fa.fit$fit,fa.fit$vcovCL),
	CIR.swimex(gd.fit$fit,gd.fit$vcovCL),
	CIR.swimex(ma.fit$fit,ma.fit$vcovCL),
	CIR.swimex(mb.fit$fit,mb.fit$vcovCL),
	CIR.swimex(su.fit$fit,su.fit$vcovCL)
	)
cir.swall <- CIRs$CIRoverall[1,1]

rownames(bCIRs) <- rownames(hCIRs) <- rownames(sCIRs) <- fCIRlab

# ------------------------------------
# run Cochran heterogeneity tests and
# calculate I2 for:
# body immersion
# head immersion
# swallowed water
# ------------------------------------
bI <- I2calc(bCIRs,cir.body)
hI <- I2calc(hCIRs,cir.head)
sI <- I2calc(sCIRs,cir.swall)

# ------------------------------------
# Enterococcus regulatory heterogeneity
# tests
# ------------------------------------

# ------------------------------------
# Enterococcus EPA 1600
# 35cfu results
# ------------------------------------
load("~/dropbox/13beaches/aim1-results/rawoutput/aim1-entero1600-35cfu-regs-body.Rdata")

# grab beach-specific CIRs and 95% CIs
CIRreg <- function(x) {
    # x: a list object returned from lbreg or mpreg, that includes $fit and $vcv
	fo <- x$fit
	vcv <- x$vcovCL
	nr <- nrow(vcv)
	lc <- c(0,1,rep(0,nr-2)) # linear combination of betas for the estimate
	est <- exp(t(lc)%*%fo[,1])
	se <- sqrt( t(lc)%*%vcv%*%lc )
	lb <- exp(log(est)-1.96*se)
	ub <- exp(log(est)+1.96*se)
	res <- c(est,lb,ub)
	return(res)
}
fitlist <- list(hufit,sifit,wpfit,wefit,avfit,bofit,dhfit,edfit,fafit,gdfit,mafit,mbfit,sufit)
entero1600.cirs <- t(sapply(fitlist,CIRreg))

# Malibu and Goddard do not have actual estimates -- remove them
entero1600.cirs <- entero1600.cirs[-c(10:11),]

# run tests for Entero 1600
entero1600I <- I2calc(entero1600.cirs,cir.all[3,1])



# ------------------------------------
# Enterococcus EPA 1611 qPCR
# 470cce results
# ------------------------------------
load("~/dropbox/13beaches/aim1-results/rawoutput/aim1-enteroQPCR-470cce-regs-body.Rdata")
fitlist <- list(hufit,sifit,wpfit,wefit,avfit,bofit,dhfit,edfit,fafit,gdfit,mafit,mbfit,sufit)
enteroQPCR.cirs <- t(sapply(fitlist,CIRreg))

# Malibu and Mission Bay do not have actual estimates -- remove them
enteroQPCR.cirs <- enteroQPCR.cirs[-c(11:12),]

# run tests for Entero 1611
enteroQPCRI <- I2calc(enteroQPCR.cirs,cir.all[3,1])


# ------------------------------------
# Save the results to a file
# ------------------------------------
save(bI,sI,hI,entero1600I,enteroQPCRI,file="~/dropbox/13beaches/aim1-results/rawoutput/aim1-homogeneity-tests.RData")








