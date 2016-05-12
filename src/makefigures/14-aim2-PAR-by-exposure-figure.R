
# --------------------------------------
# 14-aim2-PAR-by-exposure-figure.R
# ben arnold (benarnold@berkeley.edu)
#
# description:
# plot the population attributable fraction (PAF)
# that corresponds with different %s of exposure
# in the population for swimming exposure and
# for swimming in contaminated water exposure
#
# --------------------------------------

# --------------------------------------
# input files:
#	aim2-PARswimex-diar.Rdata
#	aim2-PARswimex-gi.Rdata
#	aim2-PARswimex-dailygi.Rdata
#
# output files:
#	aim2-PAR-swimex.pdf
#	aim2-PAR-entero1600noswim.pdf
# --------------------------------------

# --------------------------------------
# preamble
# --------------------------------------

rm(list=ls())
library(RColorBrewer)
library(scales)

# source the base functions
source("~/13beaches/src/aim1/0-aim1-base-functions.R")
source("~/13beaches/src/aim2/0-aim2-base-functions.R")

# --------------------------------------
# PAF and PAR functions (NOT USED)
# calculate the PAF given an RR and
# a proportion exposed
# --------------------------------------

PAFcalc <- function(Pe,RR) {
	# calculate the population attributable fraction
	# under a counterfactual where the exposure is removed
	# see section 4.7 of Jewell 2004 Statistics for Epidemiology
	# RR : relative risk associated with exposure
	# Pe : proportion of the population exposed
	PAF <- Pe*(RR-1) / (1 + Pe*(RR-1))
	return(PAF)

}

PARcalc <- function(Pe,p0,p1) {
	# calculate the population attributable risk
	# under a counterfactual where the exposure is removed
	# see section 4.7 of Jewell 2004 Statistics for Epidemiology
	# p0 : risk in the unexposed
	# p1 : risk in the exposed
	# Pe : proportion of the population exposed
	PAR <- p1*Pe + p0*(1-Pe) - p0
	return(PAR)
}

PARcalc2 <- function(Pe,p0,RR) {
	# calculate the population attributable risk
	# for an exposure variable with three levels (0, 1, 2)
	# under a counterfactual where exposure level 2 is re-assigned 
	# to level 0 (leaving level 1 the same)
	# This is akin to assigning swimmers with Entero >35 CFU to
	# be non-swimmers
	# see section 4.7 of Jewell 2004 Statistics for Epidemiology
	# for the basic sketch of how to solve for this, though this
	# is a slightly more complicated case
	# note that PAR = P(obs) - P(counterfactual)
	# where P(obs) = Pe_0*p0 + Pe_1*p1 + Pe_2*p2
	# and P(count) = (Pe_0 + Pe_2)*p0 + Pe_1*p1
	# Thus, PAR = Pe_2(p2-p0) = Pe_2*p0(RR-1)
	# 
	# Pe   : proportion of the population exposed to level 2
	# p0   : risk in the unexposed
	# RR   : CIR for exposure level 2 versus unexposed (level 0)
	PAR <- Pe*p0*(RR-1)
	return(PAR)
}

# --------------------------------------
# ALL SWIM EXPOSURE, PAR FIGURE
# --------------------------------------

# --------------------------------------
# load the swim exposure regression output
# and PAR output (to get total Ns)
# --------------------------------------
load("~/dropbox/13beaches/aim2-results/rawoutput/aim2-PARswimex-diar.Rdata")
load("~/dropbox/13beaches/aim2-results/rawoutput/aim2-PARswimex-gi.Rdata")
load("~/dropbox/13beaches/aim2-results/rawoutput/aim2-PARswimex-dailygi.Rdata")


# --------------------------------------
# Grab the empirical exposure levels to
# get actual PARs to plot points
# --------------------------------------

# empirical levels of exposure in the population
pe <- N.body[,2]/N.total[,2]
pe

# --------------------------------------
# Grab PAR estimates and 95% CIs
# --------------------------------------

# Diarrhea PAR estimates and bootstrap 95% confidence intervals
pars.diar <- c(PARswimex.diar.0to4$bootest[1], PARswimex.diar.5to10$bootest[1], PARswimex.diar.11plus$bootest[1])*1000
pars.diar.lb <- c(PARswimex.diar.0to4$boot95lb[1], PARswimex.diar.5to10$boot95lb[1], PARswimex.diar.11plus$boot95lb[1])*1000
pars.diar.ub <- c(PARswimex.diar.0to4$boot95ub[1], PARswimex.diar.5to10$boot95ub[1], PARswimex.diar.11plus$boot95ub[1])*1000

# Gi illness PAR estimates and bootstrap 95% confidence intervals
pars.gi <- c(PARswimex.gi.0to4$bootest[1], PARswimex.gi.5to10$bootest[1], PARswimex.gi.11plus$bootest[1])*1000
pars.gi.lb <- c(PARswimex.gi.0to4$boot95lb[1], PARswimex.gi.5to10$boot95lb[1], PARswimex.gi.11plus$boot95lb[1])*1000
pars.gi.ub <- c(PARswimex.gi.0to4$boot95ub[1], PARswimex.gi.5to10$boot95ub[1], PARswimex.gi.11plus$boot95ub[1])*1000

# Days missed of daily activities PAR estimates and bootstrap 95% confidence intervals
pars.dailygi <- c(PARswimex.dailygi.0to4$bootest[1], PARswimex.dailygi.5to10$bootest[1], PARswimex.dailygi.11plus$bootest[1])*1000
pars.dailygi.lb <- c(PARswimex.dailygi.0to4$boot95lb[1], PARswimex.dailygi.5to10$boot95lb[1], PARswimex.dailygi.11plus$boot95lb[1])*1000
pars.dailygi.ub <- c(PARswimex.dailygi.0to4$boot95ub[1], PARswimex.dailygi.5to10$boot95ub[1], PARswimex.dailygi.11plus$boot95ub[1])*1000

# --------------------------------------
# Extrapolations from the estimates
# assuming a constant CIR irrespective
# of Prob(exposure)
# --------------------------------------
diar.par100    <- pars.diar/pe[2:4]
gi.par100      <- pars.gi/pe[2:4]
dailygi.par100 <- pars.dailygi/pe[2:4]


# print estimates for 100% and 50% exposure
diar.par100
diar.par100/2

gi.par100
gi.par100/2

dailygi.par100
dailygi.par100/2

# --------------------------------------
# Plot the swim exposure PAR estimates
# --------------------------------------

pdf("~/dropbox/13beaches/aim2-results/figs/aim2-PAR-swimex.pdf",width=4.5,height=12)
lo <- layout(mat=matrix(1:3,nrow=3,ncol=1))
op <- par(mar=c(4,4,4,2)+0.1,cex.axis=1.2)
cols <- brewer.pal(5,"Dark2")

	ytics <- seq(0,25,by=5)
plot(1,1,type="n",
	bty="n",
	xaxt="n",xlab="",xlim=c(0,1),
	yaxt="n",ylab="",ylim=range(ytics),
	las=1
	)
	axis(1,at=seq(0,1,by=0.1),labels=seq(0,100,by=10))
	axis(2,at=ytics,las=1)
	mtext("Proportion of Beachgoers Exposed to Swimming (%)",side=1,line=3)
	mtext("Diarrhea",side=3,line=1,at=-0.14,adj=0,cex=0.8)
	mtext("Episodes per 1,000 Attributable to Swimming",side=3,line=0,at=-0.14,adj=0,cex=0.8)
	mtext("A)",side=3,line=2,at=-0.14,adj=0,font=1,cex=1)

	segments(x0=0,x1=1,y0=0,lty=1,col="gray90")

	segments(x0=rep(0,3),x1=rep(1,3),y0=rep(0,3),y1=diar.par100,lty=2,lwd=1,col=alpha(cols,0.5))

	segments(x0=pe[2:4],y0=pars.diar.lb,y1=pars.diar.ub,lwd=2,col=cols[1:3])
	points(pe[2:4],pars.diar,pch=21,bg="white",lwd=2,cex=1.4,col=cols[1:3])

	text(pe[2],pars.diar[1],"Ages 0 - 4",pos=2,col=cols[1])
	text(pe[3],pars.diar[2],"Ages 5 - 10",pos=4,col=cols[2])
	text(pe[4],pars.diar[3],"Ages >10",pos=2,col=cols[3])
	
	ytics <- seq(0,35,by=5)
plot(1,1,type="n",
	bty="n",
	xaxt="n",xlab="",xlim=c(0,1),
	yaxt="n",ylab="",ylim=range(ytics),
	las=1
	)
	axis(1,at=seq(0,1,by=0.1),labels=seq(0,100,by=10))
	axis(2,at=ytics,las=1)
	mtext("Proportion of Beachgoers Exposed to Swimming (%)",side=1,line=3)
	mtext("GI illness",side=3,line=1,at=-0.14,adj=0,cex=0.8)
	mtext("Episodes per 1,000 Attributable to Swimming",side=3,line=0,at=-0.14,adj=0,cex=0.8)
	mtext("B)",side=3,line=2,at=-0.14,adj=0,font=1,cex=1)
	
	segments(x0=0,x1=1,y0=0,lty=1,col="gray90")

	segments(x0=rep(0,3),x1=rep(1,3),y0=rep(0,3),y1=gi.par100,lty=2,lwd=1,col=alpha(cols,0.5))

	segments(x0=pe[2:4],y0=pars.gi.lb,y1=pars.gi.ub,lwd=2,col=cols[1:3])
	points(pe[2:4],pars.gi,pch=21,bg="white",lwd=2,cex=1.4,col=cols[1:3])

	text(pe[2],pars.gi[1],"Ages 0 - 4",pos=2,col=cols[1])
	text(pe[3],pars.gi[2],"Ages 5 - 10",pos=4,col=cols[2])
	text(pe[4],pars.gi[3],"Ages >10",pos=2,col=cols[3])

	ytics <- seq(0,50,by=10)
plot(1,1,type="n",
	bty="n",
	xaxt="n",xlab="",xlim=c(0,1),
	yaxt="n",ylab="",ylim=range(ytics),
	las=1
	)
	axis(1,at=seq(0,1,by=0.1),labels=seq(0,100,by=10))
	axis(2,at=ytics,las=1)
	mtext("Proportion of Beachgoers Exposed to Swimming (%)",side=1,line=3)
	mtext("Missed Daily Activities (work, school, vacation)",side=3,line=1,at=-0.14,adj=0,cex=0.8)
	mtext("Days per 1,000 Attributable to Swimming",side=3,line=0,at=-0.14,adj=0,cex=0.8)
	mtext("C)",side=3,line=2,at=-0.14,adj=0,font=1,cex=1)
	
	segments(x0=0,x1=1,y0=0,lty=1,col="gray90")

	segments(x0=rep(0,3),x1=rep(1,3),y0=rep(0,3),y1=dailygi.par100,lty=2,lwd=1,col=alpha(cols,0.5))

	segments(x0=pe[2:4],y0=pars.dailygi.lb,y1=pars.dailygi.ub,lwd=2,col=cols[1:3])
	points(pe[2:4],pars.dailygi,pch=21,bg="white",lwd=2,cex=1.4,col=cols[1:3])

	text(pe[2],pars.dailygi[1],"Ages 0 - 4",pos=2,col=cols[1])
	text(pe[3],pars.dailygi[2],"Ages 5 - 10",pos=4,col=cols[2])
	text(pe[4],pars.dailygi[3],"Ages >10",pos=2,col=cols[3])

par(op)
dev.off()


# --------------------------------------
# SWIM EXPOSURE >35 CFU, PAR FIGURE
# --------------------------------------

# --------------------------------------
# load the swim exposure regression output
# and PAR output (to get total Ns)
# --------------------------------------
load("~/dropbox/13beaches/aim2-results/rawoutput/aim2-PARentero1600noswim-diar.Rdata")
load("~/dropbox/13beaches/aim2-results/rawoutput/aim2-PARentero1600noswim-gi.Rdata")
load("~/dropbox/13beaches/aim2-results/rawoutput/aim2-PARentero1600noswim-dailygi.Rdata")

# --------------------------------------
# Grab the empirical exposure levels to
# get actual PARs to plot points
# --------------------------------------

# empirical levels of exposure in the population
pe.noswim <- c(N.age0to4[3,2]/sum(N.age0to4[,2]),N.age5to10[3,2]/sum(N.age5to10[,2]),N.age11plus[3,2]/sum(N.age11plus[,2]))
pe.noswim



# --------------------------------------
# Grab PAR estimates and 95% CIs
# --------------------------------------

# Diarrhea PAR estimates and bootstrap 95% confidence intervals
pars.noswim.diar <- c(AR35noswim.diar.0to4$bootest[1], AR35noswim.diar.5to10$bootest[1], AR35noswim.diar.11plus$bootest[1])*1000
pars.noswim.diar.lb <- c(AR35noswim.diar.0to4$boot95lb[1], AR35noswim.diar.5to10$boot95lb[1], AR35noswim.diar.11plus$boot95lb[1])*1000
pars.noswim.diar.ub <- c(AR35noswim.diar.0to4$boot95ub[1], AR35noswim.diar.5to10$boot95ub[1], AR35noswim.diar.11plus$boot95ub[1])*1000

# Gi illness PAR estimates and bootstrap 95% confidence intervals
pars.noswim.gi <- c(AR35noswim.gi.0to4$bootest[1], AR35noswim.gi.5to10$bootest[1], AR35noswim.gi.11plus$bootest[1])*1000
pars.noswim.gi.lb <- c(AR35noswim.gi.0to4$boot95lb[1], AR35noswim.gi.5to10$boot95lb[1], AR35noswim.gi.11plus$boot95lb[1])*1000
pars.noswim.gi.ub <- c(AR35noswim.gi.0to4$boot95ub[1], AR35noswim.gi.5to10$boot95ub[1], AR35noswim.gi.11plus$boot95ub[1])*1000

# Days missed of daily activities PAR estimates and bootstrap 95% confidence intervals
pars.noswim.dailygi <- c(AR35noswim.dailygi.0to4$bootest[1], AR35noswim.dailygi.5to10$bootest[1], AR35noswim.dailygi.11plus$bootest[1])*1000
pars.noswim.dailygi.lb <- c(AR35noswim.dailygi.0to4$boot95lb[1], AR35noswim.dailygi.5to10$boot95lb[1], AR35noswim.dailygi.11plus$boot95lb[1])*1000
pars.noswim.dailygi.ub <- c(AR35noswim.dailygi.0to4$boot95ub[1], AR35noswim.dailygi.5to10$boot95ub[1], AR35noswim.dailygi.11plus$boot95ub[1])*1000

# --------------------------------------
# Extrapolations from the estimates
# assuming a constant CIR irrespective
# of Prob(exposure)
# --------------------------------------
diar.noswim.par100    <- pars.noswim.diar/pe.noswim
gi.noswim.par100      <- pars.noswim.gi/pe.noswim 
dailygi.noswim.par100 <- pars.noswim.dailygi/pe.noswim

# print estimates for 100% and 50% exposure
diar.noswim.par100
diar.noswim.par100/2

gi.noswim.par100
gi.noswim.par100/2

dailygi.noswim.par100
dailygi.noswim.par100/2

# --------------------------------------
# Plot the swim exposure >35 CFU
# PAR estimates
# A) Diarrhea
# B) GI illness
# C) Missed daily activities (omitted)
# --------------------------------------

pdf("~/dropbox/13beaches/aim2-results/figs/aim2-PAR-entero1600noswim.pdf",width=12,height=5)
lo <- layout(mat=matrix(1:2,nrow=1,ncol=2))
op <- par(mar=c(4,3,4,5)+0.1,cex.axis=1.2)
# cols <- brewer.pal(5,"Dark2")
cbPalette <- c("#000000","#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cols <- cbPalette[c(3,4,5)]
xtics <- seq(0,0.5,by=0.1)

# justification of titles relative to the Y-axis
xjust <- -0.04
cexscale <- 1

	ytics <- seq(0,30,by=5)
plot(1,1,type="n",
	bty="n",
	xaxt="n",xlab="",xlim=range(xtics),
	yaxt="n",ylab="",ylim=range(ytics),
	las=1
	)
	axis(1,at=xtics,labels=xtics*100)
	axis(2,at=ytics,las=1)
	mtext("Proportion of Beachgoers Exposed >35 Enterococcus CFU/100ml (%)",side=1,line=3,cex=cexscale)
	mtext("Episodes of Diarrhea per 1,000",side=3,line=1,at=xjust,adj=0,cex=cexscale)
	mtext("Attributable to Exposure Above Regulatory Guidelines",side=3,line=0,at=xjust,adj=0,cex=cexscale)
	mtext("a",side=3,line=1,at=xjust-0.04,adj=0,font=2,cex=1.4)

	# segments(x0=0,x1=max(xtics),y0=0,lty=1,col="gray90")

	segments(x0=rep(0,3),x1=rep(0.5,3),y0=rep(0,3),y1=diar.noswim.par100/2,lty=2,lwd=2,col=alpha(cols,1))

	segments(x0=pe.noswim,y0=pars.noswim.diar.lb,y1=pars.noswim.diar.ub,lwd=2,col=cols[1:3])
	points(pe.noswim,pars.noswim.diar,pch=21,bg="white",lwd=2,cex=1.4,col=cols[1:3])

	mtext("Ages 0 - 4",side=4,line=-0.25,at=diar.noswim.par100[1]/2,col=cols[1],cex=cexscale,las=1)
	mtext("Ages 5 - 10",side=4,line=-0.25,at=diar.noswim.par100[2]/2,col=cols[2],cex=cexscale,las=1)
	mtext("Ages >10",side=4,line=-0.25,at=diar.noswim.par100[3]/2,col=cols[3],cex=cexscale,las=1)
	
	# ytics <- seq(0,14,by=2)
plot(1,1,type="n",
	bty="n",
	xaxt="n",xlab="",xlim=range(xtics),
	yaxt="n",ylab="",ylim=range(ytics),
	las=1
	)
	axis(1,at=xtics,labels=xtics*100)
	axis(2,at=ytics,las=1)
	mtext("Proportion of Beachgoers Exposed >35 Enterococcus CFU/100ml (%)",side=1,line=3,cex=cexscale)
	mtext("Episodes of Gastrointestinal illness per 1,000",side=3,line=1,at=xjust,adj=0,cex=cexscale)
	mtext("Attributable to Exposure Above Regulatory Guidelines",side=3,line=0,at=xjust,adj=0,cex=cexscale)
	mtext("b",side=3,line=1,at=xjust-0.04,adj=0,font=2,cex=1.4)
	
	# segments(x0=0,x1=max(xtics),y0=0,lty=1,col="gray90")

	segments(x0=rep(0,3),x1=rep(0.5,3),y0=rep(0,3),y1=gi.noswim.par100/2,lty=2,lwd=2,col=alpha(cols,1))

	segments(x0=pe.noswim,y0=pars.noswim.gi.lb,y1=pars.noswim.gi.ub,lwd=2,col=cols[1:3])
	points(pe.noswim,pars.noswim.gi,pch=21,bg="white",lwd=2,cex=1.4,col=cols[1:3])

	mtext("Ages 0 - 4",side=4,line=-0.25,at=gi.noswim.par100[1]/2,col=cols[1],cex=cexscale,las=1)
	mtext("Ages 5 - 10",side=4,line=-0.25,at=gi.noswim.par100[2]/2,col=cols[2],cex=cexscale,las=1)
	mtext("Ages >10",side=4,line=-0.25,at=gi.noswim.par100[3]/2,col=cols[3],cex=cexscale,las=1)
	

	# ytics <- seq(0,14,by=2)
# plot(1,1,type="n",
	# bty="n",
	# xaxt="n",xlab="",xlim=range(xtics),
	# yaxt="n",ylab="",ylim=range(ytics),
	# las=1
	# )
	# axis(1,at=xtics,labels=xtics*100)
	# axis(2,at=ytics,las=1)
	# mtext("Proportion of Beachgoers Exposed >35 Enterococcus CFU/100ml (%)",side=1,line=3,cex=cexscale)
	# mtext("Days of Missed Daily Activities per 1,000",side=3,line=1,at=xjust,adj=0,cex=cexscale)
	# mtext("Attributable to Exposure Above Regulatory Guidelines",side=3,line=0,at=xjust,adj=0,cex=cexscale)
	# mtext("C)",side=3,line=2,at=xjust,adj=0,font=1,cex=1.2)
	
	# segments(x0=0,x1=max(xtics),y0=0,lty=1,col="gray90")

	# segments(x0=rep(0,3),x1=rep(0.5,3),y0=rep(0,3),y1=dailygi.noswim.par100/2,lty=2,lwd=1,col=alpha(cols,0.5))

	# segments(x0=pe.noswim,y0=pars.noswim.dailygi.lb,y1=pars.noswim.dailygi.ub,lwd=2,col=cols[1:3])
	# points(pe.noswim,pars.noswim.dailygi,pch=21,bg="white",lwd=2,cex=1.4,col=cols[1:3])

	# text(pe.noswim[1],pars.noswim.dailygi[1],"Ages 0 - 4",pos=2,col=cols[1],font=2,cex=cexscale)
	# text(pe.noswim[2],pars.noswim.dailygi[2]+2,"Ages 5 - 10",pos=4,col=cols[2],font=2,cex=cexscale)
	# text(pe.noswim[3],pars.noswim.dailygi[3]+0.5,"Ages >10",pos=2,col=cols[3],font=2,cex=cexscale)

par(op)
dev.off()



