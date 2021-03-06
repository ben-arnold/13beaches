# --------------------------------------
# 15-aim1-sens-swim-exposure-summary-figures.R
# ben arnold (benarnold@berkeley.edu)
#
# description:
# plot CIRs associated with swim exposure
# using summary regression output
#
# this script is identical to
# 2-aim1-swim-exposure-summary-figures.R
# but summarizes estimates for GI illness
# as the outcome (rather than diarrhea)
#
# --------------------------------------

# --------------------------------------
# input files:
#	aim1-sens-swim-exposure-regs-body-GI-illness.RData
#	aim1-sens-swim-exposure-regs-head-GI-illness.RData
#	aim1-sens-swim-exposure-regs-swall-GI-illness.RData
#
# output files:
#	aim1-sens-swim-exposure-GI-illness-byage.pdf
# --------------------------------------

# --------------------------------------
# preamble
# --------------------------------------

rm(list=ls())
library(RColorBrewer)

# --------------------------------------
# Load the body immersion output 
# and rename common objects
# --------------------------------------
load("~/dropbox/13beaches/aim1-results/rawoutput/aim1-sens-swim-exposure-regs-body-GI-illness.RData")
Ns.body   <- Ns
CIs.body  <- CIs
CIRs.body <- CIRs

# stratified beach output for forest plots
blist.body <- list(hu.fit,si.fit,wp.fit,we.fit,av.fit,bo.fit,dh.fit,ed.fit,fa.fit,gd.fit,ma.fit,mb.fit,su.fit)


# --------------------------------------
# Load the head immersion output 
# and rename common objects
# --------------------------------------
load("~/dropbox/13beaches/aim1-results/rawoutput/aim1-sens-swim-exposure-regs-head-GI-illness.RData")
Ns.head   <- Ns
CIs.head  <- CIs
CIRs.head <- CIRs

# stratified beach output for forest plots
blist.head <- list(hu.fit,si.fit,wp.fit,we.fit,av.fit,bo.fit,dh.fit,ed.fit,fa.fit,gd.fit,ma.fit,mb.fit,su.fit)

# --------------------------------------
# Load the swallowed water output 
# and rename common objects
# --------------------------------------
load("~/dropbox/13beaches/aim1-results/rawoutput/aim1-sens-swim-exposure-regs-swall-GI-illness.RData")
Ns.swall   <- Ns
CIs.swall  <- CIs
CIRs.swall <- CIRs

# stratified beach output for forest plots
blist.swall <- list(hu.fit,si.fit,wp.fit,we.fit,av.fit,bo.fit,dh.fit,ed.fit,fa.fit,gd.fit,ma.fit,mb.fit,su.fit)


# --------------------------------------
# Summary of Swim Exposure Risk
# Pooled Estimates, Stratified by Age
# --------------------------------------
swimexlabs <- c("Non Swimmers","Body Immersion","Head Immersion","Swallowed Water")

# Ns
N.all    <- rbind(Ns.noswim$Noverall[1,],Ns.body$Noverall[1,],Ns.head$Noverall[1,],Ns.swall$Noverall[1,])
N.0to4   <- rbind(Ns.noswim$Noverall[2,],Ns.body$Noverall[2,],Ns.head$Noverall[2,],Ns.swall$Noverall[2,])
N.5to10  <- rbind(Ns.noswim$Noverall[3,],Ns.body$Noverall[3,],Ns.head$Noverall[3,],Ns.swall$Noverall[3,])
N.11plus <- rbind(Ns.noswim$Noverall[4,],Ns.body$Noverall[4,],Ns.head$Noverall[4,],Ns.swall$Noverall[4,])
rownames(N.all) <- rownames(N.0to4) <- rownames(N.5to10) <- rownames(N.11plus) <- swimexlabs

# Cumulative Incidence Estimates
ci.all    <- rbind(CIs.noswim$CIoverall[1,],CIs.body$CIoverall[1,],CIs.head$CIoverall[1,],CIs.swall$CIoverall[1,])
ci.0to4   <- rbind(CIs.noswim$CIoverall[2,],CIs.body$CIoverall[2,],CIs.head$CIoverall[2,],CIs.swall$CIoverall[2,])
ci.5to10  <- rbind(CIs.noswim$CIoverall[3,],CIs.body$CIoverall[3,],CIs.head$CIoverall[3,],CIs.swall$CIoverall[3,])
ci.11plus <- rbind(CIs.noswim$CIoverall[4,],CIs.body$CIoverall[4,],CIs.head$CIoverall[4,],CIs.swall$CIoverall[4,])
rownames(ci.all) <- rownames(ci.0to4) <- rownames(ci.5to10) <- rownames(ci.11plus) <- swimexlabs

# Cumulative Incidence Ratios
cir.all    <- rbind(rep(NA,3),CIRs.body$CIRoverall[1,],CIRs.head$CIRoverall[1,],CIRs.swall$CIRoverall[1,])
cir.0to4   <- rbind(rep(NA,3),CIRs.body$CIRoverall[2,],CIRs.head$CIRoverall[2,],CIRs.swall$CIRoverall[2,])
cir.5to10  <- rbind(rep(NA,3),CIRs.body$CIRoverall[3,],CIRs.head$CIRoverall[3,],CIRs.swall$CIRoverall[3,])
cir.11plus <- rbind(rep(NA,3),CIRs.body$CIRoverall[4,],CIRs.head$CIRoverall[4,],CIRs.swall$CIRoverall[4,])
rownames(cir.all) <- rownames(cir.0to4) <- rownames(cir.5to10) <- rownames(cir.11plus) <- swimexlabs


pdf("~/dropbox/13beaches/aim1-results/figs/aim1-sens-swim-exposure-GI-illness-byage.pdf",width=17,height=5)
op <- par(mar=c(4,7,8,0)+0.1,xpd=FALSE)
# cols <- c(brewer.pal(9,"YlGn")[6],brewer.pal(9,"YlOrRd")[5:7])
cols <- c(brewer.pal(9,"BuGn")[7],brewer.pal(9,"YlGnBu")[6:8])
ytics <- seq(0,140,by=20)
# set up an empty plot
MidPts <- barplot(1:4,names.arg=NA,border=NA,col=NA,
	ylim=range(ytics),ylab="",yaxt="n",
	las=1,bty="n"
	)
	segments(x0=0,x1=max(MidPts+0.5),y0=ytics,lty=2,lwd=1,col="gray80")
	segments(x0=mean(MidPts[3:4]),y0=min(ytics)-10,y1=max(ytics)+30,lwd=2,col="gray80")
	axis(2,at=ytics,las=1)
	mtext("GI illness\nIncidence\nper 1000",side=2,line=3,las=1)
	
	# calculate X coordinates relative to the mid points for each group
	xspan <- 0.37
	xplus <- c(-xspan, -xspan/3, xspan/3, xspan)  # evenly distribute 4 datapoints around each midpoint
	
	
	x0to4   <- xplus+MidPts[1]
	x5to10  <- xplus+MidPts[2]
	x11plus <- xplus+MidPts[3]
	xall    <- xplus+MidPts[4]
	allxs   <- c(x0to4,x5to10,x11plus,xall)  # for table and quartile labels in header/footer
	labx <- MidPts[1]-xspan*1.5  # for left-hand labels in the header/footer

	# plot age 0 to 4 estimates
	segments(x0=x0to4,y0=ci.0to4[,"CIlb"]*1000,y1=ci.0to4[,"CIub"]*1000,lwd=2,col=cols)
	points(x0to4,ci.0to4[,"CI"]*1000,pch=16,cex=1.25,lwd=2,col=cols)
	text(x=x0to4,y=ci.0to4[,"CI"]*1000,sprintf("%1.0f",ci.0to4[,"CI"]*1000),pos=4,cex=0.7,col=cols,font=2)
	
	# plot age 5 to 10 estimates
	segments(x0=x5to10,y0=ci.5to10[,"CIlb"]*1000,y1=ci.5to10[,"CIub"]*1000,lwd=2,col=cols)
	points(x5to10,ci.5to10[,"CI"]*1000,pch=16,cex=1.25,lwd=2,col=cols)
	text(x=x5to10,y=ci.5to10[,"CI"]*1000,sprintf("%1.0f",ci.5to10[,"CI"]*1000),pos=4,cex=0.7,col=cols,font=2)
	
	# plot age > 10 estimates
	segments(x0=x11plus,y0=ci.11plus[,"CIlb"]*1000,y1=ci.11plus[,"CIub"]*1000,lwd=2,col=cols)
	points(x11plus,ci.11plus[,"CI"]*1000,pch=16,cex=1.25,lwd=2,col=cols)
	text(x=x11plus,y=ci.11plus[,"CI"]*1000,sprintf("%1.0f",ci.11plus[,"CI"]*1000),pos=4,cex=0.7,col=cols,font=2)
	
	# plot all age estimates
	segments(x0=xall,y0=ci.all[,"CIlb"]*1000,y1=ci.all[,"CIub"]*1000,lwd=2,col=cols)
	points(xall,ci.all[,"CI"]*1000,pch=16,cex=1.25,lwd=2,col=cols)
	text(x=xall,y=ci.all[,"CI"]*1000,sprintf("%1.0f",ci.all[,"CI"]*1000),pos=4,cex=0.7,col=cols,font=2)
	
	
	# print header labels
	mtext(c("Ages\n0 to 4","Ages\n5 to 10","Ages\n>10","All\nAges"),at=MidPts,side=3,line=5.5  )
	
	mtext("Swim\nExposure",side=3,line=3.5,at=labx,adj=1,col="gray30",cex=0.8)
	mtext(c("Non\nSwimmers","Body\nImmersion","Head\nImmersion","Swallowed\nWater"),side=3,line=3.5,at=allxs,col=cols,cex=0.8,font=1)
	
	
	# Print adjusted CIRs and 95% CIs (formatted)
	cirform <- function(cirs) {
		paste(sprintf("%1.2f",cirs),sep="")
	}
	circiform <- function(circi) {
		apply(circi,1,function(x) paste("(",sprintf("%1.2f",x["CIRlb"]),", ",sprintf("%1.2f",x["CIRub"]),")",sep="") )
	}
	mtext("Adjusted CIR",side=3,line=2,at=labx,adj=1,cex=0.8,col="gray30")
	mtext(c("ref",cirform(cir.0to4[2:4,1]),"ref",cirform(cir.5to10[2:4,1]),"ref",cirform(cir.11plus[2:4,1]),"ref",cirform(cir.all[2:4,1])),side=3,line=2,at=allxs,cex=0.75)
	
	mtext("(95% CI)",side=3,line=1,at=labx,adj=1,cex=0.8,col="gray30")
	mtext(c("",circiform(cir.0to4[2:4,2:3]),"",circiform(cir.5to10[2:4,2:3]),"",circiform(cir.11plus[2:4,2:3]),"",circiform(cir.all[2:4,2:3])),side=3,line=1,at=allxs,cex=0.7)
	

	# print footer labels
	mtext(c("Non\nSwimmers","Body\nImmersion","Head\nImmersion","Swallowed\nWater"),side=1,line=1,at=allxs,col=cols,cex=0.8,font=1)
	
	# print table with Ns
	mtext("Incident GI illness Cases",side=1,line=2,at=labx,adj=1,cex=0.8,col="gray30")
	ns <- c(N.0to4[,1],N.5to10[,1],N.11plus[,1],N.all[,1])
	mtext(  format(ns,big.mark=","),side=1,line=2,at=allxs+0.03,adj=1,cex=0.75    )
	
	mtext("Population At Risk",side=1,line=3,at=labx,adj=1,cex=0.8,col="gray30")
	Ns <- c(N.0to4[,2],N.5to10[,2],N.11plus[,2],N.all[,2])
	mtext(  format(Ns,big.mark=","),side=1,line=3,at=allxs+0.03,adj=1,cex=0.75    )
	
par(op)
dev.off()





