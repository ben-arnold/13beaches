# --------------------------------------
# 12-aim1-sens-entero1600-35cfu-length-of-followup-figure.R
# ben arnold (benarnold@berkeley.edu)
#
# description:
# plot CIRs associated with swim exposure
# using summary regression output
# with varying lengths of follow-up
# from 1-10 days
#
# --------------------------------------

# --------------------------------------
# input files:
#	aim1-sens-entero1600-35cfu-length-of-follow-up.RData
#
# output files:
#	aim1-sens-entero1600-35cfu-length-of-follow-up.pdf
# --------------------------------------

# --------------------------------------
# preamble
# --------------------------------------

rm(list=ls())
library(RColorBrewer)

# --------------------------------------
# Load the length of follow-up sensitivity
# output
# --------------------------------------
load("~/dropbox/13beaches/aim1-results/rawoutput/aim1-sens-entero1600-35cfu-length-of-follow-up.RData")


# --------------------------------------
# Plotting function for the results
# --------------------------------------

fusensplot <- function(cir.all,N.all,main,labtext=FALSE) {
	# plotting function to be repeated for stratified and overall conditions
	# cir.all : CIR estimates and 95% CIs (see the actual "cir.all" object)
	# N.all   : number of cumulative incident cases and individuals at risk
	# main    : title for the plot
	# labtext   : logical. include a label for the primary outcome?
	

	op <- par(mar=c(10,22,2,1)+0.1)
	col  <- brewer.pal(9,"YlGnBu")[7]
	ncol <- brewer.pal(9,"YlGnBu")[6]
	scol <- brewer.pal(9,"YlGnBu")[8]
	ytics <- c(0.7,1,1.2,1.5,2)
	plot(1:10,cir.all[1,],type="n",bty="n",
		xaxt="n",xlab="",xlim=c(1,10),
		yaxt="n",ylab="",log="y",ylog=T,ylim=range(ytics),
		las=1,main=""
		)
		axis(1,at=1:10,cex.axis=1.25)
		axis(2,at=ytics,cex.axis=1.25,las=1)
		segments(x0=0,x1=10,y0=ytics[-c(2)],col="gray80",lty=2)
		segments(x0=0,x1=10,y0=1,col="gray40",lty=2,lwd=2)
		mtext(main,side=3,line=0,cex=1.5)
		mtext("Adjusted CIR\nDiarrhea",side=2,line=3,at=1.5,adj=1,cex=1,las=1)
	
		# plot the estimates
		segments(x0=1:10,y0=cir.all[2,],y1=cir.all[3,],col=col,lwd=2)
		points(1:10,cir.all[1,],pch=c(rep(16,9),21),bg="white",col=col,lwd=2,cex=2)
		
		if (labtext==TRUE) text(10,1.75,"The primary analysis used 10 days of follow-up",cex=0.9,adj=1)
		
		# make a table below
		mtext("Days of Follow-up",side=1,line=1,at=0,adj=1)
		
		mtext("Enterococcus <=35 CFU/100ml",side=1,line=2.5,at=0,adj=1,cex=0.8,col=ncol,font=2)
		mtext("At Risk",side=1,line=3.5,at=0,adj=1,cex=0.8,col=ncol)
		mtext(format(N.all[3,1],big.mark=","),side=1,line=3.5,at=1+0.15,cex=0.8,adj=1,col="gray40")
		mtext("Cumulative Incident Cases",side=1,line=4.5,at=0,adj=1,cex=0.8,col=ncol)
		mtext(format(N.all[1,],big.mark=","),side=1,line=4.5,at=1:10+0.15,cex=0.8,adj=1,col="gray40")
		
		
		mtext("Enterococcus >35 CFU/100ml",side=1,line=6,at=0,adj=1,cex=0.8,col=scol,font=2)
		mtext("At Risk",side=1,line=7,at=0,adj=1,cex=0.8,col=scol)
		mtext(format(N.all[4,1],big.mark=","),side=1,line=7,at=1+0.15,cex=0.8,adj=1,col="gray40")
		mtext("Cumulative Incident Cases",side=1,line=8,at=0,adj=1,cex=0.8,col=scol)
		mtext(format(N.all[2,],big.mark=","),side=1,line=8,at=1:10+0.15,cex=0.8,adj=1,col="gray40")
	par(op)

}

# --------------------------------------
# Plot the results, 3 panel figure
# --------------------------------------
pdf("~/dropbox/13beaches/aim1-results/figs/aim1-sens-entero1600-35cfu-length-of-follow-up.pdf",width=7,height=12)
	lo <- layout(mat=matrix(1:3,nrow=3,ncol=1))
	fusensplot(cir.ps,N.ps,"Identified Point Source",labtext=TRUE)
	fusensplot(cir.nps,N.nps,"No Point Source")
	fusensplot(cir.all,N.all,"All Conditions")
dev.off()








