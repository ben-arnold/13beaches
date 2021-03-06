# --------------------------------------
# 9-aim1-sens-entero-negcontrol-noswim-figures.R
# ben arnold (benarnold@berkeley.edu)
#
# description:
# plot Enterococcus
# negative control analysis results
#
# --------------------------------------

# --------------------------------------
# input files:
#	13beaches-wq.csv
#	aim1-sens-entero1600-Quartile-regs-noswim.Rdata
#	aim1-sens-entero1600-35cfu-regs-noswim.Rdata
#	aim1-sens-enteroQPCR-Quartile-regs-noswim-negcontrol.Rdata
#	aim1-sens-enteroQPCR-470cce-regs-noswim-negcontrol.Rdata
#
# output files:
#	aim1-sens-entero1600-negcontrol-noswim.pdf
# --------------------------------------

# --------------------------------------
# preamble
# --------------------------------------

rm(list=ls())
library(RColorBrewer)


# --------------------------------------
# Load the water quality dataset to
# get the mid-points of the Entero
# concentrations by each quartile
# --------------------------------------
wq <- read.csv("~/dropbox/13beaches/data/final/13beaches-wq.csv")

# --------------------------------------
# Recode the EPA 1600 values at 
# Mission Bay using the Enterolert
# averages
# --------------------------------------
wq$avgdyentero1600[wq$beach=="Mission Bay"] <- wq$avgdyenteroELT[wq$beach=="Mission Bay"]


minQs <- tapply(wq$avgdyentero1600,wq$qavgdyentero1600,function(x) min(x,na.rm=T))
maxQs <- tapply(wq$avgdyentero1600,wq$qavgdyentero1600,function(x) max(x,na.rm=T))
labQs <- paste(sprintf("%1.0f",10^minQs)," to ",sprintf("%1.0f",10^maxQs),sep="")
labQs <- paste("Q",1:4,"\n(",labQs,")",sep="")

rngQs <- paste(sprintf("%1.0f",round(10^minQs)),"-",sprintf("%1.0f",floor(10^maxQs)),sep="")
rngQs[4] <-paste(">",sprintf("%1.0f",floor(10^maxQs))[3],sep="")

midQs <- minQs + (maxQs-minQs)/2


# --------------------------------------
# open up a larger plot device for a
# figure that combines the quartile
# results and the >35cfu results
# --------------------------------------

pdf("~/dropbox/13beaches/aim1-results/figs/aim1-sens-entero1600-negcontrol-noswim.pdf",width=14,height=5)
lo <- layout(mat=matrix(1:2,nrow=1,ncol=2))


# --------------------------------------
# Load and Plot Entero 1600 Quartile Results
# only plot All Ages because the data
# in the child subgroups are too sparse
# for meaningful estimates. 
# there are too few children who were
# non-swimmers
# --------------------------------------

load("~/dropbox/13beaches/aim1-results/rawoutput/aim1-sens-entero1600-Quartile-regs-noswim.Rdata")



op <- par(mar=c(5,10,7,0)+0.1,xpd=TRUE)
cols <- brewer.pal(9,"YlGnBu")[8:5]
ytics <- seq(0,50,by=10)

# set up an empty plot
MidPts <- barplot(1:1,names.arg=NA,border=NA,col=NA,
	ylim=range(ytics),ylab="",yaxt="n",
	las=1,bty="n"
	)
	segments(x0=0.15,x1=max(MidPts+0.5),y0=ytics,lty=2,lwd=1,col="gray80")
	axis(2,at=ytics,las=1)
	mtext("Diarrhea\nIncidence\nper 1000",side=2,line=3,las=1)
	
	# calculate X coordinates relative to the mid points for each group
	xspan <- 0.37
	xplus <- c(-xspan, -xspan/3, xspan/3, xspan)  # evenly distribute 4 datapoints around each midpoint
	
	
	xall  <- xplus+MidPts[1]  # for table and quartile labels in header/footer
	allxs <- xall
	labx  <- MidPts[1]-xspan*1.5  # for left-hand labels in the header/footer

	# plot all age estimates
	segments(x0=allxs,y0=ci.all[,"CIlb"]*1000,y1=ci.all[,"CIub"]*1000,lwd=2,col=cols)
	points(xall,ci.all[,"CI"]*1000,pch=16,cex=1.75,lwd=2,col=cols)
	text(x=xall,y=ci.all[,"CI"]*1000,sprintf("%1.0f",ci.all[,"CI"]*1000),pos=4,cex=0.7,col=cols,font=2)
	
	
	# print header labels
	mtext(mtext(expression(paste(italic("Enterococcus")," Quartile")),side=3,line=5.5,cex=1))
	mtext(expression(paste(italic("Enterococcus")," Quartile")),side=3,line=4,at=labx,adj=1,col=cols[3],cex=0.8)
	mtext(c("Q1","Q2","Q3","Q4"),side=3,line=4,at=allxs,col=cols,cex=0.9,font=2)
	
	mtext("Range (CFU/100ml)",side=3,line=3,at=labx,adj=1,col=cols[3],cex=0.8)
	mtext(rngQs,side=3,line=3,at=allxs,col=cols,cex=0.8)
	
	# Print adjusted CIRs and 95% CIs (formatted)
	cirform <- function(cirs) {
		paste(sprintf("%1.2f",cirs),sep="")
	}
	circiform <- function(circi) {
		apply(circi,1,function(x) paste("(",sprintf("%1.2f",x["CIRlb"]),", ",sprintf("%1.2f",x["CIRub"]),")",sep="") )
	}
	mtext("Adjusted CIR",side=3,line=2,at=labx,adj=1,cex=0.8,col="gray30")
	mtext(c("ref",cirform(cir.all[,1])),side=3,line=2,at=allxs,cex=0.75)
	
	mtext("(95% CI)",side=3,line=1,at=labx,adj=1,cex=0.8,col="gray30")
	mtext(c("",circiform(cir.all[,2:3])),side=3,line=1,at=allxs,cex=0.7)
	

	# print footer labels
	mtext(c("Q1","Q2","Q3","Q4"),side=1,line=0.25,at=allxs,col=cols,cex=0.9,font=2)
	
	# print table with Ns
	mtext("Incident Diarrhea Cases",side=1,line=2.25,at=labx,adj=1,cex=0.8,col="gray30")
	ns <- c(N.all[,1])
	mtext(  format(ns,big.mark=","),side=1,line=2.25,at=allxs+0.03,adj=1,cex=0.75    )
	
	mtext("Population At Risk",side=1,line=3.25,at=labx,adj=1,cex=0.8,col="gray30")
	Ns <- c(N.all[,2])
	mtext(  format(Ns,big.mark=","),side=1,line=3.25,at=allxs+0.03,adj=1,cex=0.75    )


# make an overall title (commented out)
# mtext("Non-Swimmer Negative Control Analysis",side=3,line=8,cex=2,at=max(MidPts),adj=0)

# --------------------------------------
# Load and Plot Entero 1600 35 cfu Results
# only plot All Ages because the data
# in the child subgroups are too sparse
# for meaningful estimates. 
# there are too few children who were
# non-swimmers
# --------------------------------------

load("~/dropbox/13beaches/aim1-results/rawoutput/aim1-sens-entero1600-35cfu-regs-noswim.Rdata")


cols <- brewer.pal(9,"YlGn")[c(8,6)]
ytics <- seq(0,50,by=10)
# set up an empty plot
MidPts <- barplot(1:1,names.arg=NA,border=NA,col=NA,
	ylim=range(ytics),ylab="",yaxt="n",
	las=1,bty="n"
)
segments(x0=0.15,x1=max(MidPts+0.5),y0=ytics,lty=2,lwd=1,col="gray80")
axis(2,at=ytics,las=1)
mtext("Diarrhea\nIncidence\nper 1000",side=2,line=3,las=1)

# calculate X coordinates relative to the mid points for each group
xspan <- 0.2
xplus <- c(-xspan, xspan)  # evenly distribute 2 datapoints around each midpoint

xall    <- xplus+MidPts[1]
allxs   <- xall 
labx <- MidPts[1]-xspan*3  # for left-hand labels in the header/footer


# plot all conditions estimates
segments(x0=xall,y0=ci.all[,"CIlb"]*1000,y1=ci.all[,"CIub"]*1000,lwd=2,col=cols)
points(xall,ci.all[,"CI"]*1000,pch=16,cex=1.5,lwd=2,col=cols)
text(x=xall,y=ci.all[,"CI"]*1000,sprintf("%1.0f",ci.all[,"CI"]*1000),pos=4,cex=0.7,col=cols,font=2)

# print header labels
mtext(mtext(expression(paste(italic("Enterococcus")," >35 CFU/100ml")),side=3,line=5.5,cex=1  ))

mtext(expression(paste(italic("Enterococcus")," Category")),side=3,line=4,at=labx,adj=1,col=cols[2],cex=0.8)
mtext(expression(""<= 35,""> 35),side=3,line=4,at=allxs,col=cols,cex=1,font=2)
mtext("CFU/100ml",side=3,line=3,at=allxs,col=cols,cex=0.75)


# Print adjusted CIRs and 95% CIs (formatted)
cirform <- function(cirs) {
	paste(sprintf("%1.2f",cirs),sep="")
}
circiform <- function(circi) {
	paste("(",sprintf("%1.2f",circi[1]),", ",sprintf("%1.2f",circi[2]),")",sep="")
}
mtext("Adjusted CIR",side=3,line=2,at=labx,adj=1,cex=0.8,col="gray30")
mtext(c("ref",cirform(cir.all[1])),side=3,line=2,at=allxs,cex=0.75)

mtext("(95% CI)",side=3,line=1,at=labx,adj=1,cex=0.8,col="gray30")
mtext(c("",circiform(cir.all[2:3])),side=3,line=1,at=allxs,cex=0.7)


# print footer labels
mtext(expression(""<= 35,""> 35),side=1,line=0.25,at=allxs,col=cols,cex=1,font=2)
mtext("CFU/100ml",side=1,line=1.25,at=allxs,col=cols,cex=0.75)

# print table with Ns
mtext("Incident Diarrhea Cases",side=1,line=2.25,at=labx,adj=1,cex=0.8,col="gray30")
ns <- c(N.all[,1])
mtext(  format(ns,big.mark=","),side=1,line=2.25,at=allxs+0.06,adj=1,cex=0.75    )

mtext("Population At Risk",side=1,line=3.25,at=labx,adj=1,cex=0.8,col="gray30")
Ns <- c(N.all[,2])
mtext(  format(Ns,big.mark=","),side=1,line=3.25,at=allxs+0.06,adj=1,cex=0.75    )

par(op)
dev.off()



# --------------------------------------
# Load the water quality dataset to
# get the ranges of the Entero qPCR
# quartiles
# --------------------------------------

minQs <- tapply(wq$avgdyenteropcr,wq$qavgdyenteropcr,function(x) min(x,na.rm=T))
maxQs <- tapply(wq$avgdyenteropcr,wq$qavgdyenteropcr,function(x) max(x,na.rm=T))
labQs <- paste(sprintf("%1.0f",10^minQs)," to ",sprintf("%1.0f",10^maxQs),sep="")
labQs <- paste("Q",1:4,"\n(",labQs,")",sep="")

rngQs <- paste(sprintf("%1.0f",round(10^minQs)),"-",sprintf("%1.0f",floor(10^maxQs)),sep="")
rngQs[4] <-paste(">",sprintf("%1.0f",floor(10^maxQs))[3],sep="")

midQs <- minQs + (maxQs-minQs)/2


# --------------------------------------
# open up a larger plot device for a
# figure that combines the quartile
# results and the >470cce results
# --------------------------------------

pdf("~/dropbox/13beaches/aim1-results/figs/aim1-sens-enteroQPCR-negcontrol-noswim.pdf",width=14,height=5)
lo <- layout(mat=matrix(1:2,nrow=1,ncol=2))


# --------------------------------------
# Load and Plot Entero qPCR Quartile Results
# only plot All Ages because the data
# in the child subgroups are too sparse
# for meaningful estimates. 
# there are too few children who were
# non-swimmers
# --------------------------------------

load("~/dropbox/13beaches/aim1-results/rawoutput/aim1-sens-enteroQPCR-Quartile-regs-noswim-negcontrol.Rdata")



op <- par(mar=c(5,10,7,0)+0.1,xpd=TRUE)
cols <- brewer.pal(9,"YlOrBr")[8:5]
ytics <- seq(0,60,by=10)

# set up an empty plot
MidPts <- barplot(1:1,names.arg=NA,border=NA,col=NA,
	ylim=range(ytics),ylab="",yaxt="n",
	las=1,bty="n"
	)
	segments(x0=0.15,x1=max(MidPts+0.5),y0=ytics,lty=2,lwd=1,col="gray80")
	axis(2,at=ytics,las=1)
	mtext("Diarrhea\nIncidence\nper 1000",side=2,line=3,las=1)
	
	# calculate X coordinates relative to the mid points for each group
	xspan <- 0.37
	xplus <- c(-xspan, -xspan/3, xspan/3, xspan)  # evenly distribute 4 datapoints around each midpoint
	
	
	xall  <- xplus+MidPts[1]  # for table and quartile labels in header/footer
	allxs <- xall
	labx  <- MidPts[1]-xspan*1.5  # for left-hand labels in the header/footer

	# plot all age estimates
	segments(x0=allxs,y0=ci.all[,"CIlb"]*1000,y1=ci.all[,"CIub"]*1000,lwd=2,col=cols)
	points(xall,ci.all[,"CI"]*1000,pch=16,cex=1.75,lwd=2,col=cols)
	text(x=xall,y=ci.all[,"CI"]*1000,sprintf("%1.0f",ci.all[,"CI"]*1000),pos=4,cex=0.7,col=cols,font=2)
	
	
	# print header labels
	mtext(mtext(expression(paste(italic("Enterococcus")," qPCR Quartile")),side=3,line=5.5,cex=1))
	mtext(expression(paste(italic("Enterococcus")," qPCR Quartile")),side=3,line=4,at=labx,adj=1,col=cols[3],cex=0.8)
	mtext(c("Q1","Q2","Q3","Q4"),side=3,line=4,at=allxs,col=cols,cex=0.9,font=2)
	
	mtext("Range (CCE/100ml)",side=3,line=3,at=labx,adj=1,col=cols[3],cex=0.8)
	mtext(rngQs,side=3,line=3,at=allxs,col=cols,cex=0.8)
	
	# Print adjusted CIRs and 95% CIs (formatted)
	cirform <- function(cirs) {
		paste(sprintf("%1.2f",cirs),sep="")
	}
	circiform <- function(circi) {
		apply(circi,1,function(x) paste("(",sprintf("%1.2f",x["CIRlb"]),", ",sprintf("%1.2f",x["CIRub"]),")",sep="") )
	}
	mtext("Adjusted CIR",side=3,line=2,at=labx,adj=1,cex=0.8,col="gray30")
	mtext(c("ref",cirform(cir.all[,1])),side=3,line=2,at=allxs,cex=0.75)
	
	mtext("(95% CI)",side=3,line=1,at=labx,adj=1,cex=0.8,col="gray30")
	mtext(c("",circiform(cir.all[,2:3])),side=3,line=1,at=allxs,cex=0.7)
	

	# print footer labels
	mtext(c("Q1","Q2","Q3","Q4"),side=1,line=0.25,at=allxs,col=cols,cex=0.9,font=2)
	
	# print table with Ns
	mtext("Incident Diarrhea Cases",side=1,line=2.25,at=labx,adj=1,cex=0.8,col="gray30")
	ns <- c(N.all[,1])
	mtext(  format(ns,big.mark=","),side=1,line=2.25,at=allxs+0.03,adj=1,cex=0.75    )
	
	mtext("Population At Risk",side=1,line=3.25,at=labx,adj=1,cex=0.8,col="gray30")
	Ns <- c(N.all[,2])
	mtext(  format(Ns,big.mark=","),side=1,line=3.25,at=allxs+0.03,adj=1,cex=0.75    )


# make an overall title (commented out)
# mtext("Non-Swimmer Negative Control Analysis",side=3,line=8,cex=2,at=max(MidPts),adj=0)

# --------------------------------------
# Load and Plot Entero 1600 35 cfu Results
# only plot All Ages because the data
# in the child subgroups are too sparse
# for meaningful estimates. 
# there are too few children who were
# non-swimmers
# --------------------------------------

load("~/dropbox/13beaches/aim1-results/rawoutput/aim1-sens-enteroQPCR-470cce-regs-noswim-negcontrol.Rdata")


cols <- brewer.pal(9,"Purples")[c(8,6)]
ytics <- seq(0,60,by=10)
# set up an empty plot
MidPts <- barplot(1:1,names.arg=NA,border=NA,col=NA,
	ylim=range(ytics),ylab="",yaxt="n",
	las=1,bty="n"
)
segments(x0=0.15,x1=max(MidPts+0.5),y0=ytics,lty=2,lwd=1,col="gray80")
axis(2,at=ytics,las=1)
mtext("Diarrhea\nIncidence\nper 1000",side=2,line=3,las=1)

# calculate X coordinates relative to the mid points for each group
xspan <- 0.2
xplus <- c(-xspan, xspan)  # evenly distribute 2 datapoints around each midpoint

xall    <- xplus+MidPts[1]
allxs   <- xall 
labx <- MidPts[1]-xspan*3  # for left-hand labels in the header/footer


# plot all conditions estimates
segments(x0=xall,y0=ci.all[,"CIlb"]*1000,y1=ci.all[,"CIub"]*1000,lwd=2,col=cols)
points(xall,ci.all[,"CI"]*1000,pch=16,cex=1.5,lwd=2,col=cols)
text(x=xall,y=ci.all[,"CI"]*1000,sprintf("%1.0f",ci.all[,"CI"]*1000),pos=4,cex=0.7,col=cols,font=2)

# print header labels
mtext(mtext(expression(paste(italic("Enterococcus")," qPCR >470 CCE/100ml")),side=3,line=5.5,cex=1  ))

mtext(expression(paste(italic("Enterococcus")," qPCR Category")),side=3,line=4,at=labx,adj=1,col=cols[2],cex=0.8)
mtext(expression(""<= 470,""> 470),side=3,line=4,at=allxs,col=cols,cex=1,font=2)
mtext("CCE/100ml",side=3,line=3,at=allxs,col=cols,cex=0.75)


# Print adjusted CIRs and 95% CIs (formatted)
cirform <- function(cirs) {
	paste(sprintf("%1.2f",cirs),sep="")
}
circiform <- function(circi) {
	paste("(",sprintf("%1.2f",circi[1]),", ",sprintf("%1.2f",circi[2]),")",sep="")
}
mtext("Adjusted CIR",side=3,line=2,at=labx,adj=1,cex=0.8,col="gray30")
mtext(c("ref",cirform(cir.all[1])),side=3,line=2,at=allxs,cex=0.75)

mtext("(95% CI)",side=3,line=1,at=labx,adj=1,cex=0.8,col="gray30")
mtext(c("",circiform(cir.all[2:3])),side=3,line=1,at=allxs,cex=0.7)


# print footer labels
mtext(expression(""<= 470,""> 470),side=1,line=0.25,at=allxs,col=cols,cex=1,font=2)
mtext("CCE/100ml",side=1,line=1.25,at=allxs,col=cols,cex=0.75)

# print table with Ns
mtext("Incident Diarrhea Cases",side=1,line=2.25,at=labx,adj=1,cex=0.8,col="gray30")
ns <- c(N.all[,1])
mtext(  format(ns,big.mark=","),side=1,line=2.25,at=allxs+0.06,adj=1,cex=0.75    )

mtext("Population At Risk",side=1,line=3.25,at=labx,adj=1,cex=0.8,col="gray30")
Ns <- c(N.all[,2])
mtext(  format(Ns,big.mark=","),side=1,line=3.25,at=allxs+0.06,adj=1,cex=0.75    )

par(op)
dev.off()





