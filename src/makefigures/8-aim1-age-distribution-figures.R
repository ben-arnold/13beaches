



# --------------------------------------
# 8-aim1-age-distribution-figures.R
# ben arnold (benarnold@berkeley.edu)
#
# description:
# plot the distributions of age in the
# cohorts
#
# --------------------------------------

# --------------------------------------
# input files:
#    13beaches-analysis.csv
#
# output files:
#    aim1-age-distributions.pdf
# --------------------------------------


# --------------------------------------
# preamble
# --------------------------------------

rm(list=ls())
library(RColorBrewer)
library(scales)


# --------------------------------------
# load the individual participant dataset
# --------------------------------------

d <- read.csv("~/dropbox/13beaches/data/final/13beaches-analysis.csv")



# --------------------------------------
# make the plots for water quality
# samples
# --------------------------------------

pdf("~/dropbox/13beaches/aim1-results/figs/aim1-age-distributions.pdf",width=6,height=12)

# common plot parameters and layout
xtics <- seq(0,100,by=10)
histbreaks <- seq(0,max(d$age,na.rm=T)+1,by=1)
ytics <- seq(0,2000,by=200)
cols <- c(brewer.pal(9,"YlGnBu")[7], brewer.pal(9,"YlOrBr")[4]   )
lo <- layout(mat=matrix(1:3,nrow=3,ncol=1))

# All People
hist(d$age,breaks=histbreaks,
	col="white",border="gray50",
	xlab="Age, years",xlim=range(xtics),xaxt="n",
	ylab="",ylim=range(ytics),yaxt="n",
	main=""
	)
	axis(1,at=xtics)
	axis(2,at=ytics,las=1)
	mtext("Frequency",side=3,line=0,at=-15,adj=0,cex=0.7)
	mtext("All Beachgoers",side=3,line=2,cex=1.2,font=2,adj=0)
	mtext("A)",side=3,line=2,at=-10,font=2,cex=1.5)

# Non-swimmers
hist(d$age[d$anycontact=="No"],breaks=histbreaks,
	col="white",border="gray50",
	xlab="Age, years",xlim=range(xtics),xaxt="n",
	ylab="",ylim=range(ytics),yaxt="n",
	main=""
	)
	axis(1,at=xtics)
	axis(2,at=ytics,las=1)
	mtext("Frequency",side=3,line=0,at=-15,adj=0,cex=0.7)
	mtext("Non-Swimmers",side=3,line=2,cex=1.2,font=2,adj=0)
	mtext("B)",side=3,line=2,at=-10,font=2,cex=1.5)

# Body Immersion Swimmers
hist(d$age[d$bodycontact=="Yes"],breaks=histbreaks,
	col="white",border="gray50",
	xlab="Age, years",xlim=range(xtics),xaxt="n",
	ylab="",ylim=range(ytics),yaxt="n",
	main=""
	)
	axis(1,at=xtics)
	axis(2,at=ytics,las=1)
	mtext("Frequency",side=3,line=0,at=-15,adj=0,cex=0.7)
	mtext("Body Immersion Swimmers",side=3,line=2,cex=1.2,font=2,adj=0)
	mtext("C)",side=3,line=2,at=-10,font=2,cex=1.5)
	
	
dev.off()




