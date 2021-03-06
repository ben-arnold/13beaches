# --------------------------------------
# 8-aim1-descriptive-stat.R
# ben arnold (benarnold@berkeley.edu)
#
# description:
# summarize population characteristics
# overall and at each beach
#
# --------------------------------------


# --------------------------------------
# preamble
# --------------------------------------

rm(list=ls())


# --------------------------------------
# load the analysis dataset
# --------------------------------------

d <- read.csv("~/dropbox/13beaches/data/final/13beaches-analysis.csv")

# convert ID variables from factors to strings
d$hhid <- as.character(d$hhid)
d$indid <- as.character(d$indid)


# --------------------------------------
# reorder and rename the race variable
# --------------------------------------
d$race <- factor(d$race,levels=c("white","non-white, hispanic","white, hispanic","black","asian","american indian","multiple races","other","missing"),labels=c("White/caucasian","Non-White, Hispanic","White, Hispanic","African American","Asian","American Indian","Multiple Races","Other","Missing"))

# combine all hispanic individuals into the same category
d$race[d$race=="White, Hispanic"] <- "Non-White, Hispanic"
d$race <- factor(d$race)
levels(d$race) <- c(levels(d$race)[1],"Hispanic",levels(d$race)[3:8])


# --------------------------------------
# summarize variables
# --------------------------------------


# summarize the data and output it in a
# formatted text matrix

sumtab <- function(d) {
	# create a summary table
	# d is a dataset (possibly subset to different ages)
	
	# Number of individuals
	Ntot <- nrow(d)
	
	# GI illness at enrollment
	Ngibase <- length(d$gibase[d$gibase=="Yes"])
	Pgibase <- 100*mean(ifelse(d$gibase=="Yes",1,0))
	# Individuals at risk of GI illness
	Natrisk <- nrow(d)-Ngibase
	Patrisk <- 100-Pgibase
	# Incident Diarrhea within 3 days (not used - drop to shorten table)
# 	Ngi3 <- sum(d$gici3[d$gibase!="Yes"])
# 	Pgi3 <- 100*mean(d$diarrheaci3[d$gibase!="Yes"])
	# Incident Diarrhea within 10 days
	Ngi10 <- sum(d$diarrheaci10[d$gibase!="Yes"])
	Pgi10 <- 100*mean(d$diarrheaci10[d$gibase!="Yes"])
	
	
	# Age median[iqr]
	medage <- median(d$age,na.rm=TRUE)
	p25age <- quantile(d$age,0.25,na.rm=T)
	p75age <- quantile(d$age,0.75,na.rm=T)
	agesum <- paste(sprintf("%1.0f",medage)," (",sprintf("%1.0f",p25age),",",sprintf("%1.0f",p75age),")",sep="")
	
	# Female
	Nfem <- length(d$sex[d$sex=="Female"])
	Pfem <- 100*mean(ifelse(d$sex=="Female",1,0))
	
	# Race
	Nrace <- table(d$race)
	Prace <- 100*Nrace/sum(Nrace)
	
	# no water contact
	Nnoswim <- length(d$anycontact[d$anycontact=="No"])
	Pnoswim <- 100*mean(ifelse(d$anycontact=="No",1,0))
	# any water contact
	Nanyc <- length(d$anycontact[d$anycontact=="Yes"])
	Panyc <- 100*mean(ifelse(d$anycontact=="Yes",1,0))
	# body immersion
	Nbodyim <- length(d$bodycontact[d$bodycontact=="Yes"])
	Pbodyim <- 100*mean(ifelse(d$bodycontact=="Yes",1,0))
	# head immersion
	Nheadim <- length(d$headunder[d$headunder=="Yes"])
	Pheadim <- 100*mean(ifelse(d$headunder=="Yes",1,0))
	# swallowed water
	Nswallw <- length(d$swallwater[d$swallwater=="Yes"])
	Pswallw <- 100*mean(ifelse(d$swallwater=="Yes",1,0))
	
	# time spent in water among swimmers, hours, median[iqr]
  d$watertime[d$bodycontact=="No" & d$headunder=="No" & d$swallwater=="No"] <- NA
	medwt <- median(d$watertime/60,na.rm=TRUE)
	p25wt <- quantile(d$watertime/60,0.25,na.rm=T)
	p75wt <- quantile(d$watertime/60,0.75,na.rm=T)
	wtsum <- paste(sprintf("%1.1f",medwt)," (",sprintf("%1.1f",p25wt),",",sprintf("%1.1f",p75wt),")",sep="")
	
	# time spent in water, hours, category
	# flag as missing individuals with anycontact="Yes" with no
	# estimate of the time spent in the water
	watcat <- cut(d$watertime/60,c(-1,1,2,3,4,5,24))
	watcat[d$bodycontact=="No" & d$headunder=="No" & d$swallwater=="No"] <- NA
	watcat <- factor(watcat,levels=c(levels(watcat),"Missing"))
	watcat[(d$bodycontact=="Yes"|d$headunder=="Yes"|d$swallwater=="Yes") & is.na(d$watertime)==TRUE] <- "Missing"
	watcat[d$anycontact=="No"] <- NA
	watcat <- factor(watcat,levels=levels(watcat),labels=c("0 -- 1","1.1 -- 2","2.1 -- 3","3.1 -- 4","4.1 -- 5",">5","Missing"))
	Nwatcat <- table(watcat)
	Pwatcat <- 100*Nwatcat/sum(Nwatcat)
	
	
	
	# collate the estimates together
	rowlab <- c("Number of Participants","Gastrointestinal illness at enrollment","Individuals at risk for gastrointestinal illness","Incident diarrhea within 10 days","Age in years", "Female","Race/ethnicity",paste("~~~",levels(d$race),sep=""),	"No water contact","Any water contact","Body immersion","Head immersion","Swallowed water", "Hours spent in the water","Hours spent in the water (cat)",paste("~~~",levels(watcat),sep=""))
	
	Ns <- c(Ntot,Ngibase,Natrisk,Ngi10,NA,Nfem,NA,Nrace,Nnoswim,Nanyc,Nbodyim,Nheadim,Nswallw,NA,NA,Nwatcat)
	Ns <- format(Ns,big.mark=",")
	
	Ps <- c(NA,Pgibase,NA,Pgi10,NA,Pfem,NA,Prace,Pnoswim,Panyc,Pbodyim,Pheadim,Pswallw,NA,NA,Pwatcat)
	Ps <- sprintf("%1.1f",Ps)
	
	Ms <- c(rep(NA,4),agesum,rep(NA,15),wtsum,rep(NA,8))
	
	res <- matrix(c(rowlab,Ns,Ps,Ms), nrow=length(rowlab),ncol=4)

	return(res)
	
}


# calculate summary statistics overall and for each age group
tab.all    <- sumtab(d)
tab.0to4   <- sumtab(d[d$agestrat=="(0, 4]",])
tab.5to10  <- sumtab(d[d$agestrat=="(4, 10]",])
tab.11plus <- sumtab(d[d$agestrat==">10",])



# combine estimates into a single table
tab.res <- cbind(tab.all,tab.0to4[,2:4],tab.5to10[,2:4],tab.11plus[,2:4])
colnames(tab.res) <- c("Var",paste(rep(c("Allages","Age0to4","Age5to10","Age11plus"),c(3,3,3,3)),c("N","%","Median (IQR)")))

# replace NA with ""
for(i in 1:ncol(tab.res)) {
  nai <- grep("NA",tab.res[,i])
  naii <- is.na(tab.res[,i])
  tab.res[nai,i] <- ""
  tab.res[naii==TRUE,i] <-""
}


# --------------------------------------
# output results for formatting in TeX
# --------------------------------------

save(tab.res,file="~/dropbox/13beaches/aim1-results/rawoutput/aim1-descriptive-stat-tables.Rdata")



