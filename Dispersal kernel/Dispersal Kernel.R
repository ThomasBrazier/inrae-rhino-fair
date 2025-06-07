#===============================#
#            COLONY
#     Effective dispersal kernel
#===============================#

# DIRECTORIES
# Working directories
# New users need to configure correct pathways
if (Sys.info()["sysname"]=="Darwin") {
  wd="/Volumes/Samsung_T5/INRA RHINO/R/Dispersal kernel/"
  setwd(wd)
  datadir="/Volumes/Samsung_T5/INRA RHINO/R/Assignation/04_data/outputs/"
  datadirThu="/Volumes/Samsung_T5/INRA RHINO/R/AssignationThu/04_data/outputs/"
  graphdir="/Volumes/Samsung_T5/INRA RHINO/R/Graph/"
  Rdir="/Volumes/Samsung_T5/INRA RHINO/R/"
}else {
  if (Sys.info()["sysname"]=="Windows"){
    wd="E:/INRA RHINO/R/Dispersal kernel/"
    setwd(wd)
    datadir="E:/INRA RHINO/R/Assignation/04_data/outputs/"
    datadirThu="E:/INRA RHINO/R/AssignationThu/04_data/outputs/"
    graphdir="E:/INRA RHINO/R/Graph/"
    Rdir="E:/INRA RHINO/R/"
  }
}
#---------------------------------#
# HISTOGRAM OF OBSERVED DISPERSAL DISTANCES
# AND FITTED DISPERSAL DISTANCE KERNEL
#---------------------------------

# RELATIVE PROBABILITY OF DISPERSAL EVENT ~ DISPERSAL DISTANCE
# df=read.table("distancesDist.txt",header=TRUE)
# dfExp=df[which(df$dist=="H0"),]
# dfObs=df[which(df$dist=="Obs"),]
# # Vecteur de valeurs : freq relative pour chaque distance
# freqRel=data.frame(f=rep(NA,66),d=seq(0,65,1))
# for (i in 1:nrow(freqRel)) {
#   freqRel[i,1]=length(which(round(dfObs$d)==freqRel[i,2]))/nrow(dfObs)
# }
# 
# ggplot(data=dfObs, aes(x=d)) +
#   geom_histogram(aes(y = (..count..)/sum(..count..)),binwidth=1,size=4)+
#   stat_function(fun=dexp,geom="line",size=1,args = (mean=0.49))+
#   ggtitle(paste("Distribution of offspring-father distances\n(n=",length(dyadsObs$distance),")"))+
#   xlab("Geographic distance (km)") + ylab("Dyads relative frequency")+
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank(),
#         plot.title = element_text(color="black", size=34, face="bold.italic",hjust = 0.5),
#         axis.title.x = element_text(color="black", size=30),
#         axis.title.y = element_text(color="black", size=30),
#         axis.text=element_text(size=26, colour="black"),
#         legend.key = element_rect(fill = "white", size = 1),
#         legend.key.height = unit(2,"line"),
#         legend.key.width = unit(5,"line"),
#         legend.text=element_text(size=30),
#         legend.title=element_text(size=30),
#         legend.position='none')

#---------------------------------#
# EXPONENTIAL DISTRIBUTION WITH PARAMETER ESTIMATED FROM MASTERBAYES
#---------------------------------

x=seq(0,100,by=0.01)

# Probability density function
lambda=c(0.8413,0.4704699,0.2828)
y2=function( x ){lambda[2]*exp(-1*lambda[2]*x)}
y1=function( x ){lambda[1]*exp(-1*lambda[1]*x)}
y3=function( x ){lambda[3]*exp(-1*lambda[3]*x)}
curve(y2,add=F,from=0,to=60,xlab="",ylab="")
curve(y1,add=T,from=0,to=60,xlab="",ylab="",lty="dashed")
curve(y3,add=T,from=0,to=60,xlab="",ylab="",lty="dashed")
axis(1)
axis(2)
title(xlab="Dispersal distance (km)", ylab="Probability of dispersal event")

# Cumulative density function
lambda=0.682301
y2=function( x ){1-exp(-1*lambda*x)}
curve(y2,add=F,from=0,to=100,xlab="",ylab="")
axis(1)
axis(2)
title(xlab="Dispersal distance (km)", ylab="Cumulative probability of dispersal event")

library(distr)
#Exponential
a=0.4704699
a=0.1003581
DistExp=AbscontDistribution(d = function( x ){(x/(a^2))*exp(-x/(a))},low1=0,up1=100)
plot(DistExp)
par(mfrow=c(1,1))


##############################################
#       PICARDY
##############################################
#---------------------------------#
# DISTRIBUTION OF DISTANCES AND FITTED DISPERSAL KERNEL
#---------------------------------

dyadsObsSelect=read.table(paste(datadir,"dyadsObsSelect.txt",sep=""),header=TRUE)
# Uses the fitdistRplus package
library(fitdistrplus)
# Two methods can be used :
# 1/ Maximum likelihood
# 2/ maximizing goodness-of-fit estimation
fit.exp=fitdist(data=dyadsObsSelect$distance,distr="exp",method="mle")
fit.gamma=fitdist(data=dyadsObsSelect$distance,distr="gamma",method="mle",lower = c(0, 0),
                  start = list(shape = 1, scale = 1))
fit.weibull=fitdist(data=dyadsObsSelect$distance,distr="weibull",method="mle",lower = c(0, 0),
                    start = list(shape = 1, scale = 1))

denscomp(list(fit.exp,fit.gamma,fit.weibull))
gofstat(list(fit.exp,fit.gamma,fit.weibull))
# Bootstrap procedure for estimates
summary(bootdist(fit.exp,bootmethod = "nonparam", niter = 1000))
summary(bootdist(fit.gamma,bootmethod = "nonparam", niter = 1000))
summary(bootdist(fit.weibull,bootmethod = "nonparam", niter = 1000))



##############################################
#       THURINGIA
##############################################
#---------------------------------#
# DISTRIBUTION OF DISTANCES AND FITTED DISPERSAL KERNEL
#---------------------------------

dyadsObsSelect=read.table(paste(datadirThu,"dyadsObsSelect.txt",sep=""),header=TRUE)
# Uses the fitdistRplus package
library(fitdistrplus)
# Two methods can be used :
# 1/ Maximum likelihood
# 2/ maximizing goodness-of-fit estimation
fit.exp=fitdist(data=dyadsObsSelect$distance,distr="exp",method="mle")
fit.gamma=fitdist(data=dyadsObsSelect$distance,distr="gamma",method="mle",lower = c(0, 0),
                  start = list(shape = 1, scale = 1))
fit.weibull=fitdist(data=dyadsObsSelect$distance,distr="weibull",method="mle",lower = c(0, 0),
                    start = list(shape = 1, scale = 1))

denscomp(list(fit.exp,fit.gamma,fit.weibull))
gofstat(list(fit.exp,fit.gamma,fit.weibull))
# Bootstrap procedure for estimates
summary(bootdist(fit.exp,bootmethod = "nonparam", niter = 1000))
summary(bootdist(fit.gamma,bootmethod = "nonparam", niter = 1000))
summary(bootdist(fit.weibull,bootmethod = "nonparam", niter = 1000))




#---------------------------------#
# POWER LAW
#---------------------------------

#---------------------------------#
# d'apres Script PL fourni par Eric Petit - NEGATIVE EXPONENTIAL
#---------------------------------
# #representer le dispersal kernel à partir de la distance de dispersion moyenne
# x=seq(0,100,by=0.01)
# # negative exponential
# a=11.3/2
# 
# library(distr)
# #Exponential
# DistExp=AbscontDistribution(d = function( x ){(x/(a^2))*exp(-x/(a))},low1=0,up1=100)
# # plot(DistExp)
# 
# #Probability for 80
# p(DistExp)(60)
# p(DistExp)(80)
# 
# 
# par(mfrow=c(1,1))
# y=d(DistExp)(x)
# plot(y~x,type="l",axes=F,xlab="",ylab="")
# axis(1)
# axis(2)
# title(xlab="Dispersal distance (km)", ylab="Relative probability of dispersal event")
# 
# #line obtained from fitdistrplus with fairon data
# y2=function( x ){0.3777125*exp(-1*0.3777125*x)}
# curve(y2,add=F,from=0,to=100,xlab="",ylab="")
# axis(1)
# axis(2)
# title(xlab="Dispersal distance (km)", ylab="Probability of dispersal event")

