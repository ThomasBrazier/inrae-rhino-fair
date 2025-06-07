#===============================#
#            COLONY
#     Effective dispersal kernel
#===============================#

wd="C:/Users/tbrazier/Desktop/R/Distance kernel"
wd="/Volumes/TRANSCEND16/R/Distance kernel"
wd="F:/R/Distance kernel"
setwd(wd)


#---------------------------------#
# EXPONENTIAL POWER
#---------------------------------
#representer le dispersal kernel à partir de la distance de dispersion moyenne
x=seq(0,100,by=0.01)
d=11.3 # distance de dispersion moyenne
b=0.5
tau=1
a=11.3/2


library(distr)
#Exponential
DistExp=AbscontDistribution(d = function( x ){(x/(a^2))*exp(-x/(a))},low1=0,up1=100)
plot(DistExp)

#Probability for 80
p(DistExp)(60)
p(DistExp)(80)


par(mfrow=c(1,1))
y=d(DistExp)(x)
plot(y~x,type="l",axes=F,xlab="",ylab="")
axis(1)
axis(2)
title(xlab="Dispersal distance (km)", ylab="Relative probability of dispersal event")

#line obtained from fitdistrplus with fairon data
y2=function( x ){0.3777125*exp(-1*0.3777125*x)}
curve(y2,add=F,from=0,to=100,xlab="",ylab="")
axis(1)
axis(2)
title(xlab="Dispersal distance (km)", ylab="Probability of dispersal event")



#---------------------------------#
# POWER LAW
#---------------------------------

#---------------------------------#
# d'apres Script PL fourni par Eric Petit - NEGATIVE EXPONENTIAL
#---------------------------------
#representer le dispersal kernel à partir de la distance de dispersion moyenne
x=seq(0,100,by=0.01)
# negative exponential
a=11.3/2

library(distr)
#Exponential
DistExp=AbscontDistribution(d = function( x ){(x/(a^2))*exp(-x/(a))},low1=0,up1=100)
plot(DistExp)

#Probability for 80
p(DistExp)(60)
p(DistExp)(80)


par(mfrow=c(1,1))
y=d(DistExp)(x)
plot(y~x,type="l",axes=F,xlab="",ylab="")
axis(1)
axis(2)
title(xlab="Dispersal distance (km)", ylab="Relative probability of dispersal event")

#line obtained from fitdistrplus with fairon data
y2=function( x ){0.3777125*exp(-1*0.3777125*x)}
curve(y2,add=F,from=0,to=100,xlab="",ylab="")
axis(1)
axis(2)
title(xlab="Dispersal distance (km)", ylab="Probability of dispersal event")
