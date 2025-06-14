#===============================#
#            COLONY
#     Effective dispersal kernel
#===============================#

# Need the distancesDist.txt output file as dataset

# clear global environment: remove all variables
rm(list=ls(all=TRUE))
# library(rstudioapi)
# Get the directory of the file & set working directory
# filedir=dirname(rstudioapi::getSourceEditorContext()$path)
# setwd(filedir)

#----------------------------------------------------------#
# Loading packages
# Check if packages are installed, install if necessary
source("Sources/packages.R")



#=========================================================
# Mean dispersal distance
#=========================================================
##############################################
#       PICARDY
##############################################
#---------------------------------#
# DISTRIBUTION OF DISTANCES AND FITTED DISPERSAL KERNEL
#---------------------------------

dyadsObsSelect = read.table("AssignationPic/dyadsObsSelect_gametic.txt", header=TRUE)

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
mean(unlist(bootdist(fit.exp,bootmethod = "nonparam", niter = 1000)$estim))
summary(bootdist(fit.exp,bootmethod = "nonparam", niter = 1000))
summary(bootdist(fit.gamma,bootmethod = "nonparam", niter = 1000))
summary(bootdist(fit.weibull,bootmethod = "nonparam", niter = 1000))



##############################################
#       THURINGIA
##############################################
#---------------------------------#
# DISTRIBUTION OF DISTANCES AND FITTED DISPERSAL KERNEL
#---------------------------------

dyadsObsSelect=read.table(paste("Data/Thu/dyadsObsSelect.txt",sep=""),header=TRUE)
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
mean(unlist(bootdist(fit.exp,bootmethod = "nonparam", niter = 1000)$estim))
summary(bootdist(fit.exp,bootmethod = "nonparam", niter = 1000))
summary(bootdist(fit.gamma,bootmethod = "nonparam", niter = 1000))
summary(bootdist(fit.weibull,bootmethod = "nonparam", niter = 1000))


#---------------------------------#
# MOMENTS OF DISPERSAL KERNELS
#---------------------------------
# A kurtosis > 3 means leptokurtosis in the distribution, with excess of short/long distances over intermediate distances

# "The kurtosis of any univariate normal distribution is 3. It is common to compare the kurtosis of a distribution to this value.
# Distributions with kurtosis less than 3 are said to be platykurtic, although this does not imply the distribution is "flat-topped" as is sometimes stated.
# Rather, it means the distribution produces fewer and less extreme outliers than does the normal distribution.
# An example of a platykurtic distribution is the uniform distribution, which does not produce outliers. Distributions with kurtosis greater than 3 are said to be leptokurtic.
# An example of a leptokurtic distribution is the Laplace distribution, which has tails that asymptotically approach zero more slowly than a Gaussian,
# and therefore produces more outliers than the normal distribution.
# It is also common practice to use an adjusted version of Pearson's kurtosis, the excess kurtosis, which is the kurtosis minus 3,
# to provide the comparison to the normal distribution. Some authors use "kurtosis" by itself to refer to the excess kurtosis.
# For clarity and generality, however, this article follows the non-excess convention and explicitly indicates where excess kurtosis is meant. " (Wikipedia)
# distPic=read.table("AssignationPic/outputs/distancesDist.txt",header=TRUE)
# distThu=read.table("AssignationThu/outputs/distancesDist.txt",header=TRUE)

# # Selecting data
# Pic = distPic$d[distPic$d > 0 & distPic$dist == "ObsSelect"]
# Thu = distThu$d[distThu$d > 0 & distThu$dist == "ObsSelect"]
# # Figures
# hist(Pic, breaks = 20)
# hist(Thu, breaks = 20)
# 
# #---------------------------------------------------------
# # SKEWNESS
# #---------------------------------------------------------
# ?skewness
# # Skewness estimated with 'moments' library
# skewness(Pic)
# skewness(Thu)
# 
# # Estimate skewness with resampling
# nboot = 1000
# boot = numeric(nboot)
# for (i in 1:nboot) {
#   boot[i] = skewness(sample(Pic, replace = TRUE))
# }
# hist(boot, breaks = 20)
# mean(boot)
# quantile(boot, 0.025)
# quantile(boot, 0.975)
# 
# 
# nboot = 1000
# boot = numeric(nboot)
# for (i in 1:nboot) {
#   boot[i] = skewness(sample(Thu, replace = TRUE))
# }
# hist(boot, breaks = 20)
# mean(boot)
# quantile(boot, 0.025)
# quantile(boot, 0.975)
# 
# #---------------------------------------------------------
# # KURTOSIS
# #---------------------------------------------------------
# ?kurtosis
# # Kurtosis estimated with 'moments' library
# kurtosis(Pic)
# kurtosis(Thu)
# # Hence excess of kurtosis is simply kurtosis - 3
# kurtosis(Pic) - 3
# kurtosis(Thu) - 3
# 
# # Estimate excess kurtosis with resampling
# nboot = 1000
# boot = numeric(nboot)
# for (i in 1:nboot) {
#   boot[i] = kurtosis(sample(Pic, replace = TRUE)) - 3
# }
# hist(boot, breaks = 20)
# mean(boot)
# quantile(boot, 0.025)
# quantile(boot, 0.975)
# 
# nboot = 1000
# boot = numeric(nboot)
# for (i in 1:nboot) {
#   boot[i] = kurtosis(sample(Thu, replace = TRUE)) - 3
# }
# hist(boot, breaks = 20)
# mean(boot)
# quantile(boot, 0.025)
# quantile(boot, 0.975)


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



#===============================#
#            NATAL DISPERSAL
#===============================#
dataPic = read.table("Data/Pic/uniqueGenotypesWithInfo.txt", header = TRUE)
dataThu = read.table("Data/Thu/uniqueGenotypesWithInfo.txt", header = TRUE)

# Where are juveniles sampled in next years?
# If only one colony sampled multiple times -> no natal dispersal
# If more than one colony sampled -> natal dispersal
# We can make the ratio, i.e. the natal dispersal rate
dataPic$idcol[dataPic$ageWhenFirstCaptur == "Juv" & dataPic$sexe == "M"]

apply(dataPic[dataPic$ageWhenFirstCaptur == "Juv" & dataPic$sexe == "M", 22:29], 1, function(x) length(x[!is.na(x)]))
(ntimesampled = sum(apply(dataPic[dataPic$ageWhenFirstCaptur == "Juv" & dataPic$sexe == "M", 22:29], 1, function(x) length(x[!is.na(x)])) > 1))
# 307 juveniles have been sampled more than one time
# involving 56 males

sapply(apply(dataPic[dataPic$ageWhenFirstCaptur == "Juv" & dataPic$sexe == "M", 22:29], 1, function(x) unique(x[!is.na(x)])), length)
(ncolsampled = sum(sapply(apply(dataPic[dataPic$ageWhenFirstCaptur == "Juv" & dataPic$sexe == "M", 22:29], 1, function(x) unique(x[!is.na(x)])), length) > 1))
# 18 juveniles have been sampled in more than one colony
# involving 2 males

(nataldisprate = ncolsampled/ntimesampled)
# Natal dispersal rate in Picardy is very low = 0.0586
# For male only, natal dispersal rate = 0.036




dataThu$idcol[dataThu$ageWhenFirstCaptur == "Juv" & dataThu$sexe == "M"]

apply(dataThu[dataThu$ageWhenFirstCaptur == "Juv" & dataThu$sexe == "M", 22:27], 1, function(x) length(x[!is.na(x)]))
(ntimesampled = sum(apply(dataThu[dataThu$ageWhenFirstCaptur == "Juv" & dataThu$sexe == "M", 22:27], 1, function(x) length(x[!is.na(x)])) > 1))
# 475 juveniles have been sampled more than one time
# involving 146 males

sapply(apply(dataThu[dataThu$ageWhenFirstCaptur == "Juv" & dataThu$sexe == "M", 22:27], 1, function(x) unique(x[!is.na(x)])), length)
(ncolsampled = sum(sapply(apply(dataThu[dataThu$ageWhenFirstCaptur == "Juv" & dataThu$sexe == "M", 22:27], 1, function(x) unique(x[!is.na(x)])), length) > 1))
# 16 juveniles have been sampled in more than one colony
# 8 males

(nataldisprate = ncolsampled/ntimesampled)
# Natal dispersal rate in Thuringia is even lower = 0.0337
# For males, natal dispersal rate = 0.055


#===============================#
#            END
#===============================#