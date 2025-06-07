#===============================#
#            MasterBayes
#  Full probability approach for pedigree inference
#         and parameters estimates
#===============================#

# clear global environment: remove all variables
rm(list=ls(all=TRUE))
# library(rstudioapi)
# Get the directory of the file & set working directory
# filedir=dirname(rstudioapi::getSourceEditorContext()$path)
# setwd(filedir)

#----------------------------------------------------------#
# Loading packages
# Check if packages are installed, install if necessary
source("../Sources/packages.R")


#----------------------------------#
# IMPORTANT
# MasterBayes is no longer available for install with CRAN
# and manual install of older versions can fail
# Hence the difficutly to re-run the analysis from the beginning
#
# All the commands, formulas and hyperparameters necessary top re-run MasterBayes are commmented and tabulated.
# Moreover, the models inferred with MasterBayes are saved in `MasterBayes\` and can be loaded and interpreted in `MasterBayes\MasterBayes.R`.


# Uncomment the MasterBayes commands if you want to re-run from scratch
# but you will have to find a solution to install MasterBayes
#
#----------------------------------#
      # install_version("MasterBayes",
      #                 version = "2.58",
      #                 upgrade = "never")
      # library(MasterBayes)



#-------------------------------#
# Aide MasterBayes
#-------------------------------
# vignette("Tutorial","MasterBayes")

#-------------------------------#
# MCMCPed function
#-------------------------------
# According to MasterBayes manual :
# Markov chain Monte Carlo methods for estimating the joint posterior distribution of a pedigree and
# the parameters that predict its structure using genetic and non-genetic data. These parameters can
# be associated with covariates of fecundity such as a sexually selected trait or age, or can be associated
# with spatial or heritable traits that relate parents to specific offspring. Population size, allele
# frequencies, allelic dropout rates, and stochastic genotyping error rates can also be simultaneously
# estimated.

##########################
# BEGIN : PARAMETER ESTIMATES FOR PICARDY
info=read.table("Picardy/uniqueGenotypesWithInfo79.txt",header=TRUE,stringsAsFactors = FALSE)
coord=read.table("Picardy/coordpic.txt",header=TRUE)
# IMPROVING DATA SET WITH SUPPLEMENTARY INFOS
# Add coordinates
info$long=rep(NA,nrow(info))
info$lat=rep(NA,nrow(info))
for (i in 1:nrow(info)) {
  info[i,30:31]=unlist(coord[which(as.character(coord[,1])==as.character(info[i,18])),2:3])
}
# CONVERSION DES COORDONNEES AU SYSTEME CARTESIEN (km)
cartesian=as.data.frame(info[,30:31])
colnames(cartesian)=c("long","lat")
library(sp)

install_version("rgdal", version = "1.6-7", update = "never")
library(rgdal)

coordinates(cartesian) <- c("long", "lat") # Creates a SpatialPoints object where x
# and y are the coordinates, assuming that "x" is the horizontal coordinate
# and "y" is the vertical one
proj4string(cartesian) <- CRS("+proj=longlat + ellps=WGS84") # This sets the
# projection of the data, assuming that you are using the WGS84 ellipsoid
# (if you do not know for sure, it is a safe assumption)
cartesian=spTransform(cartesian,CRS="+proj=utm +zone=27 +datum=WGS84") # Use
# whatever UTM zone you need.
cartesian=as.data.frame(coordinates(cartesian)) # This converts the SpatialPoints
info[,30:31]=cartesian/1000 # coordinates are integrated in the MasterBayes package in cartesian (km)

# Logical for offspring
info$offspring=rep(NA,nrow(info))
for (i in 1:nrow(info)) {
  if(info$ageWhenFirstCaptur[i]=="Juv"){
    info$offspring[i]=1
  }else{
    info$offspring[i]=0
  }
}
# Indicating cohort for multigenerational analyses
# only for offspring
# parents are -1 (birth year is unknwon)
info$cohort=rep(NA,nrow(info))
for (i in 1:nrow(info)) {
  if(info$offspring[i]==1){
    info$cohort[i]=(info$yearWhenFirstCaptur[i]-min(info$yearWhenFirstCaptur))+1
  }else{
    info$cohort[i]=0
  }
}
# Sexe information at the correct format
for(i in 1:nrow(info)){
  if(info[i,19]=='F'){
    info[i,19]=as.character("Female")
  }else{
    info[i,19]=as.character("Male")
  }
}

##### BEGINNING OF THE MODEL
# gc(reset=TRUE) # freeing memory allocation : MCMCPed need a lot of memory
# #####
# # Year 2013 ONLY (memory limitation)
# # Selection of individuals captured in 2013
# Gen=info[which(info$yearWhenFirstCaptur==2013),1:17]
# GdP=GdataPed(Gen[,2:17],Gen[,1],marker.type="MSW")
# # paternity is to be modelled as a function of distance
# # between offspring and male territories
# var1=expression(varPed(c("lat","long"), gender="Male",
#                         relational="OFFSPRING"))
# # individuals from the offspring generation are excluded as parents
# res1=expression(varPed("offspring", restrict=0))
# # individuals from the offspring generation or before (cohort equal or inferior) are excluded as parents
# # res1=expression(varPed("cohort", restrict=">"))
# 
# # mothers not from the offspring territory are excluded
# # offspring colonies provide information for the mother estimates
# res2=expression(varPed("col", gender="Female", relational="OFFSPRING",
#                         restrict="=="))
# 
# # ALLELIC FREQUENCIES
# # can be incorporated in the model
# 
# # START PARAMETERS
# # E1 : dropout rates
# # E2 : misstyping errors
# locus=read.table("Picardy/Input Locus 2.txt",header=TRUE)
# startParam=startPed(E1=0.05, estE1=TRUE, E2=0.05, estE2=TRUE,USsire=2000,estUSsire=TRUE)
# # startParam=startPed(G=Gen[,2:17], id=Gen[,1], E1=locus[2,], estE1=TRUE, E2=locus[3,], estE2=TRUE,estG=FALSE,
# #                     estUSsire=TRUE,USsire=2000)
# 
# Phen=info[which(info$yearWhenFirstCaptur==2013),c(1,18:19,21,30:33)]
# colnames(Phen)=c("id","col","sex","year","long","lat","offspring","cohort")
# # USsire = TRUE for unsampled fathers
# PdP=PdataPed(formula=list(var1,res1,res2), data=Phen,USsire=TRUE)
# tP=tunePed(beta=30)
# 
# model1=MCMCped(PdP=PdP, GdP=GdP, sP=startParam, tP=tP, nitt=1010000, thin=100, burnin=10000, jointP=FALSE)
# # Same model with joint distribution of pedigree (needs a lot of computational time)
# # model1joint=MCMCped(PdP=PdP, GdP=GdP, sP=startParam, tP=tP, nitt=110000, thin=1, burnin=10000, jointP=TRUE)
# gc(reset=TRUE) # freeing memory allocation : MCMCPed need a lot of memory
# save(model1,file="Picardy/model1")


##### ONCE IT IS WORKING FOR 2013 INDIVIDUALS
# IMPROVING WITH INDIVIDUALS OF 2013-2016 SAMPLINGS
# use of timevar according to Hadfield
# timevar is a possibitlity offered by MasterBayes to treat with longitudinal data and multigenerational pedigrees

# Females has only one offspring per year
# 
# gc(reset=TRUE) # freeing memory allocation : MCMCPed need a lot of memory
# #####
    # # All individuals
    # Gen=info[,1:17]
    # GdP=GdataPed(Gen[,2:17],Gen[,1],marker.type="MSW")
    # # paternity is to be modelled as a function of distance
    # # between offspring and male territories
    # var1=expression(varPed(c("lat","long"),lag=c(0,0), gender="Male",
    #                        relational="OFFSPRING"))
    # # individuals from the offspring generation are excluded as parents
    # res1=expression(varPed("offspring", restrict=0))
    # # individuals from the offspring generation or before (cohort equal or inferior) are excluded as parents
    # res2=expression(varPed("year",relational="OFFSPRING",restrict=">="))
    # # mothers not from the offspring territory are excluded
    # # offspring colonies provide information for the mother estimates
    # res3=expression(varPed("col", gender="Female", relational="OFFSPRING",restrict="=="))
    # # ALLELIC FREQUENCIES
    # # can be incorporated in the model
    # #☺ or calculated by MasterByes with the function extractA (default)
    # # locus=read.table("Picardy/Input Locus.txt",header=TRUE)
    # # START PARAMETERS
    # # E1 : dropout rates
    # # E2 : misstyping errors
    # startParam=startPed(E1=0.05, estE1=TRUE, E2=0.05, estE2=TRUE,USsire=2000,estUSsire=TRUE)
    # Phen=info[,c(1,18:19,21,30:33)]
    # colnames(Phen)=c("id","col","sex","year","long","lat","offspring","cohort")
    # # USsire = TRUE for unsampled fathers
    # PdP=PdataPed(formula=list(var1,res1,res3), data=Phen,timevar=Phen$year,USsire=TRUE,USdam=FALSE)
    # tP=tunePed(beta=30)
    # # Number of iterations 1010000
    # # Thinning of 100 for reducing autocorrelation
    # # burnin of 10000
    # model1=MCMCped(PdP=PdP, GdP=GdP, sP=startParam, tP=tP, nitt=1020000, thin=100, burnin=20000, jointP=FALSE)
    # gc(reset=TRUE)
    # save(model1,file="Picardy/model1")




##### PARAMETERS ESTIMATES
load(file="Picardy/model1")
summary(model1$beta)

# The mean shape of the negative exponential kernel function
mean(abs(model1$beta))
# The credible interval 95%
quantile(abs(model1$beta), c(0.025, 0.975))

# The mean dispersal distance is the inverse of Beta
1 / mean(abs(model1$beta))
# The credible interval 95%
1 / quantile(abs(model1$beta), c(0.025, 0.975))


output_pic = data.frame(population = "Picardy",
                   beta_mean = mean(abs(model1$beta)),
                   beta_lower_95 = quantile(abs(model1$beta), c(0.025)),
                   beta_upper_95 = quantile(abs(model1$beta), c(0.975)),
                   distance_mean = 1 / mean(abs(model1$beta)),
                   distance_lower_95 = 1 / quantile(abs(model1$beta), c(0.975)),
                   distance_upper_95 = 1 / quantile(abs(model1$beta), c(0.025)))

# Genotyping errors
summary(model1$E1)
summary(model1$E2)
# Unsampled fathers
summary(model1$USsire)
# posterior distribution of pedigree/ marginal distribution of parents
# one matrix per offspring : parents in lines and columns
# in each cell, the number of iterations where parents are find as parents of the offspring
summary(model1$P)

##### DIAGNOSTICS FOR CONVERGENCE
# see Hamra et al. 2013
# Trace + autocorr + density plots
plot(model1$beta)
# Genotyping errors
plot(model1$E1)
plot(model1$E2)
# Unsampled fathers
plot(model1$USsire)
# Autocorrelation of the estimated parameter
autocorr(model1$beta)
autocorr(model1$USsire)

# Gelman-Rubin diagnostic for 3 chains
# Hamra et al. 2013
# "assess model convergence. This diagnostic involves running
# multiple MCMC procedures, specifying the same
# model and prior information, from different starting
# values and comparing the variance within each chain
# with the variance between chains.
# Lack of model convergence
# is indicated when the variance between chains
# is larger than the variance within chains."

    # startParam=startPed(beta=-1,E1=0.05, estE1=TRUE, E2=0.05, estE2=TRUE,USsire=2000,estUSsire=TRUE)
    # model1b=MCMCped(PdP=PdP, GdP=GdP, sP=startParam, tP=tP, nitt=1020000, thin=100, burnin=20000,
    #                 jointP=FALSE,verbose=FALSE)
    # save(model1b,file="Picardy/model1b")
    # startParam=startPed(beta=-2,E1=0.05, estE1=TRUE, E2=0.05, estE2=TRUE,USsire=2000,estUSsire=TRUE)
    # model1c=MCMCped(PdP=PdP, GdP=GdP, sP=startParam, tP=tP, nitt=1020000, thin=100, burnin=20000,
    #                 jointP=FALSE,verbose=FALSE)
    # save(model1c,file="Picardy/model1c")
    # startParam=startPed(beta=0,E1=0.05, estE1=TRUE, E2=0.05, estE2=TRUE,USsire=2000,estUSsire=TRUE)
    # model1d=MCMCped(PdP=PdP, GdP=GdP, sP=startParam, tP=tP, nitt=1020000, thin=100, burnin=20000,
    #                 jointP=FALSE,verbose=FALSE)
    # save(model1d,file="Picardy/model1d")
    # startParam=startPed(beta=10,E1=0.05, estE1=TRUE, E2=0.05, estE2=TRUE,USsire=2000,estUSsire=TRUE)
    # model1e=MCMCped(PdP=PdP, GdP=GdP, sP=startParam, tP=tP, nitt=1020000, thin=100, burnin=20000,
    #                 jointP=FALSE,verbose=FALSE)
    # save(model1e,file="Picardy/model1e")
    # startParam=startPed(beta=10,E1=0.05, estE1=TRUE, E2=0.05, estE2=TRUE,USsire=2000,estUSsire=TRUE)
    # model1f=MCMCped(PdP=PdP, GdP=GdP, sP=startParam, tP=tP, nitt=1020000, thin=100, burnin=20000,
    #                 jointP=FALSE,verbose=FALSE)
    # save(model1f,file="Picardy/model1f")

load(file="Picardy/model1b")
load(file="Picardy/model1c")
load(file="Picardy/model1d")
load(file="Picardy/model1e")
load(file="Picardy/model1f")

# METRIC:
# According to Gelman & Rubin 1996
# "'potential scale reduction
#factor’ (PSRF), labelled sqrt(R), which is essentially the square root of the variances
# of the values of the scalar summary for all the simulated sequences mixed together,
# divided by the average of the variances within the separate sequences."


# So calculation of PSRF for the beta chain of estimates
PSRFbeta = sqrt(sd(c(model1b$beta,model1c$beta,model1d$beta,model1e$beta,model1f$beta))/
  mean(c(sd(model1b$beta),sd(model1c$beta),sd(model1d$beta),sd(model1e$beta),sd(model1f$beta))))
PSRFbeta

# And PSRF for USsire
PSRFUSsire = sqrt(sd(c(model1b$USsire,model1c$USsire,model1d$USsire,model1e$USsire,model1f$USsire))/
                  mean(c(sd(model1b$USsire),sd(model1c$USsire),sd(model1d$USsire),sd(model1e$USsire),sd(model1f$USsire))))
PSRFUSsire

# Gelman-Rubin index, variance between chains/variance within chains
# If index is close to 1, then convergence is assessed



##################################################################
##################################################################
##################################################################
#               THURINGIA
##################################################################
##########################
# BEGIN : PARAMETER ESTIMATES FOR GERMANY
info=read.table("Thuringia/uniqueGenotypesWithInfo79.txt",header=TRUE,stringsAsFactors = FALSE)
coord=read.table("Thuringia/coordThu.txt",header=TRUE)
# IMPROVING DATA SET WITH SUPPLEMENTARY INFOS
# Add coordinates
info$long=rep(NA,nrow(info))
info$lat=rep(NA,nrow(info))
for (i in 1:nrow(info)) {
  info[i,28:29]=unlist(coord[which(as.character(coord[,1])==as.character(info[i,18])),2:3])
}
# CONVERSION DES COORDONNEES AU SYSTEME CARTESIEN (km)
cartesian=as.data.frame(info[,28:29])
colnames(cartesian)=c("long","lat")
library(sp)
library(rgdal)
coordinates(cartesian) <- c("long", "lat") # Creates a SpatialPoints object where x
# and y are the coordinates, assuming that "x" is the horizontal coordinate
# and "y" is the vertical one
proj4string(cartesian) <- CRS("+proj=longlat + ellps=WGS84") # This sets the
# projection of the data, assuming that you are using the WGS84 ellipsoid
# (if you do not know for sure, it is a safe assumption)
cartesian=spTransform(cartesian,CRS="+proj=utm +zone=27 +datum=WGS84") # Use
# whatever UTM zone you need.
cartesian=as.data.frame(coordinates(cartesian)) # This converts the SpatialPoints
info[,28:29]=cartesian/1000 # coordinates are integrated in the MasterBayes package in cartesian (km)

# Logical for offspring
info$offspring=rep(NA,nrow(info))
for (i in 1:nrow(info)) {
  if(info$ageWhenFirstCaptur[i]=="Juv"){
    info$offspring[i]=1
  }else{
    info$offspring[i]=0
  }
}
# Indicating cohort for multigenerational analyses
# only for offspring
# parents are -1 (birth year is unknwon)
info$cohort=rep(NA,nrow(info))
for (i in 1:nrow(info)) {
  if(info$offspring[i]==1){
    info$cohort[i]=(info$yearWhenFirstCaptur[i]-min(info$yearWhenFirstCaptur))+1
  }else{
    info$cohort[i]=0
  }
}
# Sexe information at the correct format
for(i in 1:nrow(info)){
  if(info[i,19]=='F'){
    info[i,19]=as.character("Female")
  }else{
    info[i,19]=as.character("Male")
  }
}

##### BEGINNING OF THE MODEL
gc(reset=TRUE) # freeing memory allocation : MCMCPed need a lot of memory
    # #####
    # # All individuals of 2015 and 2016
    # # 2017 is not taken in this model because it seems to be an impossible configuration
    # # for MasterBayes
    # # "Error : no possible dams."
    # # Despite the loss of information, the model seemed conssitent with the results in Picardy
    # Gen=info[which(info$yearWhenFirstCaptur==2015 | info$yearWhenFirstCaptur==2016),1:17]
    # GdP=GdataPed(Gen[,2:17],Gen[,1],marker.type="MSW")
    # # paternity is to be modelled as a function of distance
    # # between offspring and male territories
    # var1=expression(varPed(c("lat","long"),lag=c(0,0), gender="Male",
    #                        relational="OFFSPRING"))
    # # individuals from the offspring generation are excluded as parents
    # res1=expression(varPed("offspring",lag=c(0,0), restrict=0))
    # # individuals from the offspring generation or before (cohort equal or inferior) are excluded as parents
    # # res2=expression(varPed("year",relational="OFFSPRING",restrict=">="))
    # # mothers not from the offspring territory are excluded
    # # offspring colonies provide information for the mother estimates
    # res3=expression(varPed("col", gender="Female",lag=c(0,0), relational="OFFSPRING",restrict="=="))
    # # ALLELIC FREQUENCIES
    # # can be incorporated in the model
    # # or calculated by MasterByes with the function extractA (default)
    # # locus=read.table("Picardy/Input Locus.txt",header=TRUE)
    # # START PARAMETERS
    # # E1 : dropout rates
    # # E2 : misstyping errors
    # startParam=startPed(E1=0.05, estE1=TRUE, E2=0.05, estE2=TRUE,USsire=200,estUSsire=TRUE)
    # Phen=info[which(info$yearWhenFirstCaptur==2015 | info$yearWhenFirstCaptur==2016),c(1,18:19,21,28:31)]
    # colnames(Phen)=c("id","col","sex","year","long","lat","offspring","cohort")
    # # USsire = TRUE for unsampled fathers
    # # res2 and res3 are not taken in account
    # # because of an incompatible configuration
    # # It is resulting in a loss of information, but results of this incomplete model seemed consistent
    # # with the complete model in Picardy
    # PdP=PdataPed(formula=list(var1,res1), data=Phen,timevar=Phen$year,USsire=TRUE,USdam=FALSE)
    # # Number of iterations 1010000
    # # Thinning of 100 for reducing autocorrelation
    # # burnin of 10000
    # model1=MCMCped(PdP=PdP, GdP=GdP, sP=startParam,nitt=1020000, thin=100, burnin=20000,
    #                jointP=FALSE)
    # gc(reset=TRUE)
    # save(model1,file="Thuringia/model1")

##### PARAMETERS ESTIMATES
load(file="Thuringia/model1")
summary(model1$beta)

# The mean shape of the negative exponential kernel function
mean(abs(model1$beta))
# The credible interval 95%
quantile(abs(model1$beta), c(0.025, 0.975))

# The mean dispersal distance is the inverse of Beta
1 / mean(abs(model1$beta))
# The credible interval 95%
1 / quantile(abs(model1$beta), c(0.025, 0.975))

output_thu = data.frame(population = "Thuringia",
                        beta_mean = mean(abs(model1$beta)),
                        beta_lower_95 = quantile(abs(model1$beta), c(0.025)),
                        beta_upper_95 = quantile(abs(model1$beta), c(0.975)),
                        distance_mean = 1 / mean(abs(model1$beta)),
                        distance_lower_95 = 1 / quantile(abs(model1$beta), c(0.975)),
                        distance_upper_95 = 1 / quantile(abs(model1$beta), c(0.025)))

output = rbind(output_pic, output_thu)

write.table(output, "../Tables/table_masterbayes.tsv", sep = "\t",
            col.names = T, row.names = F)


# Genotyping errors
summary(model1$E1)
summary(model1$E2)
# Unsampled fathers
summary(model1$USsire)
# posterior distribution of pedigree/ marginal distribution of parents
# one matrix per offspring : parents in lines and columns
# in each cell, the number of iterations where parents are find as parents of the offspring
summary(model1$P)

##### DIAGNOSTICS FOR CONVERGENCE
# see Hamra et al. 2013
# Trace + autocorr + density plots
plot(model1$beta)


# Genotyping errors
plot(model1$E1)
plot(model1$E2)
# Unsampled fathers
plot(model1$USsire)
# Autocorrelation of the estimated parameter
autocorr(model1$beta)
autocorr(model1$USsire)


# Gelman-Rubin diagnostic for 3 chains
# Hamra et al. 2013
# "assess model convergence. This diagnostic involves running
# multiple MCMC procedures, specifying the same
# model and prior information, from different starting
# values and comparing the variance within each chain
# with the variance between chains.
# Lack of model convergence
# is indicated when the variance between chains
# is larger than the variance within chains."

    # startParam=startPed(beta=-1,E1=0.05, estE1=TRUE, E2=0.05, estE2=TRUE,USsire=200,estUSsire=TRUE)
    # model1b=MCMCped(PdP=PdP, GdP=GdP, sP=startParam, tP=tP, nitt=1020000, thin=100, burnin=20000,
    #                 jointP=FALSE,verbose=FALSE)
    # save(model1b,file="Thuringia/model1b")
    # startParam=startPed(beta=-2,E1=0.05, estE1=TRUE, E2=0.05, estE2=TRUE,USsire=200,estUSsire=TRUE)
    # model1c=MCMCped(PdP=PdP, GdP=GdP, sP=startParam, tP=tP, nitt=1020000, thin=100, burnin=20000,
    #                 jointP=FALSE,verbose=FALSE)
    # save(model1c,file="Thuringia/model1c")
    # startParam=startPed(beta=0,E1=0.05, estE1=TRUE, E2=0.05, estE2=TRUE,USsire=200,estUSsire=TRUE)
    # model1d=MCMCped(PdP=PdP, GdP=GdP, sP=startParam, tP=tP, nitt=1020000, thin=100, burnin=20000,
    #                 jointP=FALSE,verbose=FALSE)
    # save(model1d,file="Thuringia/model1d")
    # startParam=startPed(beta=10,E1=0.05, estE1=TRUE, E2=0.05, estE2=TRUE,USsire=200,estUSsire=TRUE)
    # model1e=MCMCped(PdP=PdP, GdP=GdP, sP=startParam, tP=tP, nitt=1020000, thin=100, burnin=20000,
    #                 jointP=FALSE,verbose=FALSE)
    # save(model1e,file="Thuringia/model1e")
    # startParam=startPed(beta=10,E1=0.05, estE1=TRUE, E2=0.05, estE2=TRUE,USsire=200,estUSsire=TRUE)
    # model1f=MCMCped(PdP=PdP, GdP=GdP, sP=startParam, tP=tP, nitt=1020000, thin=100, burnin=20000,
    #                 jointP=FALSE,verbose=FALSE)
    # save(model1f,file="Thuringia/model1f")

load(file="Thuringia/model1b")
load(file="Thuringia/model1c")
load(file="Thuringia/model1d")
load(file="Thuringia/model1e")
load(file="Thuringia/model1f")

# METRIC:
# According to Gelman & Rubin 1996
# "'potential scale reduction
#factor’ (PSRF), labelled sqrt(R), which is essentially the square root of the variances
# of the values of the scalar summary for all the simulated sequences mixed together,
# divided by the average of the variances within the separate sequences."


# So calculation of PSRF for the beta chain of estimates
PSRFbeta = sqrt(sd(c(model1b$beta,model1c$beta,model1d$beta,model1e$beta,model1f$beta))/
                  mean(c(sd(model1b$beta),sd(model1c$beta),sd(model1d$beta),sd(model1e$beta),sd(model1f$beta))))
PSRFbeta

# And PSRF for USsire
PSRFUSsire = sqrt(sd(c(model1b$USsire,model1c$USsire,model1d$USsire,model1e$USsire,model1f$USsire))/
                    mean(c(sd(model1b$USsire),sd(model1c$USsire),sd(model1d$USsire),sd(model1e$USsire),sd(model1f$USsire))))
PSRFUSsire

# Gelman-Rubin index, variance between chains/variance within chains
# If index is close to 1, then convergence is assessed


