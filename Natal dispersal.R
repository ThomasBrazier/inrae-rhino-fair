#==========================================================#
#            COLONY
#     Natal dispersal
#     How juveniles genotyped dispersed in subsequent years
#==========================================================#

#==========================================================#
# LOADING ENVIRONMENT
#==========================================================#

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

#==========================================================#
# Import data
#==========================================================#
dataPic = read.table("Data/Pic/uniqueGenotypesWithInfo.txt", header = TRUE)
dataThu = read.table("Data/Thu/uniqueGenotypesWithInfo.txt", header = TRUE)
# Migrant fathers
MigrantFathersPic = read.table("Tables/MigrantFathersPic.txt", header = TRUE)
MigrantFathersThu = read.table("Tables/MigrantFathersThu.txt", header = TRUE)


#==========================================================#
# How juveniles genotyped dispersed in subsequent years
#==========================================================#


#------------------------------------------------#
# Picardy
#------------------------------------------------#
# How many male juveniles were found in the same colony subsequently or did they dispersed?
# Individuals
indiv = as.character(dataPic$idind[which(dataPic$ageWhenFirstCaptur == "Juv" & dataPic$sexe == "M")])
# In which colonies were sampled each male juvenile
sampledcol = apply(dataPic[which(dataPic$ageWhenFirstCaptur == "Juv" & dataPic$sexe == "M"), 22:29], 1, function(x) unique(x, na.rm =TRUE))
# To which colony were they assigned
assignedcol = as.character(dataPic$idcol[which(dataPic$ageWhenFirstCaptur == "Juv" & dataPic$sexe == "M")])
# For each individual, is there a colony different than the assigned one? Finally, is there more than one colony?
# I should filter out NAs, but simply said, how many individuals have more than 2 items in the sampledcol list
length(indiv[lapply(sampledcol, length) > 2])
# 2 individuals
length(indiv[lapply(sampledcol, length) > 2])/length(indiv)
# 0.3% of all male juvenile sampled where found only in one colony

# We can do better, how many of male juveniles sampled more than once where found in a single colony?
# Subset only males sampled multiple times
sampledcolall = dataPic[which(dataPic$ageWhenFirstCaptur == "Juv" & dataPic$sexe == "M"), 22:29]
subsample = sampledcolall[which(apply(sampledcolall, 1, function(x) length(x) - sum(is.na(x))) > 1),]
subsample = apply(subsample, 1, function(x) unique(x, na.rm =TRUE))
# Only 56 juveniles sampled multiple times
sum(lapply(subsample, length) > 2)
# 2 juveniles sampled in two different colonies
sum(lapply(subsample, length) > 2)/length(subsample)
# 3.5% of juveniles sampled multiple times were in more than one colony


#------------------------------------------------#
# Thuringia
#------------------------------------------------#
# Subset only males sampled multiple times
sampledcolall = dataThu[which(dataThu$ageWhenFirstCaptur == "Juv" & dataThu$sexe == "M"), 22:27]
subsample = sampledcolall[which(apply(sampledcolall, 1, function(x) length(x) - sum(is.na(x))) > 1),]
subsample = apply(subsample, 1, function(x) unique(x, na.rm =TRUE))
# Only 146 juveniles sampled multiple times
sum(lapply(subsample, length) > 2)
# 8 juveniles sampled in two different colonies
sum(lapply(subsample, length) > 2)/length(subsample)
# 5.5% of juveniles sampled multiple times were in more than one colony




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



#------------------------------------------------#
# END
#------------------------------------------------#
