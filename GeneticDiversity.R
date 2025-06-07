#===============================#
#            COLONY
#    Genetic diversity
#===============================#

#==========================================================
# LOADING ENVIRONMENT
#==========================================================

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


#---------------------------------
# Import data
#---------------------------------
dataPic = read.table("Data/Pic/uniqueGenotypesWithInfo.txt", header = TRUE)
dataThu = read.table("Data/Thu/uniqueGenotypesWithInfo.txt", header = TRUE)


dataPic = dataPic[,c(18, 1:17)]
dataThu = dataThu[,c(18, 1:17)]

locus = cbind(paste(dataPic[,3], dataPic[,4], sep = ":"),
                   paste(dataPic[,5], dataPic[,6], sep = ":"),
                   paste(dataPic[,7], dataPic[,8], sep = ":"),
                   paste(dataPic[,9], dataPic[,10], sep = ":"),
                   paste(dataPic[,11], dataPic[,12], sep = ":"),
                   paste(dataPic[,13], dataPic[,14], sep = ":"),
                   paste(dataPic[,15], dataPic[,16], sep = ":"),
                   paste(dataPic[,17], dataPic[,18], sep = ":"))
colnames(locus) = colnames(dataPic)[c(3, 5, 7, 9, 11, 13, 15, 17)]
ind = as.character(dataPic$idind)
population = as.character(dataPic$idcol)
Pic = df2genind(locus, ploidy = 2, ind.names = ind, pop = population, sep = ":")

nAll(Pic)
summary(Pic)
basic.stats(Pic)$overall

locus = cbind(paste(dataThu[,3], dataThu[,4], sep = ":"),
              paste(dataThu[,5], dataThu[,6], sep = ":"),
              paste(dataThu[,7], dataThu[,8], sep = ":"),
              paste(dataThu[,9], dataThu[,10], sep = ":"),
              paste(dataThu[,11], dataThu[,12], sep = ":"),
              paste(dataThu[,13], dataThu[,14], sep = ":"),
              paste(dataThu[,15], dataThu[,16], sep = ":"),
              paste(dataThu[,17], dataThu[,18], sep = ":"))
colnames(locus) = colnames(dataThu)[c(3, 5, 7, 9, 11, 13, 15, 17)]
ind = as.character(dataThu$idind)
population = as.character(dataThu$idcol)
Thu = df2genind(locus, ploidy = 2, ind.names = ind, pop = population, sep = ":")

nAll(Thu)
summary(Thu)
basic.stats(Thu)$overall

