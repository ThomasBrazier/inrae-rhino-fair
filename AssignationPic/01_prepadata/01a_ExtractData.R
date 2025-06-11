#===============================#
#            Assignation
#   Extraction des donnees de la BDD
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
source("../../Sources/packages.R")
source("fonctionsPrepaData.R")


#-------------------------------#
# Get the data
#-------------------------------#
# Read the data extracted from the Database
BDD=read.table("allInfosBDD.txt",h=T,sep="\t")

# Build the starting dataset
# Construit le jeu de donnees
# argument 1 : BDD contient toutes les donnees
# argument 2 : falseloci= max number of incorrect loci to keep the individual
# e.g. falseloci=2 means individuals with 0,1 ou 2 incomplete loci are kept
# argument 3 : withcolony=TRUE to get sampling data for each year
# argument 4 : sexmark=TRUE to keep sexual maerker dby8
unique=UniqueGenotypes(BDD,falseloci=2,withcolony=T,sexmark=F)
# le data set pret est dans uniqueGenotypesWithInfo.txt

#-------------------------------#
# Sample sizes
#-------------------------------#
effectifs=effectifsTable(unique)
View(effectifs)



