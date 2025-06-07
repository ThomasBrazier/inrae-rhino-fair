#===============================#
#            Assignation
#   Extraction des donnees de la BDD
#===============================#

# DIRECTORIES
# Working directories
# New users need to configure correct pathways
if (Sys.info()["sysname"]=="Darwin") {
  wd="/Volumes/Samsung_T5/INRA RHINO/R/AssignationThu/01_prepadata/"
  graphdir="/Volumes/Samsung_T5/INRA RHINO/R/Graph/"
  Rdir="/Volumes/Samsung_T5/INRA RHINO/R/"
  setwd(wd)
  source(paste(wd,"fonctionsPrepaData.R",sep=""))
  source(paste(wd,"fonctionsDataThu.R",sep=""))
}else {
  if (Sys.info()["sysname"]=="Windows"){
    wd="E:/INRA RHINO/R/AssignationThu/01_prepadata/"
    graphdir="E:/INRA RHINO/R/Graph/"
    Rdir="E:/INRA RHINO/R/"
    setwd(wd)
    source(paste(wd,"fonctionsPrepaData.R",sep=""))
    source(paste(wd,"fonctionsDataThu.R",sep=""))
  }
}


#-------------------------------#
# Lancer la fonction
#-------------------------------#
# Convertit le fichier All_Genotypes_FINAL_man.csv en uniqueGenotypesWithInfos.txt
unique=convertToUnique(falseloci=2,sexmark=FALSE)

#-------------------------------#
# Construction du tableau des effectifs
#-------------------------------#
effectifs=effectifsTable(unique)
View(effectifs)

#-------------------------------#
# Analyse du jeu de donnees uniqueGenotypesWithInfo
#-------------------------------#
totalF=sum(unique$sexe=="F") # nombre total de femelles
totalF
totalM=sum(unique$sexe=="M") # nombre total de males
totalM

# n total
totalF+totalM

# n juveniles
sum(unique$ageWhenFirstCaptur=="Juv")

# sexe ratio (% de males)
totalM/(totalF+totalM)*100

#-------------------------------#
# Routines de test
#-------------------------------
#source(paste(wd,"01fonction_tests.R",sep=""))
# Renseigner l'objet "unique", le nombre de tirages et si vous voulez un rapport plus detaille
tests(unique,10,verbose=F)



