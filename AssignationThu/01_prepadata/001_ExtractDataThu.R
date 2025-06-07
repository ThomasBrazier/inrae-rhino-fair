#===============================#
#            Assignation
#   Extraction des donnees de la BDD
#===============================#

# DIRECTORIES

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
#write.table(effectifs,paste(wd,"effectifs.txt",sep=""))


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
# % de males
totalM/(totalF+totalM)*100

#### Adultes uniquement
totalFAd=sum(unique$sexe=="F" & unique$ageWhenFirstCaptur=="Adult") # nombre total de femelles adultes
totalFAd 
totalMAd=sum(unique$sexe=="M" & unique$ageWhenFirstCaptur=="Adult") # nombre total de males adultes
totalMAd
# n total
totalFAd+totalMAd
# % de males adultes
totalMAd/(totalFAd+totalMAd)*100

