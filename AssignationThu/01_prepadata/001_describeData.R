#===============================#
#            Assignation
#   Description des donnees de la BDD
#===============================#

# DIRECTORIES

if (Sys.info()["sysname"]=="Darwin") {
  wd="/Volumes/Samsung_T5/INRA RHINO/R/01_prepadata/"
  graphdir="/Volumes/Samsung_T5/INRA RHINO/R/Graph/"
  Rdir="/Volumes/Samsung_T5/INRA RHINO/R/"
  setwd(wd)
}else {
  if (Sys.info()["sysname"]=="Windows"){
    wd="E:/INRA RHINO/R/Assignation/01_prepadata/"
    graphdir="E:/INRA RHINO/R/Graph/"
    Rdir="E:/INRA RHINO/R/"
    setwd(wd)
  }
}



#-------------------------------#
# Nombre d'alleles par marqueur
#-------------------------------
genotypes=read.table("uniqueGenotypesWithInfo.txt",h=T)

for (i in seq(2,18,2)) {
  cat(colnames(genotypes)[i]," : ",length(unique(as.vector(t(genotypes[,i:(i+1)]))))-1,"\n") # nombre d'elements - 1 a cause des zeros
}


#-------------------------------#
# Description des effectifs
#-------------------------------
effectifs=read.table("effectifs.txt",h=T)


# total number of juveniles
juv=length(genotypes$idind[which(genotypes$ageWhenFirstCaptur=="Juv" & genotypes$yearWhenFirstCaptur==2013)])
binom.test(c(sum(genotypes$ageWhenFirstCaptur=="Juv"),sum(genotypes$ageWhenFirstCaptur=="Adult")))

# total number of candidate fathers
fa=length(genotypes$idind[which(genotypes$ageWhenFirstCaptur=="Adult" & genotypes$sexe=="M" & genotypes$yearWhenFirstCaptur==2013)])

# total number of candidate mothers
mo=length(genotypes$idind[which(genotypes$ageWhenFirstCaptur=="Adult" & genotypes$sexe=="F" & genotypes$yearWhenFirstCaptur==2013)])

# proportion of each
juv/(juv+fa+mo)
fa/(juv+fa+mo)
mo/(juv+fa+mo)

# nombre de genotypes incomplets
count=0
for (i in 1:nrow(genotypes)) {
  if (prod(genotypes[i,2:17])==0){
    count=count+1
  }
}
count

#-------------------------------#
# Description des colonies
#-------------------------------

# mean population size of colonies
col=c('CXSGT','M244','M399','M809','M811','M813','M815','M831','M1079','M1102','M1106',
      'M1112','M1113','M1551','M1975','M1979','PP1','PP3')
colsize=c()
# for (c in col) {
#   colsize=c(colsize,length(which(genotypes$idcol==c & genotypes$yearWhenFirstCaptur==2013)))
# }
for (c in col) {
  colsize=c(colsize,length(which(genotypes$idcol==c)))
}
colsize=as.data.frame(colsize)
colsize=t(colsize)
colnames(colsize)=col
colsize
# Effectif moyen d'une colonie
meanCol=mean(colsize)
meanCol
sqrt(var(as.vector(colsize)))
sdCol=sd(colsize)
sdCol
IC=meanCol+c(-1,1)*1.96*sdCol/sqrt(length(col))
IC

#-------------------------------#
# Porportion of males
#-------------------------------

# Proportion of males par colonie
colSex=c()
colF=c()
colTot=c()
for (c in col) {
  colSex=c(colSex,binom.test(c(sum(genotypes$sexe=="F" & genotypes$idcol==c & genotypes$ageWhenFirstCaptur=="Adult"),
                                  sum(genotypes$sexe=="M" & genotypes$idcol==c & genotypes$ageWhenFirstCaptur=="Adult")))$estimate)
  colF=c(colF,sum(genotypes$sexe=="F" & genotypes$idcol==c))
  colTot=c(colTot,sum(genotypes$idcol==c))
}
binom.test(colF[1],colTot[1])
# sex-ratio sur l'ensemble de la population
binom.test(c(sum(genotypes$sexe=="F" & genotypes$ageWhenFirstCaptur=="Adult"),sum(genotypes$sexe=="M" & genotypes$ageWhenFirstCaptur=="Adult")))

#-------------------------------#
# Sex ratio
#-------------------------------
# n males/n females
(colTot-colF)/colF
mean((colTot-colF)/colF)


#-------------------------------#
# Philopatry of females
#-------------------------------
# Chercher les femelles qui ont plus d'une colonie d'echantillonnage
# i.e. les femelles qui changent de colonie durant leur vie
# pour chaque femelle, on compte le nombre de colonies presentes dans les colonnes d'echantillonnage

# selection des ech. de femelles, juveniles et adultes ensemble
femalesSampled=genotypes[which(genotypes$sexe=="F" & genotypes$idcol!="P412"),24:31]
# compte les colonies ind. par ind.
for (i in 1:nrow(femalesSampled)) {
  colSamp=c(as.character(femalesSampled[i,1]),
            as.character(femalesSampled[i,2]),
            as.character(femalesSampled[i,3]),
            as.character(femalesSampled[i,4]),
            as.character(femalesSampled[i,5]),
            as.character(femalesSampled[i,6]),
            as.character(femalesSampled[i,7]),
            as.character(femalesSampled[i,8]))
  femalesSampled$nbCol[i]=length(unique(na.omit(colSamp)))
}
# Indices de femelles qui ont change de colonie au moins une fois
which(femalesSampled$nbCol>1)
length(which(femalesSampled$nbCol>1))
