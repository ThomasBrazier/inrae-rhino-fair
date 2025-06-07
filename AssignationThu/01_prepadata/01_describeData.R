#===============================#
#            Assignation
#   Description des donnees de la BDD
#===============================#

# DIRECTORIES
# Working directories
# New users need to configure correct pathways
if (Sys.info()["sysname"]=="Darwin") {
  wd="/Volumes/Samsung_T5/INRA RHINO/R/AssignationThu/01_prepadata/"
  graphdir="/Volumes/Samsung_T5/INRA RHINO/R/Graph/"
  Rdir="/Volumes/Samsung_T5/INRA RHINO/R/"
  setwd(wd)
}else {
  if (Sys.info()["sysname"]=="Windows"){
    wd="E:/INRA RHINO/R/AssignationThu/01_prepadata/"
    graphdir="E:/INRA RHINO/R/Graph/"
    Rdir="E:/INRA RHINO/R/"
    setwd(wd)
  }
}



#-------------------------------#
# Nombre d'alleles par marqueur
#-------------------------------
genotypes=read.table("uniqueGenotypesWithInfo99.txt",h=T)

for (i in seq(2,18,2)) {
  cat(colnames(genotypes)[i]," : ",length(unique(as.vector(t(genotypes[,i:(i+1)]))))-1,"\n") # nombre d'elements - 1 a cause des zeros
}


#-------------------------------#
# Description des effectifs
#-------------------------------
effectifs=read.table("effectifs99.txt",h=T)


# total number of juveniles
juv=length(genotypes$idind[which(genotypes$ageWhenFirstCaptur=="Juv")])
juv
binom.test(c(sum(genotypes$ageWhenFirstCaptur=="Juv"),sum(genotypes$ageWhenFirstCaptur=="Adult")))

# total number of candidate fathers
fa=length(genotypes$idind[which(genotypes$ageWhenFirstCaptur=="Adult" & genotypes$sexe=="M")])
fa
# total number of candidate mothers
mo=length(genotypes$idind[which(genotypes$ageWhenFirstCaptur=="Adult" & genotypes$sexe=="F")])
mo

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
col=paste("Thu",seq(20,50,1),sep="")
colsize=c()

for (c in col) {
  colsize=c(colsize,length(which(genotypes$idcol==c)))
}
colsize=as.data.frame(colsize)
colsize=t(colsize)
colnames(colsize)=col
colsize
min(colsize)
max(colsize)

# BOOTSTRAP MEAN COLONY SIZE
boots=numeric(10000)
for (i in 1:10000){
  boots[i]=mean(sample(colsize,replace=T))
}
limBoot1=quantile(boots,c(.025,.975))
limBoot1
mboot=mean(boots)
mboot
etboot=sqrt(var(boots))
etboot

#-------------------------------#
# Proportion of males
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
