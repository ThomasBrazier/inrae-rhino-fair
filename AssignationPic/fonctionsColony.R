#===============================#
#          COLONY
#         Results
#    Analyse des resultats
#===============================#


#---------------------------------#
# MAKING THE FINAL RESULT DATASET
# COMBINING THE TWO DATASETS
#---------------------------------#

# Combiner les deux jeux de donnees
# Argument : "99" prend le jeu de donnees 9/9 loci complets
#           "79" prend le jeu de donnees 7/9 loci complets
#           "all" prend les deux jeux de donnees
combineParentageMethods=function(dataset="all",sex="M"){
  Paternities=data.frame()

  if(sex=="M"){
    if(dataset=="all"){
      # enlever les paternites en doublon (1 seul echantillon garde)
      # Priorite aux resultats 99 sur 79 quand deux paternites differentes retrouvees
      # prend tous les resultats de paternity 99
      # rajoute les paternites 79 une par une si elle ne sont pas deja dans le data frame
      Paternities=paternity2
      for (i in 1:nrow(paternity)) {
        if (!(paternity[i,1] %in% Paternities)){
          Paternities=rbind(Paternities,paternity[i,])
        }
      }
    }else{
      if(dataset=="79"){
        Paternities=paternity
      }else{
        if(dataset=="99"){
          Paternities=paternity2
        }else{
          cat("Erreur sur le choix de la methode.")
        }
      }
    }
    write.table(Paternities,"paternity.txt",col.names=FALSE,row.names=FALSE)
  }else{
    if(sex=="F"){
      if(dataset=="all"){
        # enlever les paternites en doublon (1 seul echantillon garde)
        # Priorite aux resultats 99 sur 79 quand deux paternites differentes retrouvees
        # prend tous les resultats de paternity 99
        # rajoute les paternites 79 une par une si elle ne sont pas deja dans le data frame
        Paternities=maternity2
        for (i in 1:nrow(maternity)) {
          if (!(maternity[i,1] %in% Paternities)){
            Paternities=rbind(Paternities,maternity[i,])
          }
        }
      }else{
        if(dataset=="79"){
          Paternities=maternity
        }else{
          if(dataset=="99"){
            Paternities=maternity2
          }else{
            cat("Erreur sur le choix de la methode.")
          }
        }
      }
      write.table(Paternities,"maternity.txt",col.names=FALSE,row.names=FALSE)
    }
  }
  return(Paternities)
}

#---------------------------------#
# MAKING THE FINAL RESULT DATASET
# AND EXCLUDING FALSE POSITIVE
#---------------------------------#

# Recuperation de tous les juveniles assignes a un parent
# construction du fichier des juveniles
# Recuperation des paternites
# Construction du fichier des paternites
# Exclusion des fausses assignations

# 1/ path, chemin d'acces aux resultats de Colony
# 2/ seuil d'exclusion des fausses paternites sur le critere de proportion d'assignation du bestFather
constructResults=function(path="",excl=0.5){
  # compte le nombre de runs
  # combien de fichiers contenant .Paternity dans le dossier
  nrun=length(list.files(path=path,pattern="\\.Paternity"))
  
  # ecrit un fichier des juveniles
  # concatene tous les fichiers .Paternity et .Maternity de tous les runs
  cat("Creation du fichier des juveniles.\n")
  offsp=c()
  for (i in 1:nrun) {
    offspf=read.table(paste(path,list.files(path=path,pattern="\\.Maternity")[i],sep=""),h=T,fill=TRUE)
    offspm=read.table(paste(path,list.files(path=path,pattern="\\.Paternity")[i],sep=""),h=T,fill=TRUE)
    offsp=c(offsp,offspf$OffspringID,offspm$OffspringID)
  }
  # enleve les doublons
  offsp=unique(offsp)
  cat("Ecriture du fichier des juveniles.\n")
  write.table(offsp,"Offsprings.txt",col.names=F,row.names=F)

  # ecrit un fichier des paternites
  # Dresser une liste de couples offspring-father retrouves par Colony (parente validee selon critere a definir)
  # liste de tous les juveniles avec un pere, sans doublons
  listPat=c()
  for (i in 1:nrun) {
    offspm=read.table(paste(path,list.files(path=path,pattern="\\.Paternity")[i],sep=""),h=T,fill=TRUE)
    listPat=c(listPat,offspm$OffspringID)
  }
  # enleve les doublons
  listPat=unique(listPat)
  # chaque paternite est enregistree dans une matrice avec les juveniles en y et le numero du run en x
  # et en premiere colonne l'ID du juvenile
  cat("Creation du fichier des paternites.\n")
  paternities=matrix(ncol=(nrun+1),nrow=length(listPat))
  paternities[,1]=listPat
  for (i in 1:nrun) {
    pat=read.table(paste(path,list.files(path=path,pattern="\\.Paternity")[i],sep=""),h=T,fill=TRUE)
    for (j in 1:nrow(pat)) {
      paternities[which(paternities[,1]==pat$OffspringID[j]),i+1]=pat$InferredDad1[j]
    }
  }
  cat(length(listPat),"paternites a ce stade.\n")
  
  # Exclusion des paternites en dessous du seuil fixe
  # proportion d'assignation du bestFather
  # d'abord on identifie le bestFather, celui qui revient le plus souvent
  # ensuite on calcule sa proportion d'assignation
  cat("Exclusion des paternites en dessous du seuil.\n")
  # identifier le meilleur pere : celui qui revient le plus souvent dans les assignations
  # !!!! on utilise les numeros de lignes : ils doivent concorder entre bestFather et paternities
  bestFather=c()
  for (i in 1:nrow(paternities)){
    uniq=unique(na.omit(as.character(paternities[i,2:(nrun+1)])))
    
    if(length(uniq)==0){
      bestFather[i]=NA
    }else{
      best=c(NA,0)
      for (j in 1:length(uniq)) {
        occur=length(which(paternities[i,2:(nrun+1)]==uniq[j]))
        if(occur>best[2]){
          best[1]=uniq[j]
          best[2]=occur
        }
      }
      bestFather[i]=best[1]
    }
  }
  # suppression des cellules qui ne contiennent pas le meilleur pere
  for (i in 1:nrow(paternities)) {
    for (j in 2:(nrun+1)){
      if(!is.na(paternities[i,j]) & paternities[i,j]!=bestFather[i]){
        paternities[i,j]=NA
      }
    }
  }
  # suppression des lignes qui sont en dessous du seuil d'exclusion
  # fais la liste des lignes a conserver
  index=c()
  for (i in 1:nrow(paternities)) {
    if (sum(!is.na(paternities[i,2:(nrun+1)]))/nrun > excl){
      index=c(index,i)
    }
  }
  paternities=paternities[index,]
  liste=c()
  for (i in 1:nrow(paternities)) {
    liste[i]=unique(paternities[i,2:(nrun+1)])[!is.na(unique(paternities[i,2:(nrun+1)]))]
  }
  paternities=cbind(paternities[,1],liste)
  colnames(paternities)=c("offspring","father")
  
  cat(nrow(paternities),"paternites finalement retenues avec",length(unique(liste)),"peres.\n")
  
  cat("Ecriture du fichier des paternites.\n")
  write.table(paternities,"Paternity.txt",col.names=F,row.names=F)  
  # renvoie les paternites
  return(paternities)
}


#---------------------------------#
# DISTRIBUTION OF NUMBER OF RUNS
# FOR WHICH A FATHER IS FOUND
#---------------------------------#

# Chaque pere trouve pour un juvenile apparait combien de fois
# l'idee est de voir si on a bien une distribution bimodale comme attendu d'apres simulations
# quelques bons peres revenant souvent et beaucoup de mauvais peres revenant rarement

# Tableau
# 1/ Offspring
# 2/ father
# 3/ number of occurences of that particular dyad

# 1/ path, chemin d'acces aux resultats de Colony

distBestFather=function(path=""){
  # compte le nombre de runs
  # combien de fichiers contenant .Paternity dans le dossier
  nrun=length(list.files(path=path,pattern="\\.Paternity"))
  
  tab=data.frame(offspring=NA,father=NA,occurences=NA)
  tab=tab[-1,]
  # parcourt tous les fichiers Paternity
  # si le couple offspring-father existe deja dans tab
  #       ajoute 1 a la case occurences
  # sinon ajoute le couple et initialise la case occurences a 1
  for (i in 1:nrun) {
    run=read.table(paste(path,list.files(path=path,pattern="\\.Paternity")[i],sep=""),h=T,fill=TRUE)
    # pour chaque ligne du fichier ouvert
    for (l in 1:nrow(run)) {
      index=which(tab[,1]==run$OffspringID[l] & tab[,2]==run$InferredDad1[l])
      if (length(index)==0){
        tab=rbind(tab,c(run$OffspringID[l],run$InferredDad1[l],1))
      }else{
        tab[index,3]=tab[index,3]+1
      }
    }
  }
  # Colonie de l'offsp
  info=read.table("uniqueGenotypesWithInfo.txt",h=T)
  colOffsp=c()
  for(i in 1:nrow(tab)){
    colOffsp[i]=as.character(info$idcol[which(info$idind==tab[i,1])])
  }
  # EXCLUSION DE LA COLONIE M399
  # colonie du pere
  colFather=c()
  for(i in 1:nrow(tab)){
    colFather[i]=as.character(info$idcol[which(info$idind==tab[i,2])])
  }
  # enlever les individus (pere ou offspring) venant d'une des colonies exclues
  # dyads impossibles ,"M1979","M1975","CXSGT","M1079"
  excludedCol=c("M399")
  tab=tab[which(!(colOffsp %in% excludedCol)),]
  tab=tab[which(!(colFather %in% excludedCol)),]
  
  colnames(tab)=c("offspring","father","occurences")
  write.table(tab,"distributionBestFathers.txt",col.names=T,row.names=F)  
  # renvoie les paternites
  return(tab)
}


#---------------------------------#
# DESCRIPTION OF RESULTS
#---------------------------------
# complete les resultats avec :
# 1/ colonie offsp
# 2/ colonie father
# 3/ year first captur offsp
# 4/ year first captur father
# 5/ age first captur Father
# 6/ Sexe juvenile


describeResults=function(x=paternity){
  resultsWithInfo=data.frame(offspring=x[,1],father=x[,2],colonyOffsp=rep(NA,nrow(x)),colonyFather=rep(NA,nrow(x)),
                             yearOffsp=rep(NA,nrow(x)),yearFather=rep(NA,nrow(x)),ageFather=rep(NA,nrow(x)),
                             sexOffsp=rep(NA,nrow(x)),y2013_1=rep(NA,nrow(x)),y2013_2=rep(NA,nrow(x)),
                             y2014_1=rep(NA,nrow(x)),y2014_2=rep(NA,nrow(x)),
                             y2015_1=rep(NA,nrow(x)),y2015_2=rep(NA,nrow(x)),
                             y2016_1=rep(NA,nrow(x)),y2016_2=rep(NA,nrow(x)))
  # Colonie de l'offsp
  info=read.table("uniqueGenotypesWithInfo.txt",h=T)
  for(i in 1:nrow(resultsWithInfo)){
    resultsWithInfo[i,3]=as.character(info$idcol[which(info$idind==resultsWithInfo[i,1])])
  }
  # colonie du pere
  for(i in 1:nrow(resultsWithInfo)){
    resultsWithInfo$colonyFather[i]=as.character(info$idcol[which(info$idind==resultsWithInfo[i,2])])
  }
  # enlever les individus (pere ou offspring) venant d'une des colonies exclues
  # dyads impossibles ,"M1979","M1975","CXSGT","M1079"
  excludedCol=c("M399")
  resultsWithInfo=resultsWithInfo[which(!(resultsWithInfo$colonyOffsp %in% excludedCol)),]
  resultsWithInfo=resultsWithInfo[which(!(resultsWithInfo$colonyFather %in% excludedCol)),]
  # year first captur offsp
  for(i in 1:nrow(resultsWithInfo)){
    resultsWithInfo$yearOffsp[i]=as.character(info$yearWhenFirstCaptur[which(info$idind==resultsWithInfo[i,1])])
  }
  # year first captur Father
  for(i in 1:nrow(resultsWithInfo)){
    resultsWithInfo$yearFather[i]=as.character(info$yearWhenFirstCaptur[which(info$idind==resultsWithInfo[i,2])])
  }  
  # age premiere capture Father
  for(i in 1:nrow(resultsWithInfo)){
    resultsWithInfo$ageFather[i]=as.character(info$ageWhenFirstCaptur[which(info$idind==resultsWithInfo[i,2])])
  }
  # sexe juv
  for(i in 1:nrow(resultsWithInfo)){
    resultsWithInfo$sexOffsp[i]=as.character(info$sexe[which(info$idind==resultsWithInfo[i,1])])
  }
  
  # sessions echantillonnage
  for (s in 1:8) {
    for(i in 1:nrow(resultsWithInfo)){
      resultsWithInfo[i,(8+s)]=as.character(info[which(info$idind==resultsWithInfo[i,2]),(21+s)])
    }
  }
  
  resultsWithInfo=resultsWithInfo[order(resultsWithInfo$father),]

  write.table(resultsWithInfo,"resultsWithInfo.txt",col.names=T,row.names=F)  
  
  return(resultsWithInfo)
}




#---------------------------------#
# TREATING MOTHERS
#---------------------------------#

# Recuperation de tous les juveniles assignes a un parent
# construction du fichier des juveniles
# Recuperation des paternites
# Construction du fichier des paternites
# Exclusion des fausses assignations

# 1/ path, chemin d'acces aux resultats de Colony
# 2/ seuil d'exclusion des fausses paternites sur le critere de proportion d'assignation du bestFather
constructResultsMums=function(path="",excl=0.5){
  # compte le nombre de runs
  # combien de fichiers contenant .Maternity dans le dossier
  nrun=length(list.files(path=path,pattern="\\.Maternity"))
  

  # ecrit un fichier des maternites
  # Dresser une liste de couples offspring-mother retrouves par Colony (parente validee selon critere a definir)
  # liste de tous les juveniles avec une mere, sans doublons
  listMat=c()
  for (i in 1:nrun) {
    offspm=read.table(paste(path,list.files(path=path,pattern="\\.Maternity")[i],sep=""),h=T,fill=TRUE)
    listMat=c(listMat,offspm$OffspringID)
  }
  # enleve les doublons
  listMat=unique(listMat)
  # chaque maternite est enregistree dans une matrice avec les juveniles en y et le numero du run en x
  # et en premiere colonne l'ID du juvenile
  cat("Creation du fichier des maternites.\n")
  maternities=matrix(ncol=(nrun+1),nrow=length(listMat))
  maternities[,1]=listMat
  for (i in 1:nrun) {
    mat=read.table(paste(path,list.files(path=path,pattern="\\.Maternity")[i],sep=""),h=T,fill=TRUE)
    for (j in 1:nrow(mat)) {
      maternities[which(maternities[,1]==mat$OffspringID[j]),i+1]=mat$InferredMum1[j]
    }
  }
  cat(length(listMat),"maternites a ce stade.\n")
  
  # Exclusion des paternites en dessous du seuil fixe
  # proportion d'assignation du bestMother
  # d'abord on identifie le bestMother, celui qui revient le plus souvent
  # ensuite on calcule sa proportion d'assignation
  cat("Exclusion des maternites en dessous du seuil.\n")
  # identifier le meilleur pere : celui qui revient le plus souvent dans les assignations
  # !!!! on utilise les numeros de lignes : ils doivent concorder entre bestMother et maternities
  bestMother=c()
  for (i in 1:nrow(maternities)){
    uniq=unique(na.omit(as.character(maternities[i,2:(nrun+1)])))
    
    if(length(uniq)==0){
      bestMother[i]=NA
    }else{
      best=c(NA,0)
      for (j in 1:length(uniq)) {
        occur=length(which(maternities[i,2:(nrun+1)]==uniq[j]))
        if(occur>best[2]){
          best[1]=uniq[j]
          best[2]=occur
        }
      }
      bestMother[i]=best[1]
    }
  }
  # suppression des cellules qui ne contiennent pas le meilleur pere
  for (i in 1:nrow(maternities)) {
    for (j in 2:(nrun+1)){
      if(!is.na(maternities[i,j]) & maternities[i,j]!=bestMother[i]){
        maternities[i,j]=NA
      }
    }
  }
  # suppression des lignes qui sont en dessous du seuil d'exclusion
  # fais la liste des lignes a conserver
  index=c()
  for (i in 1:nrow(maternities)) {
    if (sum(!is.na(maternities[i,2:(nrun+1)]))/nrun > excl){
      index=c(index,i)
    }
  }
  maternities=maternities[index,]
  liste=c()
  for (i in 1:nrow(maternities)) {
    liste[i]=unique(maternities[i,2:(nrun+1)])[!is.na(unique(maternities[i,2:(nrun+1)]))]
  }
  maternities=cbind(maternities[,1],liste)
  colnames(maternities)=c("offspring","mother")
  
  cat(nrow(maternities),"maternites finalement retenues avec",length(unique(liste)),"meres.\n")
  
  cat("Ecriture du fichier des maternites.\n")
  write.table(maternities,"Maternity.txt",col.names=F,row.names=F)  
  # renvoie les paternites
  return(maternities)
}


#---------------------------------#
# RESULTS FOR ASSIGNATIONS WITH INFOS
#---------------------------------#

# Recupere les assignations inferees et reconstruites par COLONY
# Critere de selection :
# - Valeur de probability > 0.98 (0.99 ?) cf simulations (bootstrap value of the mean)

# Jointure sur le juvenile (ID juv, colonne 1)
# Deux sexes : peres et meres
# Deux groupes
# - individus inferes
# - genotypes reconstruits

# Colonnes du data frame
# 1/ Offsp ID
# 2/ Parent ID
# 3/ Sexe (male ou female)
# 4/ nb of runs where ID is found
# 5/ ID inferred ($) or reconstructed by COLONY (*)
# 6/ Colonnes 1:16 du genotype

# THEN
# Selection process of the best father
# - best father = the mostly found for each offspring
# - threshold : found at least on 30% of the runs

# ARGUMENTS
# 1/ dir =  repertoire des donnees dans le dossier output
# 2/ runs = nombre de runs sur lesquels construire les resultats
# 3/ pourcentage de runs necessaires pour conserver l'assignation

AssignationWithInfos=function(dir="",runs=30,p=0.3){
  setwd(dir)
  
  assignM=matrix(NA,ncol=21,nrow=1)
  assignM=as.data.frame(assignM)
  assignM=assignM[-1,]
  colnames(assignM)=c("offspID","parentID","sex","nbRuns","obs",rep("mk",16))
  
  for (r in 1:runs) {
    # Import des donnees pour le run r
    cat("Run", r, ": ")
    lsFiles=list.files(".")
    bestConfig=read.table(lsFiles[which(grepl(paste("[a-z]",r,".BestConfig_Ordered",sep=""),lsFiles))],header=TRUE,fill=TRUE)
    dadGenotype=read.table(lsFiles[which(grepl(paste("[a-z]",r,".DadGenotype",sep=""),lsFiles))],header=TRUE,fill=TRUE,sep=",")

    # recupere les ID des fathers/mothers pour chaque offspring
    # dans le fichier output .BestConfig
    # Incremente d'un 1 run la col nbRuns si le parent est deja affecte a ce juvenile
    # renseigne si le parent est inferre ou estime
    
    
    # A chaque nouvelle ligne, compare le genotype unique avec le nouveau genotype unique pour le juv
    # Si meme genotype -> meme individu donc incremente de 1 le nb de runs
    # Sinon, si genotype different, pere different et nouvelle ligne :
    #     si pere infere -> garde le nom
    #     si pere estime -> affecte un nouveau nom
    
    cat("Ajoute les males\n")
    for (i in 1:nrow(bestConfig)) {
      # first the father
      # ne garde que les peres echantillonnes et les peres estimes a partir d'une vraie mere
      if (!((grepl("#",bestConfig[i,3]) & grepl("\\*",bestConfig[i,2])) | (grepl("\\*",bestConfig[i,2]) & is.na(bestConfig[i,3])))){
        OM=bestConfig[i,1:2]
        if(grepl("\\*",OM[1,2])){
          genotype=as.vector(t(dadGenotype[which(dadGenotype[,1]==OM[1,2]),6:7]))
        }else{
          genotype=as.vector(t(dadGenotype[which(dadGenotype[,1]==OM[1,2]),3:4]))
        }
        
        # Juvenile deja present ?
        if(as.character(OM[1,1]) %in% assignM[,1]){
          lineOffs=which(as.character(OM[1,1])==assignM[,1])
        }else{
          lineOffs=0
        }
        # juvenile n'est pas encore present OU tableau assign est vide  -> on ajoute une ligne et une description de la nouvelle paire
        if((min(lineOffs)==0 & length(lineOffs)==1) | nrow(assignM)==0){
          assignM[nrow(assignM)+1,1:21]=c(as.character(OM[1,1]),as.character(OM[1,2]),"M",1,NA,genotype)
        }else{
          #si il est present, genotype(s) uniques(s) associe(s) = genotype donc on incremente le nb de runs
          match=0
          for (j in 1:length(lineOffs)) {
            alreadySampled=assignM[lineOffs[j],6:21]
            # si les genotypes sont identiques
            if(prod(alreadySampled==genotype)!=0 & match==0){
              match=1
              assignM[lineOffs[j],4]=as.numeric(assignM[lineOffs[j],4])+1
            }
          }
          # si le genotype unique n'a pas ete trouve dans assign, alors c'est un nouvel individu, nouvelle ligne
          if(match==0) {
            assignM[nrow(assignM)+1,1:21]=c(as.character(OM[1,1]),as.character(OM[1,2]),"M",1,NA,genotype)
          }
        }
      }
    }
  }
  
  # SELECTION
  # exclut tous les genotypes qui ne sont pas retrouves sur plus de p % des runs
  cat(nrow(assignM),"individus rassembles\nExclusion :\n")
  assignM=assignM[(as.numeric(assignM[,4])/runs)>p,]
  cat(nrow(assignM),"individus conserves au seuil",p,"\n")
  
  # Une fois le tableau assemble, renomme tous les individus estimes (\*[0-9])
  # chaque genotype unique est un individu
  # ne touche qu'aux genotypes estimes
  cat("Renomme les individus (genotype unique = 1 individu)\n")
  count_estimatedParents=0 # compteur de parents estimes pour affecter un ID unique a chaque nouveau parent estime retrouve
  gRef=c()
  g=c()
  for (i in 1:nrow(assignM)) {
    if(grepl("\\*[0-9]",assignM[i,2])){
      count_estimatedParents=count_estimatedParents+1
      assignM[i,2]=paste("e",count_estimatedParents,sep="")
      gRef=assignM[i,6:21]
      # verifie toutes les autres lignes pour trouver genotypes identiques
      for(j in 1:nrow(assignM)){
        g=assignM[j,6:21]
        if(grepl("\\*[0-9]",assignM[j,2]) & paste(g,collapse=" ")==paste(gRef,collapse=" ")){
          assignM[j,2]=paste("e",count_estimatedParents,sep="")
        }
      }
    }
  }
  for (i in 1:nrow(assignM)) {
    if(grepl("e[0-9]",assignM[i,2])){
      assignM[i,5]="E"
    }else{
      assignM[i,5]="O"
    }
  }
  cat(sum(assignM[,5]=="O"),"assignations avec parents echantillonnes et",sum(assignM[,5]=="E"),"assignations avec parents estimes\n")
  
  
  # MAINTENANT AU TOUR DES FEMELLES
  assignF=matrix(NA,ncol=21,nrow=1)
  assignF=as.data.frame(assignF)
  assignF=assignF[-1,]
  colnames(assignF)=c("offspID","parentID","sex","nbRuns","obs",rep("mk",16))
  
  for (r in 1:runs) {
    # Import des donnees pour le run r
    cat("Run", r, ": ")
    lsFiles=list.files(".")
    bestConfig=read.table(lsFiles[which(grepl(paste("[a-z]",r,".BestConfig_Ordered",sep=""),lsFiles))],header=TRUE,fill=TRUE)
    mumGenotype=read.table(lsFiles[which(grepl(paste("[a-z]",r,".MumGenotype",sep=""),lsFiles))],header=TRUE,fill=TRUE,sep=",")

    cat("Ajoute les femelles\n")
    for (i in 1:nrow(bestConfig)) {
      OF=bestConfig[i,c(1,3)]
      if(!is.na(OF[1,2])){
        if(grepl("#",OF[1,2])){
          genotype=as.vector(t(mumGenotype[which(mumGenotype[,1]==OF[1,2]),6:7]))
        }else{
          genotype=as.vector(t(mumGenotype[which(mumGenotype[,1]==OF[1,2]),3:4]))
        }
        
        # Juvenile deja present ?
        if(as.character(OF[1,1]) %in% assignF[,1]){
          lineOffs=which(OF[1,1]==assignF[,1])
        }else{
          lineOffs=0
        }
        # juvenile n'est pas encore present OU tableau assign est vide  -> on ajoute une ligne et une description de la nouvelle paire
        if(min(lineOffs)==0 | nrow(assignF)==0){
          assignF[nrow(assignF)+1,1:21]=c(as.character(OF[1,1]),as.character(OF[1,2]),"F",1,NA,genotype)
        }else{
          #si il est present, genotype(s) uniques(s) associe(s) = genotype donc on incremente le nb de runs
          match=0
          for (j in 1:length(lineOffs)) {
            alreadySampled=assignF[lineOffs[j],6:21]
            if(prod(alreadySampled==genotype)!=0 & match==0){
              match=1
              assignF[lineOffs[j],4]=as.numeric(assignF[lineOffs[j],4])+1
            }
          }
          # si le genotype unique n'a pas ete trouve dans assign, alors c'est un nouvel individu, nouvelle ligne
          if(match==0) {
            assignF[nrow(assignF)+1,1:21]=c(as.character(OF[1,1]),as.character(OF[1,2]),"F",1,NA,genotype)
          }
        }
      }
    }
  }
    
  # SELECTION
  # exclut tous les genotypes qui ne sont pas retrouves sur plus de p % des runs
  cat(nrow(assignF),"individus rassembles\nExclusion :\n")
  assignF=assignF[(as.numeric(assignF[,4])/runs)>p,]
  cat(nrow(assignF),"individus conserves au seuil",p,"\n")
  
  # Une fois le tableau assemble, renomme tous les individus estimes (\*[0-9])
  # chaque genotype unique est un individu
  # ne touche qu'aux genotypes estimes
  cat("Renomme les individus (genotype unique = 1 individu)\n")
  gRef=c()
  g=c()
  for (i in 1:nrow(assignF)) {
    if(grepl("\\#[0-9]",assignF[i,2])){
      count_estimatedParents=count_estimatedParents+1
      assignF[i,2]=paste("e",count_estimatedParents,sep="")
      gRef=assignF[i,6:21]
      # verifie toutes les autres lignes pour trouver genotypes identiques
      for(j in 1:nrow(assignF)){
        g=assignF[j,6:21]
        if(grepl("\\#[0-9]",assignF[j,2]) & paste(g,collapse=" ")==paste(gRef,collapse=" ")){
          assignF[j,2]=paste("e",count_estimatedParents,sep="")
        }
      }
    }
  }
  for (i in 1:nrow(assignF)) {
    if(grepl("e[0-9]",assignF[i,2])){
      assignF[i,5]="E"
    }else{
      assignF[i,5]="O"
    }
  }  
  cat(sum(assignF[,5]=="O"),"assignations avec parents echantillonnes et",sum(assignF[,5]=="E"),"assignations avec parents estimes\n")
  
  # Reunion peres+meres
  assign=rbind(assignM,assignF)

  # Annee de naissance du juvenile
  setwd(wd)
  uniqueGenInfos=read.table("uniqueGenotypesWithInfo.txt",header=TRUE)
  assign$yearOffsp=rep(NA,nrow(assign))
  for (i in 1:nrow(assign)) {
    assign$yearOffsp[i]=uniqueGenInfos$yearWhenFirstCaptur[which(uniqueGenInfos$idind==assign$offspID[i])]
  }
  
  cat("Enregistrement du dataset\n")
  write.table(assign,file="Assignations.txt",col.names=TRUE,row.names=FALSE)
  return(assign)

  cat("Fin\n")
}

#-----------------------------#
# MEAN NUMBER OF JUVENILES PER PARENT
#-----------------------------#
Rsuccess=function(){
  # un data frame contenant pour chaque categorie (breeding, sampled/unsampled, male/female) :
  # - effectif n
  # - moyenne du succes reproductif
  # - ecart type
  # - IC 95%
  # - annee
  Rsuccess=data.frame(n=NA,year=NA,YearlyMeanReproductiveSuccess=NA,SEM=NA,lowerIC=NA,upperIC=NA)
  
  # Nombre de juveniles/an par pere s'etant reproduit
  # un vecteur contenant les peres et le nombre de juveniles sur chaque annee (2013-2016)
  ReprodSuccess=data.frame(FatherID=unique(assignations[which(assignations[,3]=="M"),2]),
                           NbOffspring2013=rep(NA,length(unique(assignations[which(assignations[,3]=="M"),2]))),
                           NbOffspring2014=rep(NA,length(unique(assignations[which(assignations[,3]=="M"),2]))),
                           NbOffspring2015=rep(NA,length(unique(assignations[which(assignations[,3]=="M"),2]))),
                           NbOffspring2016=rep(NA,length(unique(assignations[which(assignations[,3]=="M"),2]))))
  for (i in 1:nrow(ReprodSuccess)) {
    ReprodSuccess[i,2]=length(which(assignations[,2]==ReprodSuccess[i,1] & assignations[,4]=="2013"))
    ReprodSuccess[i,3]=length(which(assignations[,2]==ReprodSuccess[i,1] & assignations[,4]=="2014"))
    ReprodSuccess[i,4]=length(which(assignations[,2]==ReprodSuccess[i,1] & assignations[,4]=="2015"))
    ReprodSuccess[i,5]=length(which(assignations[,2]==ReprodSuccess[i,1] & assignations[,4]=="2016"))
  }
  
  # Succes reproductif sur les peres echantillonnes UNIQUEMENT
  # ReprodSuccessSampled=ReprodSuccess[!grepl("e[0-9]",ReprodSuccess[,1]),]
  
  # BOOTSTRAP REPRODUCTIVE SUCCESS - BREEDING MALES (SAMPLED ONLY)
  # YEARS 2013-2016
  boots=numeric(10000)
  for (i in 1:10000){
    boots[i]=mean(sample(as.vector(t(ReprodSuccess[,2:5])),replace=T),na.rm=TRUE)
  }
  limBoot1=quantile(boots,c(.025,.975))
  limBoot1
  mboot=mean(boots)
  mboot
  meanReprodSuccess=mboot
  etboot=sqrt(var(boots))
  etboot
  Rsuccess[1,]=c(nrow(ReprodSuccess),"MalesAllSampled",mboot,etboot,limBoot1[1],limBoot1[2])
  # 2013
  boots=numeric(10000)
  for (i in 1:10000){
    boots[i]=mean(sample(as.vector(t(ReprodSuccess[,2])),replace=T),na.rm=TRUE)
  }
  limBoot1=quantile(boots,c(.025,.975))
  limBoot1
  mboot=mean(boots)
  mboot
  meanReprodSuccess=mboot
  etboot=sqrt(var(boots))
  etboot
  Rsuccess[2,]=c(nrow(ReprodSuccess),"Males2013Sampled",mboot,etboot,limBoot1[1],limBoot1[2])
  # 2014
  boots=numeric(10000)
  for (i in 1:10000){
    boots[i]=mean(sample(as.vector(t(ReprodSuccess[,3])),replace=T),na.rm=TRUE)
  }
  limBoot1=quantile(boots,c(.025,.975))
  limBoot1
  mboot=mean(boots)
  mboot
  meanReprodSuccess=mboot
  etboot=sqrt(var(boots))
  etboot
  Rsuccess[3,]=c(nrow(ReprodSuccess),"Males2014Sampled",mboot,etboot,limBoot1[1],limBoot1[2])
  # 2015
  boots=numeric(10000)
  for (i in 1:10000){
    boots[i]=mean(sample(as.vector(t(ReprodSuccess[,4])),replace=T),na.rm=TRUE)
  }
  limBoot1=quantile(boots,c(.025,.975))
  limBoot1
  mboot=mean(boots)
  mboot
  meanReprodSuccess=mboot
  etboot=sqrt(var(boots))
  etboot
  Rsuccess[4,]=c(nrow(ReprodSuccess),"Males2015Sampled",mboot,etboot,limBoot1[1],limBoot1[2])
  # 2016
  boots=numeric(10000)
  for (i in 1:10000){
    boots[i]=mean(sample(as.vector(t(ReprodSuccess[,5])),replace=T),na.rm=TRUE)
  }
  limBoot1=quantile(boots,c(.025,.975))
  limBoot1
  mboot=mean(boots)
  mboot
  meanReprodSuccess=mboot
  etboot=sqrt(var(boots))
  etboot
  Rsuccess[5,]=c(nrow(ReprodSuccess),"Males2016Sampled",mboot,etboot,limBoot1[1],limBoot1[2])
  
  
  # BOOTSTRAP REPRODUCTIVE SUCCESS - BREEDING MALES (SAMPLED & ESTIMATED)
  # YEARS 2013-2016
  # boots=numeric(10000)
  # for (i in 1:10000){
  #   boots[i]=mean(sample(as.vector(t(ReprodSuccess[,2:5])),replace=T),na.rm=TRUE)
  # }
  # limBoot1=quantile(boots,c(.025,.975))
  # limBoot1
  # mboot=mean(boots)
  # mboot
  # meanReprodSuccess=mboot
  # etboot=sqrt(var(boots))
  # etboot
  # Rsuccess[6,]=c(nrow(ReprodSuccess),"MalesAllEstimated",mboot,etboot,limBoot1[1],limBoot1[2])
  # # 2013
  # boots=numeric(10000)
  # for (i in 1:10000){
  #   boots[i]=mean(sample(as.vector(t(ReprodSuccess[,2])),replace=T),na.rm=TRUE)
  # }
  # limBoot1=quantile(boots,c(.025,.975))
  # limBoot1
  # mboot=mean(boots)
  # mboot
  # meanReprodSuccess=mboot
  # etboot=sqrt(var(boots))
  # etboot
  # Rsuccess[7,]=c(nrow(ReprodSuccess),"Males2013Estimated",mboot,etboot,limBoot1[1],limBoot1[2])
  # # 2014
  # boots=numeric(10000)
  # for (i in 1:10000){
  #   boots[i]=mean(sample(as.vector(t(ReprodSuccess[,3])),replace=T),na.rm=TRUE)
  # }
  # limBoot1=quantile(boots,c(.025,.975))
  # limBoot1
  # mboot=mean(boots)
  # mboot
  # meanReprodSuccess=mboot
  # etboot=sqrt(var(boots))
  # etboot
  # Rsuccess[8,]=c(nrow(ReprodSuccess),"Males2014Estimated",mboot,etboot,limBoot1[1],limBoot1[2])
  # # 2015
  # boots=numeric(10000)
  # for (i in 1:10000){
  #   boots[i]=mean(sample(as.vector(t(ReprodSuccess[,4])),replace=T),na.rm=TRUE)
  # }
  # limBoot1=quantile(boots,c(.025,.975))
  # limBoot1
  # mboot=mean(boots)
  # mboot
  # meanReprodSuccess=mboot
  # etboot=sqrt(var(boots))
  # etboot
  # Rsuccess[9,]=c(nrow(ReprodSuccess),"Males2015Estimated",mboot,etboot,limBoot1[1],limBoot1[2])
  # # 2016
  # boots=numeric(10000)
  # for (i in 1:10000){
  #   boots[i]=mean(sample(as.vector(t(ReprodSuccess[,5])),replace=T),na.rm=TRUE)
  # }
  # limBoot1=quantile(boots,c(.025,.975))
  # limBoot1
  # mboot=mean(boots)
  # mboot
  # meanReprodSuccess=mboot
  # etboot=sqrt(var(boots))
  # etboot
  # Rsuccess[10,]=c(nrow(ReprodSuccess),"Males2016Estimated",mboot,etboot,limBoot1[1],limBoot1[2])
  # 
  #-------------------------------
  # FEMALES
  ReprodSuccessF=data.frame(FatherID=unique(assignations[which(assignations[,3]=="F"),2]),
                            NbOffspring2013=rep(NA,length(unique(assignations[which(assignations[,3]=="F"),2]))),
                            NbOffspring2014=rep(NA,length(unique(assignations[which(assignations[,3]=="F"),2]))),
                            NbOffspring2015=rep(NA,length(unique(assignations[which(assignations[,3]=="F"),2]))),
                            NbOffspring2016=rep(NA,length(unique(assignations[which(assignations[,3]=="F"),2]))))
  for (i in 1:nrow(ReprodSuccessF)) {
    ReprodSuccessF[i,2]=length(which(assignations[,2]==ReprodSuccessF[i,1] & assignations[,4]=="2013"))
    ReprodSuccessF[i,3]=length(which(assignations[,2]==ReprodSuccessF[i,1] & assignations[,4]=="2014"))
    ReprodSuccessF[i,4]=length(which(assignations[,2]==ReprodSuccessF[i,1] & assignations[,4]=="2015"))
    ReprodSuccessF[i,5]=length(which(assignations[,2]==ReprodSuccessF[i,1] & assignations[,4]=="2016"))
  }
  
  # Succes reproductif sur les peres echantillonnes UNIQUEMENT
  # ReprodSuccessSampledF=ReprodSuccessF[!grepl("e[0-9]",ReprodSuccessF[,1]),]
  
  # BOOTSTRAP REPRODUCTIVE SUCCESS - BREEDING MALES (SAMPLED ONLY)
  # YEARS 2013-2016
  boots=numeric(10000)
  for (i in 1:10000){
    boots[i]=mean(sample(as.vector(t(ReprodSuccessF[,2:5])),replace=T),na.rm=TRUE)
  }
  limBoot1=quantile(boots,c(.025,.975))
  limBoot1
  mboot=mean(boots)
  mboot
  meanReprodSuccess=mboot
  etboot=sqrt(var(boots))
  etboot
  Rsuccess[6,]=c(nrow(ReprodSuccessF),"FemalesAllSampled",mboot,etboot,limBoot1[1],limBoot1[2])
  # 2013
  boots=numeric(10000)
  for (i in 1:10000){
    boots[i]=mean(sample(as.vector(t(ReprodSuccessF[,2])),replace=T),na.rm=TRUE)
  }
  limBoot1=quantile(boots,c(.025,.975))
  limBoot1
  mboot=mean(boots)
  mboot
  meanReprodSuccess=mboot
  etboot=sqrt(var(boots))
  etboot
  Rsuccess[7,]=c(nrow(ReprodSuccessF),"Females2013Sampled",mboot,etboot,limBoot1[1],limBoot1[2])
  # 2014
  boots=numeric(10000)
  for (i in 1:10000){
    boots[i]=mean(sample(as.vector(t(ReprodSuccessF[,3])),replace=T),na.rm=TRUE)
  }
  limBoot1=quantile(boots,c(.025,.975))
  limBoot1
  mboot=mean(boots)
  mboot
  meanReprodSuccess=mboot
  etboot=sqrt(var(boots))
  etboot
  Rsuccess[8,]=c(nrow(ReprodSuccessF),"Females2014Sampled",mboot,etboot,limBoot1[1],limBoot1[2])
  # 2015
  boots=numeric(10000)
  for (i in 1:10000){
    boots[i]=mean(sample(as.vector(t(ReprodSuccessF[,4])),replace=T),na.rm=TRUE)
  }
  limBoot1=quantile(boots,c(.025,.975))
  limBoot1
  mboot=mean(boots)
  mboot
  meanReprodSuccess=mboot
  etboot=sqrt(var(boots))
  etboot
  Rsuccess[9,]=c(nrow(ReprodSuccessF),"Females2015Sampled",mboot,etboot,limBoot1[1],limBoot1[2])
  # 2016
  boots=numeric(10000)
  for (i in 1:10000){
    boots[i]=mean(sample(as.vector(t(ReprodSuccessF[,5])),replace=T),na.rm=TRUE)
  }
  limBoot1=quantile(boots,c(.025,.975))
  limBoot1
  mboot=mean(boots)
  mboot
  meanReprodSuccess=mboot
  etboot=sqrt(var(boots))
  etboot
  Rsuccess[10,]=c(nrow(ReprodSuccessF),"Females2016Sampled",mboot,etboot,limBoot1[1],limBoot1[2])
  
  
  # BOOTSTRAP REPRODUCTIVE SUCCESS - BREEDING MALES (SAMPLED & ESTIMATED)
  # YEARS 2013-2016
  # boots=numeric(10000)
  # for (i in 1:10000){
  #   boots[i]=mean(sample(as.vector(t(ReprodSuccessF[,2:5])),replace=T),na.rm=TRUE)
  # }
  # limBoot1=quantile(boots,c(.025,.975))
  # limBoot1
  # mboot=mean(boots)
  # mboot
  # meanReprodSuccess=mboot
  # etboot=sqrt(var(boots))
  # etboot
  # Rsuccess[16,]=c(nrow(ReprodSuccessF),"FemalesAllEstimated",mboot,etboot,limBoot1[1],limBoot1[2])
  # # 2013
  # boots=numeric(10000)
  # for (i in 1:10000){
  #   boots[i]=mean(sample(as.vector(t(ReprodSuccessF[,2])),replace=T),na.rm=TRUE)
  # }
  # limBoot1=quantile(boots,c(.025,.975))
  # limBoot1
  # mboot=mean(boots)
  # mboot
  # meanReprodSuccess=mboot
  # etboot=sqrt(var(boots))
  # etboot
  # Rsuccess[17,]=c(nrow(ReprodSuccessF),"Females2013Estimated",mboot,etboot,limBoot1[1],limBoot1[2])
  # # 2014
  # boots=numeric(10000)
  # for (i in 1:10000){
  #   boots[i]=mean(sample(as.vector(t(ReprodSuccessF[,3])),replace=T),na.rm=TRUE)
  # }
  # limBoot1=quantile(boots,c(.025,.975))
  # limBoot1
  # mboot=mean(boots)
  # mboot
  # meanReprodSuccess=mboot
  # etboot=sqrt(var(boots))
  # etboot
  # Rsuccess[18,]=c(nrow(ReprodSuccessF),"Females2014Estimated",mboot,etboot,limBoot1[1],limBoot1[2])
  # # 2015
  # boots=numeric(10000)
  # for (i in 1:10000){
  #   boots[i]=mean(sample(as.vector(t(ReprodSuccessF[,4])),replace=T),na.rm=TRUE)
  # }
  # limBoot1=quantile(boots,c(.025,.975))
  # limBoot1
  # mboot=mean(boots)
  # mboot
  # meanReprodSuccess=mboot
  # etboot=sqrt(var(boots))
  # etboot
  # Rsuccess[19,]=c(nrow(ReprodSuccessF),"Females2015Estimated",mboot,etboot,limBoot1[1],limBoot1[2])
  # # 2016
  # boots=numeric(10000)
  # for (i in 1:10000){
  #   boots[i]=mean(sample(as.vector(t(ReprodSuccessF[,5])),replace=T),na.rm=TRUE)
  # }
  # limBoot1=quantile(boots,c(.025,.975))
  # limBoot1
  # mboot=mean(boots)
  # mboot
  # meanReprodSuccess=mboot
  # etboot=sqrt(var(boots))
  # etboot
  # Rsuccess[20,]=c(nrow(ReprodSuccessF),"Females2016Estimated",mboot,etboot,limBoot1[1],limBoot1[2])
  # 
  return(Rsuccess)
}


#-----------------------------#
# SIM RESULTS WITH INFOS - EXPORT DATA
#-----------------------------
# give a file containing :
# - offsrping ID
# - fathers and mothers ID
# - sex
# - correct/incorrect assignation
# - genotype

empiricalResultsWithInfo=function(){
  
  paternity=constructResults(path="outputs/Colony79lociSansM399/",excl=0.3)
  paternity2=constructResults(path="outputs/Colony99loci/",excl=0.3)
  Paternities=combineParentageMethods("all",sex="M")
  descResultsM=describeResults(Paternities)
  
  maternity=constructResultsMums(path="outputs/Colony79lociSansM399/",excl=0.3)
  maternity2=constructResultsMums(path="outputs/Colony99loci/",excl=0.3)
  Maternities=combineParentageMethods("all",sex='F')
  descResultsF=describeResults(Maternities)
  
  dyads=data.frame(offspring=unique(c(descResultsM$offspring,descResultsF$offspring)),
                      father=rep(NA,length(unique(c(descResultsM$offspring,descResultsF$offspring)))),
                      mother=rep(NA,length(unique(c(descResultsM$offspring,descResultsF$offspring)))))
  for (i in 1:nrow(dyads)) {
    if(sum(descResultsM$offspring==dyads$offspring[i])>0){
      dyads$father[i]=descResultsM[which(descResultsM$offspring==dyads$offspring[i]),2]
    }
    if (sum(descResultsF$offspring==dyads$offspring[i])>0){
      dyads$mother[i]=descResultsF[which(descResultsF$offspring==dyads$offspring[i]),2]
    }
  }
  
  # Add the genotype of the parent
  # columns 5 to 20  are for father genotype
  # columns 21 to 36 are for mother genotype
  dyads[,4:35]=rep(NA,nrow(dyads))
  for (i in 1:nrow(dyads)) {
    if (!is.na(dyads$father[i])){
      dyads[i,4:19]=info[which(info[,1]==dyads$father[i]),2:17]
    }
    if (!is.na(dyads$mother[i])){
      dyads[i,20:35]=info[which(info[,1]==dyads$mother[i]),2:17]
    }
  }
  
  
  # Export
  write.table(dyads,paste("Pic_resultsWithInfo.txt",sep=""),col.names=TRUE,row.names=FALSE)
  return(dyads)
}



#-------------------------------#
#   CLEAN DIRECTORIES & REMOVE USELESS FILES
#-------------------------------
# remove results files from COLONY which are useless for us
# e.g. output files of  .GtypeData, .MidResults...
# in a specified directory (path = ...)
cleanUp=function(path=""){
  listOfErasal=c("\\.ConfigArchive","\\.BestFSFamily","\\.GtypeData","\\.MidResult","\\.BestPrntCND","\\.Distribution","\\.Ne","\\.OffGenotype","\\.AlleleFreq")
  for (i in 1:length(listOfErasal)) {
    index=list.files(path=path,pattern=listOfErasal[i])
    file.remove(paste(datadir,path,index,sep=""))
  }
}