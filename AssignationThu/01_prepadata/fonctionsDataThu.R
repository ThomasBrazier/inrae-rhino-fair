#-------------------------------#
# FONCTIONS
#-------------------------------#


#-------------------------------#
# CONVERT GERMAN FILE TO CORRECT FORMAT
#-------------------------------
# Convertit le fichier All_Genotypes_FINAL_man.csv en uniqueGenotypesWithInfos.txt
convertToUnique=function(falseloci=1,sexmark=FALSE,withcolony=TRUE){
  timeIn=Sys.time()
  
  infile=read.table("All_Genotypes_FINAL_man.csv",header=TRUE,sep=",")
  
  # individus presents en double
  # on recupere une liste des individus uniques
  index=sort(unique(infile$indID))
  
  # Structure : Names (idind), Alleles, colonie (idcol), sexe, age premiere capture (juv,adult), annee de premiere capture
  uniqueGenotypesWithInfo=matrix(nrow=length(index),ncol=29)
  colnames(uniqueGenotypesWithInfo)=c("idind","dby8","dby8b","rha101","rha101b","rha109","rha109b","rha4","rha4b","rha7",
                                      "rha7b","rhc108","rhc108b","rhc3","rhc3b","rhd102","rhd102b","rhd103","rhd103b",
                                      "idcol","sexe","ageWhenFirstCaptur","yearWhenFirstCaptur","col2015_1","col2015_2",
                                      "col2016_1","col2016_2","col2017_1","col2017_2")
  uniqueGenotypesWithInfo=as.data.frame(uniqueGenotypesWithInfo)

  # ajoute les infos ind par ind (ligne par ligne)
  for (i in 1:length(index)) {
    l=min(which(infile$indID==index[i]))
    uniqueGenotypesWithInfo$idind[i]=infile$indID[l]
    if(infile$sess[l]==1){
      uniqueGenotypesWithInfo$ageWhenFirstCaptur[i]="Adult"
    }else{
      uniqueGenotypesWithInfo$ageWhenFirstCaptur[i]="Juv"
    }
    uniqueGenotypesWithInfo$idcol[i]=as.character(infile$pop[l])
    uniqueGenotypesWithInfo$yearWhenFirstCaptur[i]=infile$year[l]
    # recree un vecteur du genotype pour le mettre au bon format
    g=infile[l,7:15]
    vec=c()
    for (j in 1:9) {
      if(is.na(g[j])){
        vec[seq(1,17,2)[j]:seq(2,18,2)[j]]=0
      }else{
        vec[seq(1,17,2)[j]]=as.numeric(substr(as.character(unlist(g[j])),1,3))
        vec[seq(2,18,2)[j]]=as.numeric(substr(as.character(unlist(g[j])),5,7))
      }
    }
    uniqueGenotypesWithInfo[i,2:19]=vec
    
  }
  #-------------------------------#
  # Reconstruire historique d'echantillonnage
  #-------------------------------
  # pour chaque individu, on cherche quelles annees et quelles sessions il a ete capture
  for (i in 1:nrow(uniqueGenotypesWithInfo)) {
    index=(which(infile$indID==uniqueGenotypesWithInfo$idind[i]))
    for (j in 1:length(index)) {
      if(infile$year[index[j]]==2015){
        if(1==infile$sess[index[j]]){
          uniqueGenotypesWithInfo[i,24]=as.character(infile$pop[index[j]])
        }else{
          uniqueGenotypesWithInfo[i,25]=as.character(infile$pop[index[j]])
        }
      }
      if(infile$year[index[j]]==2016){
        if(1==infile$sess[index[j]]){
          uniqueGenotypesWithInfo[i,26]=as.character(infile$pop[index[j]])
        }else{
          uniqueGenotypesWithInfo[i,27]=as.character(infile$pop[index[j]])
        }
      }
      if(infile$year[index[j]]==2017){
        if(1==infile$sess[index[j]]){
          uniqueGenotypesWithInfo[i,28]=as.character(infile$pop[index[j]])
        }else{
          uniqueGenotypesWithInfo[i,29]=as.character(infile$pop[index[j]])
        }
      }
    }
  }
  
  
  #-------------------------------#
  # Enlever les donnees manquantes
  #-------------------------------
  
  # Option 1
  # On garde l'individu quand (9-"falseloci") loci ont ete correctement genotypes : 7/9 d'apres Jan (2017)
  # grouper les alleles deux a deux (produit des deux alleles)
  # si 1 des alleles=0 -> FALSE (0), si  les 2 alleles != 0 -> TRUE (1)
  # si le nombre de loci faux pour une ligne > ou = "falseloci" alors le genotype est rejete
  loci=uniqueGenotypesWithInfo[,seq(2,18,2)]*uniqueGenotypesWithInfo[,seq(3,19,2)]
  uniqueGenotypesWithInfo=uniqueGenotypesWithInfo[(rowSums(apply(loci,c(1,2),function(x) 0 %in% x),na.rm=TRUE))<=falseloci,]
  

  #-------------------------------#
  # Identifier les femelles
  #-------------------------------
  
  # D'abord nettoyer le jeu de donnees en enlevant les ind sans marqueur sexuel (dbY8 ou dbY8b == 0)
  # en prenant la precaution de verifier si il existe des valeurs nulles dans les colonnes concernees
  if(sum(uniqueGenotypesWithInfo$dby8==0)!=0 | sum(uniqueGenotypesWithInfo$dbyb8==0)!=0){
    uniqueGenotypesWithInfo=uniqueGenotypesWithInfo[-which(uniqueGenotypesWithInfo$dby8==0 | uniqueGenotypesWithInfo$dby8b==0),]
  }
  
  # Option 1
  # Si l'individu n'a aucun de ces trois alleles, c'est une femelle (Zarsozo-Lacoste 2017)
  #colf=col[col$dby8!=128&col$dby8!=131&col$dby8!=133,]
  # uniqueGenotypesWithInfo$sexe[which(uniqueGenotypesWithInfo$dby8!=128&uniqueGenotypesWithInfo$dby8!=131&uniqueGenotypesWithInfo$dby8!=133)]="F"
  
  # Option 2
  # Les femelles ont une difference absolue entre leurs deux alleles inferieure a 10
  uniqueGenotypesWithInfo$sexe[which(abs(uniqueGenotypesWithInfo$dby8-uniqueGenotypesWithInfo$dby8b)<10)]="F"
  
  #-------------------------------#
  # Identifier les males
  #------------------------------- 
  ###pour ne prendre que les males
  
  # Option 1
  # Si l'individu possede un de ces trois alleles, c'est un male (Zarsozo-Lacoste 2017)
  #colm=col[col$dby8==128|col$dby8==131|col$dby8==133,]
  # uniqueGenotypesWithInfo$sexe[which(uniqueGenotypesWithInfo$dby8==128|uniqueGenotypesWithInfo$dby8==131|uniqueGenotypesWithInfo$dby8==133)]="M"
  
  # Option 2
  uniqueGenotypesWithInfo$sexe[which(abs(uniqueGenotypesWithInfo$dby8-uniqueGenotypesWithInfo$dby8b)>=10)]="M"
  
  #-------------------------------#
  # Enlever les dates d'echantillonnage
  #------------------------------- 
  if (withcolony==F) {
    uniqueGenotypesWithInfo=uniqueGenotypesWithInfo[,-c(24:29)]
  }
  
  #-------------------------------#
  # Enlever le marqueur sexuel dbY8
  #-------------------------------
  if (sexmark==F) {
    uniqueGenotypesWithInfo=uniqueGenotypesWithInfo[,-c(2,3)]
  }
  
  #-------------------------------#
  # Ecriture du fichier
  #-------------------------------
  write.table(uniqueGenotypesWithInfo,paste(wd,"uniqueGenotypesWithInfo.txt",sep=""),row.names=F,col.names=T,quote=F,sep="\t")
  
  timeOut=Sys.time()
  cat("Fichier des genotypes uniques construit en",difftime(timeOut,timeIn),"secondes.")
  
  return(uniqueGenotypesWithInfo)
  
}