excluded=function(){
  unique=read.table("uniqueGenotypesWithInfo.txt",h=T)
  
  # i. Excluded Maternities
  # I write lists of individuals that cannot be the mother of the focal offspring
  # ie, all females born the same or subsequent years for a given offspring
  # The lists for the two options (Cohorts and Alltogether) are not the same because the offspring lists are not the same
  # As explained in the ColonyUserGuide (p. 25), each row is made of 
  # The ID of the focal offspring, the number of excluded potential mothers, followed by their IDs

  # on travaille sur les juveniles
 
    index=which(unique$ageWhenFirstCaptur=="Juv")
  for (i in 1:length(index)){
    temp=as.vector(unique[index[i],1])
    
    y=unique[index[i],21]
    col=as.character(unique$idcol[index[i]])
    
# Les femelles nees en 2016 ne sont pas prises en compte !!!!
      # Exclusion des femelles nees apres ou meme annee y
      if (y!=2016) {
        for (j in y:2015) {
          temp=c(temp,as.vector(unique[which(unique$sexe=="F" & unique$idcol==col & unique$ageWhenFirstCaptur=="Juv" & unique$yearWhenFirstCaptur==j),1]))
        }
        itself=which(temp==as.vector(unique[index[i],1]))
        if (length(itself)==2) temp=temp[-itself[2]]
      }
      
      # exclusion de toutes les femelles issues d'une autre colonie
      for (j in 2013:2015) {
        temp=c(temp,as.vector(unique[which(unique$sexe=="F" & unique$idcol!=col & unique$yearWhenFirstCaptur==j),1]))
      }
      itself=which(temp==as.vector(unique[index[i],1]))
      if (length(itself)==2) temp=temp[-itself[2]]
      
      # Exclusion des femelles non juveniles d'une autre colonie nees en 2016
      temp=c(temp,as.vector(unique[which(unique$sexe=="F" & unique$idcol!=col & unique$yearWhenFirstCaptur==2016 & unique$ageWhenFirstCaptur!="Juv"),1]))
      itself=which(temp==as.vector(unique[index[i],1]))
      if (length(itself)==2) temp=temp[-itself[2]]

    
    # temp=c(temp,as.vector(unique[which(unique$sexe=="F" & unique$idcol!=col),1]))
    # itself=which(temp==as.vector(unique[index[i],1]))
    # if (length(itself)==2) temp=temp[-itself[2]]
    
    # enlever les doublons
    #temp=c(temp[1],unique(temp[-1]))
    
    n=length(temp)-1
    temp=c(temp[1],n,temp[-1])
    write.table(as.matrix(t(temp)),"Rhino_ExcludedMothers.txt",row.names=F,col.names=F,quote=F,sep=" ",append=T)
  }
  
   
  # ii. Excluded Paternities
  # I write lists of individuals that cannot be the father of the focal offspring
  # ie, all males born the same or subsequent years for a given offspring
  # The lists for the two options (Cohorts and Alltogether) are not the same because the offspring lists are not the same
  # As explained in the ColonyUserGuide (p. 25), each row is made of 
  # The ID of the focal offspring, the number of excluded potential fathers, followed by their IDs
  
  index=which(unique$ageWhenFirstCaptur=="Juv")
  for (i in 1:length(index)){
    y=unique[index[i],21]
    if (y!=2016) {
      temp=as.vector(unique[index[i],1])
      for (j in y:2015) {
        temp=c(temp,as.vector(unique[which(unique$sexe=="M" & unique$ageWhenFirstCaptur=="Juv" & unique$yearWhenFirstCaptur==j),1]))
      }
      itself=which(temp==as.vector(unique[index[i],1]))
      if (length(itself)==2) temp=temp[-itself[2]]
      n=length(temp)-1
      temp=c(temp[1],n,temp[-1])
      write.table(as.matrix(t(temp)),"Rhino_ExcludedFathers.txt",row.names=F,col.names=F,quote=F,sep=" ",append=T)
    }
  }
  
  # iii. Excluded Maternal Sibships
  # Note that this list are the same for the two options (Cohorts and Alltogether)
  # I write lists of offspring that cannot be the maternal sibs of the focal offspring
  # ie, all juveniles born the same year for a given offspring
  # As explained in the ColonyUserGuide (p. 25), each row is made of 
  # The ID of the focal offspring, the number of excluded maternal sibs, followed by their IDs
  
  for (i in 1:length(levels(as.factor(unique$yearWhenFirstCaptur))))	{
    OFS=unique[which(unique$ageWhenFirstCaptur=="Juv" & unique$yearWhenFirstCaptur==levels(as.factor(unique$yearWhenFirstCaptur))[i]),1]
    temp=as.data.frame(matrix(,length(OFS),length(OFS)))
    temp[,1]=as.character(OFS,quote=F)
    for (j in 1:(length(OFS)))	{
      temp[j,c(2:length(OFS))]=temp[-j,1]
    }
    numExc=rep(length(OFS)-1,length(OFS))
    temp=cbind(temp[,1],numExc,temp[,-1])
    write.table(temp,"Rhino_ExcludedMaternalSibs.txt",row.names=F,col.names=F,quote=F,sep=" ",append=T)
  }
  
  
}