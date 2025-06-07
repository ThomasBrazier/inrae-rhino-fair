#===============================#
#     Mantel by permutations
#     
#===============================#

#------------------
# MANTEL TEST

manteltest = function(nperm=9999){
  
}


#------------------
# CORRELOGRAMME
# correlogramme sur les distances genetiques et geographiques (km, non log)
# fournit IC, mantel r (permutation) + error bars et classes de distance
# sur distances reelles, correlation genetique en fonction de chaque couple de colonie (distance)
# d'apres Diniz-Filho et al. 2013
# matrice de distances genetiques G
# matrice de distances geog. D

# D divisee en plusieurs sous-matrices, chacune comprenant paires de populations dans un intervalle de distances
# sous-matrices Wk decrivent les populations connectees par un lien binaire
# 0 non connecte / 1 connecte
# k parametre de distance class -> nombre de classes
# statistiquement et idealement, meme nombre de paires de distance dans chaque classe.

# prend en arguments :
# 1/ matrice de distances genetiques
# 2/ matrice de distances geographiques
# 3/ nombre de classes de distance
# 4/ nombre de permutations
mantel.correlogram=function(G,D,nclass=5,nperm=9999){
  library(ecodist)
  
  # division de D en sous matrices Wk avec k le nombre de classes
  k=nclass
  r=c() # stockage des r pour chaque classe
  up=c() # upper IC
  low=c() # lower IC
  meandist=c()
  p=c() # two tailed p value under hypothesis r=0
  upH0=c() # borne sup IC95% de HO
  lowH0=c() # borne inf IC95% de HO
  
  # meme nombre de paires de distance dans chaque classe
  # tri des paires de distance != 0 par ordre croissant
  # divise la liste des distances triees en k classes de meme taille
  distance=which(D!=0)
  # liste des indices par ordre croissant de distance
  d=D[distance]
  d=sort(d)
  index=c()
  for (j in 1:length(d)) {
    index=c(index,which(D==d[j]))
  }
  index=unique(index)
  W=D[index] # distances geo classees par ordre croissant dans la matrice D
  V=G[index] # distances gen de meme
  boundaries=ceiling(seq(1,length(distance)+1,length.out=k+1)) # limites des classes en fonction du nombre de paires de distance
  # calcul de Mantel r pour chaque classe k
  for (i in 1:k) {
    # Wk=W[boundaries[i]:(boundaries[i+1]-1)]
    # Vk=V[boundaries[i]:(boundaries[i+1]-1)]
    # # passage des vecteurs Wk t Vk en matrices de dissimilarites
    # # matrices triangulaires inferieures, sans la diagonale
    # dimMatrix=
    # Wdist=matrix(0,nrow=dimMatrix,ncol=dimMatrix) # initialise matrices vides
    # Vdist=matrix(0,nrow=dimMatrix,ncol=dimMatrix) # initialise matrices vides
    # Ajouter les valeurs dans les demi matrices inferieures
    meandist[i]=mean(c(W[boundaries[i]],W[(boundaries[i+1]-1)]))
    
    Wk=matrix(1,dim(D)[1],dim(D)[2])
    Wk[index[boundaries[i]:(boundaries[i+1]-1)]]=0
    
    # bootstrap pour avoir estimation du r de mantel avec IC ?
    # calcul du r de Mantel avec test maison
    #r[i]=manteltest(G,W,nperm=10000)$r
    
    # utilisation du package ecodist
    res=mantel(as.dist(G)~as.dist(Wk),nperm = 10000)
    r[i]=res[1]
    up[i]=res[5]
    low[i]=res[6]
    p[i]=res[2]

    # Calcul de H0 pour chaque classe de distance
    # Hypothese : distance genetique au hasard, parmi n valeurs tirees au hasard
    # ou n = nombre de couples dans chaque classe
    boots=numeric(1000)
    for (l in 1:1000){
      H0=matrix(1,dim(D)[1],dim(D)[2])
      H0[sample(index,boundaries[2],replace=TRUE)]=0
      boots[l]=mantel(as.dist(G)~as.dist(H0),nperm = 10000)[1]
    }
    upH0[i]=quantile(boots,c(.975))
    lowH0[i]=quantile(boots,c(.025))
  }

  corr=data.frame(r=r,dist.class=meandist,upperIC=up,lowerIC=low,p=p,upH0=upH0,lowH0=lowH0)
  return(corr)
}

