#==========================================================#
#            The geometry of gametic dispersal
#               Figures for Article
#==========================================================#

#==========================================================#
# LOADING ENVIRONMENT
#==========================================================#

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


#==========================================================#
# FIGURE 1. SIMULATIONS - Assignment to the best parent ----
#==========================================================#
# Figure 1. Assignment frequency of the inferred parent according to the number of runs in simulated datasets for 2 groups:
#   black dot when the true parent was found or grey dot when the wrong parent was found (n = 3892 simulated offspring).
# (a) Paternities and (b) maternities in a simulated dataset with allelic frequencies of Picardy or Thuringia (respectively c and d).
# For each number of runs, 100 frequencies were estimated on r runs randomly sampled.
# The dashed line is the exclusion threshold of 0.3 or 0.25 according to simulations.

#-----------------------------#
# Assignation frequencies of the best father & mother in Pic79
#-----------------------------#
distOcc=read.table("AssignationPic/outputs/occurencesBestFatherIncompleteGen.txt",
                   header = T)
colnames(distOcc)=c("runs","occurence","father")
dfT=distOcc[which(distOcc$father==TRUE),1:2]
dfF=distOcc[which(distOcc$father==FALSE),1:2]
colnames(dfT)=c("runs","occurence")
colnames(dfF)=c("runs","occurence")
pMPic=ggplot(data=dfT, aes(x=runs, y=occurence)) +
  geom_point(colour="Black",size=1)+
  geom_point(data=dfF,colour="darkgrey",size=1) +
  geom_hline(yintercept=0.3,linetype="dashed",size=1) +
  scale_y_continuous(breaks=c(0,0.2,0.3,0.4,0.6),limits=c(0,0.7)) +
  # annotate(geom="text", x=1.2, y=0.7, label="(a)",size=6) +
  # ylim(0,0.7) +
  # ggtitle("Assignation frequencies of the best father\n",
  # subtitle=paste("n=1200",sep = "")) +
  xlab("") + ylab("Assignment frequency") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=14,hjust = 0.5),
        axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=14),
        axis.text=element_text(size=14, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        #legend.key.height = unit(2,"line"),
        #legend.key.width = unit(5,"line"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14))
pMPic

distOcc=read.table("AssignationPic/outputs/occurencesBestMotherIncompleteGen.txt",
                   header = T)
colnames(distOcc)=c("runs","occurence","father")
dfT=distOcc[which(distOcc$father==TRUE),1:2]
dfF=distOcc[which(distOcc$father==FALSE),1:2]
colnames(dfT)=c("runs","occurence")
colnames(dfF)=c("runs","occurence")
pFPic=ggplot(data=dfT, aes(x=runs, y=occurence)) +
  geom_point(colour="Black",size=1)+
  geom_point(data=dfF,colour="darkgrey",size=1) +
  geom_hline(yintercept=0.25,linetype="dashed",size=1) +
  scale_y_continuous(breaks=c(0,0.2,0.25,0.4,0.6),limits=c(0,0.7)) +
  # annotate(geom="text", x=1.2, y=0.7, label="(b)",size=6) +
  # ylim(0,0.7) +
  # ggtitle("Assignation frequencies of the best father\n",
  # subtitle=paste("n=1200",sep = "")) +
  xlab("") + ylab(" ") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=14,hjust = 0.5),
        axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=14),
        axis.text=element_text(size=14, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        #legend.key.height = unit(2,"line"),
        #legend.key.width = unit(5,"line"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14))
pFPic

#-----------------------------#
# Assignation frequencies of the best father & mother in Thu79
#-----------------------------#
distOcc=read.table("AssignationThu/outputs/occurencesBestFatherIncompleteGen.txt",
                   header = TRUE)
colnames(distOcc)=c("runs","occurence","father")
dfT=distOcc[which(distOcc$father==TRUE),1:2]
dfF=distOcc[which(distOcc$father==FALSE),1:2]
colnames(dfT)=c("runs","occurence")
colnames(dfF)=c("runs","occurence")
dfT$category = "True father"
dfF$category = "Wrong father"
df = rbind(dfT, dfF)
pMThu=ggplot(data=df, aes(x=runs, y=occurence, group = category, colour = category)) +
  geom_point(size=1)+
  # geom_point(data=dfF,colour="darkgrey",size=1) +
  geom_hline(yintercept=0.25,linetype="dashed",size=1) +
  scale_y_continuous(breaks=c(0,0.2,0.25,0.4,0.6),limits=c(0,0.7)) +
  # annotate(geom="text", x=1.2, y=0.7, label="(c)",size=6) +
  # ylim(0,0.7) +
  # ggtitle("Assignation frequencies of the best father\n",
  # subtitle=paste("n=1200",sep = "")) +
  xlab("Number of runs") + ylab("Assignment frequency") +
  scale_colour_manual(name='Assignment',
                     breaks=c('True father', 'Wrong father'),
                     values=c('True father'='black', 'Wrong father'='darkgrey'))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=14,hjust = 0.5),
        axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=14),
        axis.text=element_text(size=14, colour="black"),
        legend.key = element_rect(fill = "white", colour = "white", size = 3),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(2,"line"),
        legend.text=element_text(size=10),
        legend.title=element_text(size=10))
pMThu

distOcc=read.table("AssignationThu/outputs/occurencesBestMotherIncompleteGen.txt",
                   header = T)
colnames(distOcc)=c("runs","occurence","father")
dfT=distOcc[which(distOcc$father==TRUE),1:2]
dfF=distOcc[which(distOcc$father==FALSE),1:2]
colnames(dfT)=c("runs","occurence")
colnames(dfF)=c("runs","occurence")
pFThu=ggplot(data=dfT, aes(x=runs, y=occurence)) +
  geom_point(colour="Black",size=1)+
  geom_point(data=dfF,colour="darkgrey",size=1) +
  geom_hline(yintercept=0.25,linetype="dashed",size=1) +
  scale_y_continuous(breaks=c(0,0.2,0.25,0.4,0.6),limits=c(0,0.7)) +
  # annotate(geom="text", x=1.2, y=0.7, label="(d)",size=6) +
  # ylim(0,0.7) +
  # ggtitle("Assignation frequencies of the best father\n",
  # subtitle=paste("n=1200",sep = "")) +
  xlab("Number of runs") + ylab(" ") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=14,hjust = 0.5),
        axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=14),
        axis.text=element_text(size=14, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(2,"line"),
        legend.text=element_text(size=10),
        legend.title=element_text(size=10))
pFThu

# ggarrange(pMPic,pFPic,pMThu,pFThu,widths=1:1,heights=1:1,labels="AUTO")
ggarrange(pMPic,pMThu, widths = c(1, 1.5), labels="AUTO")

# Save for publication
ggsave("Figures/Fig1.jpeg",
       dpi=320,units="cm",width=23,height=8,
       create.dir = T)



#==========================================================#
# FIGURE 2. Offspring-father mating dispersal distances ----
#==========================================================#
# Figure 2. Distribution of the offspring-father mating dispersal distances in France (a) and Germany (b).
# Black dots are distances observed for accepted paternities (n = 82 offspring and 17 fathers in a and n = 62 offspring and 17 fathers in b).
# The grey ribbon is the confidence interval at 95% for expected distances under a random mating hypothesis
# between candidate fathers and true offspring, estimated with resampling methods.
# The black curve is the theoretical cumulative density function of distances of the selected model with FITDSITRPLUS,
# with parameters fitted on individual distances (dashed curves are confidence interval at 95%).
# The grey curve is the negative exponential cumulative density function of distances with the β parameter estimated by MASTERBAYES (dashed curves are confidence interval at 95%).

#-----------------------------#
# Distribution of offspring-father distances PicAll
#-----------------------------#
dyadsH0=read.table("AssignationPic/outputs/dyadsH0.txt", header=TRUE)
dyadsObsSelect=read.table("AssignationPic/outputs/dyadsObsSelect.txt", header=TRUE)
icH0=read.table("AssignationPic/outputs/icH0.txt", header=TRUE)

##### Fitted distribution are estimated in the Dispersal kernel directory
# Cumulative density function of the fitted Weibull distribution with fitdistrplus
paramPic = read.table("Tables/ParamPic.txt", header = T)
z=seq(0,100,by=0.01)
yfit=function( z ){pweibull(z, paramPic$mean[2], paramPic$mean[3])*(1-paramPic$mean[1]) + paramPic$mean[1]}
yfitInf=function( z ){pweibull(z, paramPic$lower[2], paramPic$lower[3])*(1-paramPic$upper[1]) + paramPic$upper[1]}
yfitSup=function( z ){pweibull(z, paramPic$upper[2], paramPic$upper[3])*(1-paramPic$lower[1]) + paramPic$lower[1]}

# Cumulative density function from MasterBayes estimates
lambda2=0.189627
lambda2Inf=0.3530
lambda2Sup=0.1280
ybayes=function( z ){1-exp(-1*lambda2*z)}
ybayesInf=function( z ){1-exp(-1*lambda2Inf*z)}
ybayesSup=function( z ){1-exp(-1*lambda2Sup*z)}


# DYADS RELATIVE FREQUENCY IN FUNCTION OF DISTANCE

df=data.frame(d=c(sort(dyadsH0$distance),sort(dyadsObsSelect$distance)),
              freq=c(seq(1,length(dyadsH0$distance),1)/length(dyadsH0$distance),seq(1,length(dyadsObsSelect$distance),1)/length(dyadsObsSelect$distance)),
              dist=c(rep("H0",length(dyadsH0$distance)),rep("ObsSelect",length(dyadsObsSelect$distance)))
)
df1=df[which(df$dist=="H0"),]
df3=df[which(df$dist=="ObsSelect"),]

pPic=ggplot(data=df1, aes(x=d, y=freq)) +
  stat_function(fun = yfit)+
  stat_function(fun = yfitInf,lty=2)+
  stat_function(fun = yfitSup,lty=2)+
  # stat_function(fun = yfitGamma,size=1.5,lty=3)+
  stat_function(fun = ybayes,colour="darkgrey")+
  stat_function(fun = ybayesInf,colour="darkgrey",lty=2)+
  stat_function(fun = ybayesSup,colour="darkgrey",lty=2)+
  # geom_point(aes(color="expected"),size=4)+
  geom_ribbon(data=icH0,aes(x=distance,y=meanFreq,ymin=lowerIC, ymax=upperIC),alpha=0.5) +
  # geom_point(data=df2,aes(color="observed"),size=4) +
  geom_point(data=df3,aes(color="observed with selection")) +
  # geom_line(data=icH0,aes(x=distance,y=lowerIC),color="grey",size=1.5) +
  # geom_line(data=icH0,aes(x=distance,y=upperIC),color="grey",size=1.5) +
  scale_colour_manual(values=c("black"))+
  # annotate(geom="text", x=0, y=1, label="(a)",size=6) +
  ylim(0,1) + xlim(0,100) +
  # scale_x_continuous(breaks=c(0,round(meanDispersal,digits=1),20,40,60),limits=c(0,60)) +
  # ggtitle(paste("Distribution of offspring-father distances\n(n=",nrow(df3),")"))+
  xlab("") + ylab("Dyads cumulative frequency") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=14),
        axis.text=element_text(size=14, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(5,"line"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14),
        legend.position='none')
pPic

#-----------------------------#
# Distribution of offspring-father distances ThuAll
#-----------------------------#
dyadsH0=read.table("AssignationThu/outputs/dyadsH0.txt",header=TRUE)
dyadsObsSelect=read.table("AssignationThu/outputs/dyadsObsSelect.txt",header=TRUE)
icH0=read.table("AssignationThu/outputs/icH0.txt",header=TRUE)

##### Fitted distribution are estimated in the Dispersal kernel directory
# Cumulative density function from fitdistrplus and Weibull distribution
paramThu = read.table("Tables/ParamThu.txt", header = T)
z=seq(0,100,by=0.01)
# yfit=function( z ){plnorm(z, paramThu$mean[2], paramThu$mean[3])*(1-paramThu$mean[1]) + paramThu$mean[1]}
# yfitInf=function( z ){plnorm(z, paramThu$lower[2], paramThu$lower[3])*(1-paramThu$upper[1]) + paramThu$upper[1]}
# yfitSup=function( z ){plnorm(z, paramThu$upper[2], paramThu$upper[3])*(1-paramThu$lower[1]) + paramThu$lower[1]}
yfit=function( z ){pweibull(z, paramThu$mean[2], paramThu$mean[3])*(1-paramThu$mean[1]) + paramThu$mean[1]}
yfitInf=function( z ){pweibull(z, paramThu$lower[2], paramThu$lower[3])*(1-paramThu$upper[1]) + paramThu$upper[1]}
yfitSup=function( z ){pweibull(z, paramThu$upper[2], paramThu$upper[3])*(1-paramThu$lower[1]) + paramThu$lower[1]}

# Cumulative density function from MasterBayes estimates
lambda2b=0.1293
lambda2Infb=0.111
lambda2Supb=0.1502
ybayesb=function( z ){1-exp(-1*lambda2b*z)}
ybayesInfb=function( z ){1-exp(-1*lambda2Infb*z)}
ybayesSupb=function( z ){1-exp(-1*lambda2Supb*z)}


# DYADS RELATIVE FREQUENCY IN FUNCTION OF DISTANCE

df=data.frame(d=c(sort(dyadsH0$distance),sort(dyadsObsSelect$distance)),
              freq=c(seq(1,length(dyadsH0$distance),1)/length(dyadsH0$distance),seq(1,length(dyadsObsSelect$distance),1)/length(dyadsObsSelect$distance)),
              dist=c(rep("H0",length(dyadsH0$distance)),rep("ObsSelect",length(dyadsObsSelect$distance)))
)
df1=df[which(df$dist=="H0"),]
df3=df[which(df$dist=="ObsSelect"),]

pThu=ggplot(data=df1, aes(x=d, y=freq)) +
  stat_function(fun = yfit)+
  stat_function(fun = yfitInf,lty=2)+
  stat_function(fun = yfitSup,lty=2)+
  # stat_function(fun = yfitGamma,size=1.5,lty=3)+
  stat_function(fun = ybayesb,colour="darkgrey")+
  stat_function(fun = ybayesInfb,colour="darkgrey",lty=2)+
  stat_function(fun = ybayesSupb,colour="darkgrey",lty=2)+
  # geom_point(aes(color="expected"),size=4)+
  geom_ribbon(data=icH0,aes(x=distance,y=meanFreq,ymin=lowerIC, ymax=upperIC),alpha=0.5) +
  # geom_point(data=df2,aes(color="observed"),size=4) +
  geom_point(data=df3,aes(color="observed with selection")) +
  # geom_line(data=icH0,aes(x=distance,y=lowerIC),color="grey",size=1.5) +
  # geom_line(data=icH0,aes(x=distance,y=upperIC),color="grey",size=1.5) +
  scale_colour_manual(values=c("black"))+
  # annotate(geom="text", x=0, y=1, label="(b)",size=6) +
  ylim(0,1) + xlim(0,100) +
  # scale_x_continuous(breaks=c(0,round(meanDispersal,digits=1),20,40,60),limits=c(0,60)) +
  # ggtitle(paste("Distribution of offspring-father distances\n(n=",nrow(df3),")"))+
  xlab("Geographic distance (km)") + ylab("Dyads cumulative frequency") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=14),
        axis.text=element_text(size=14, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(5,"line"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14),
        legend.position='none')
pThu

ggarrange(pPic,pThu,heights=1:1,align="v",nrow=2,labels="AUTO")

# Save in png for review and jpeg for publication
# ggsave("Figures/Fig2.png",
#        device="png",dpi=320,units="cm",width=19,height=26,
#        create.dir = T)
ggsave("Figures/Fig2.jpeg",
       dpi=320,units="cm",width=19,height=26,
       create.dir = T)


#==========================================================#
# FIGURE 3. Dispersal kernel ----
#==========================================================#
# Figure 3. Competing models of gametic dispersal kernels fitted on strictly positive dispersal distances in Picardy (a) and Thuringia (b).
# Histograms show densities of empirical gametic dispersal distances (i.e. natal + mating dispersal).
# In each case, the selected model is in solid line and rejected models are in dashed lines.
#-----------------------------#
# Dispersal kernel for Pic and Thu
#-----------------------------#
distPic = read.table("AssignationPic/outputs/distancesDistPic_gametic.txt", header = TRUE)
distThu = read.table("AssignationThu/outputs/distancesDistThu_gametic.txt", header = TRUE)

distPic = distPic[which(distPic$dist == "ObsSelect"),]
distThu = distThu[which(distThu$dist == "ObsSelect"),]

#-----------------------------#
# Distribution of offspring-father distances Picardy
#-----------------------------#
dyadsH0=read.table("AssignationPic/outputs/dyadsH0.txt",header=TRUE)
dyadsObsSelect=read.table("AssignationPic/outputs/dyadsObsSelect.txt",header=TRUE)
icH0=read.table("AssignationPic/outputs/icH0.txt",header=TRUE)

##### Fitted distribution are estimated in the Dispersal kernel directory
# Empirical probability density function
paramPicKernel = read.table("Tables/paramPicKernel GAMETIC.txt", header = T)
z=seq(0,100,by=0.01)

# Making functions for distribution kernel
library(gamlss.dist)
library(MASS)
library(fitdistrplus)
library(gnorm)
shape_lehnen = 1/16.77
expfit_lehnen = function( z ){dexp(z, shape_lehnen)}
weibullfit = function( z ){dweibull(z, paramPicKernel$param1[3], paramPicKernel$param2[3])}

linesize = 0.7

# Dispersal Kernel (PDF)
# df = data.frame(d = dyadsObsSelect$distance[dyadsObsSelect$distance > 0])
df = data.frame(d = distPic$d[distPic$d > 0])

histPic = ggplot(data=df, aes(x=d)) +
  geom_histogram(aes(y = ..density..), fill = "white", color = 'black', bins = 20) +
  stat_function(fun = weibullfit, color = "black", linetype = "solid", size = linesize) +
  stat_function(fun = expfit_lehnen, color = "black", linetype = "dashed", size = linesize) +
  # geom_density(color="black", size = 1) +
  ylim(0,0.07) +
  xlim(0,100) +
  # scale_x_continuous(breaks=c(0,round(meanDispersal,digits=1),20,40,60),limits=c(0,60)) +
  # ggtitle(paste("Distribution of offspring-father distances\n(n=",nrow(df3),")"))+
  xlab("Geographic distance (km)") + ylab("Density") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=14),
        axis.text=element_text(size=14, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(2,"line"),
        legend.text=element_text(size=10),
        legend.title=element_text(size=10, face = "bold"))
histPic

#-----------------------------#
# Distribution of offspring-father distances ThuAll
#-----------------------------#
dyadsH0=read.table("AssignationThu/outputs/dyadsH0.txt",header=TRUE)
dyadsObsSelect=read.table("AssignationThu/outputs/dyadsObsSelect.txt",header=TRUE)
icH0=read.table("AssignationThu/outputs/icH0.txt",header=TRUE)

##### Fitted distribution are estimated in the Dispersal kernel directory
# Empirical probability density function
paramThuKernel = read.table("Tables/ParamThuKernel GAMETIC.txt", header = T)
z=seq(0,100,by=0.01)

# Making functions for distribution kernel
shape_lehnen = 1/16.77
expfit_lehnen = function( z ){dexp(z, shape_lehnen)}
weibullfit = function( z ){dweibull(z, paramThuKernel$param1[3], paramThuKernel$param2[3])}

# Dispersal Kernel (PDF)
# df = data.frame(d=dyadsObsSelect$distance[dyadsObsSelect$distance > 0])
df = data.frame(d=distThu$d[distThu$d > 0])

histThu = ggplot(data=df, aes(x=d)) +
  geom_histogram(aes(y = ..density..), fill = "white", color = 'black', bins = 20) +
  stat_function(fun = weibullfit, color = "black", linetype = "solid", size = linesize) +
  stat_function(fun = expfit_lehnen, color = "black", linetype = "dashed", size = linesize) +
  # geom_density(color="black", size = 1) +
  ylim(0,0.07) +
  xlim(0,100) +
  xlab("Geographic distance (km)") + ylab("Density") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=14),
        axis.text=element_text(size=14, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(2,"line"),
        legend.text=element_text(size=10),
        legend.title=element_text(size=10, face = "bold"))

histThu

ggarrange(histPic,histThu,heights=1:1,align="v",nrow=1, ncol =2,labels="AUTO", common.legend = TRUE, legend = "right")

ggsave("Figures/Fig3.png",
       device="png",dpi=320,units="cm",width=26,height=10,
       create.dir = T)
ggsave("Figures/Fig3.jpeg",
       dpi=320,units="cm",width=26,height=10,
       create.dir = T)





#==========================================================#
# FIGURE S1. Dispersal kernel ----
#==========================================================#
# Figure S1. Competing models of gametic dispersal kernels fitted on strictly positive dispersal distances in Picardy (a) and Thuringia (b).
# Histograms show densities of empirical gametic dispersal distances (i.e. natal + mating dispersal).
# In each case, the selected model is in solid line and rejected models are in dashed lines.
#-----------------------------#
# Dispersal kernel for Pic and Thu
#-----------------------------#
distPic = read.table("AssignationPic/outputs/distancesDistPic_gametic.txt", header = TRUE)
distThu = read.table("AssignationThu/outputs/distancesDistThu_gametic.txt", header = TRUE)

distPic = distPic[which(distPic$dist == "ObsSelect"),]
distThu = distThu[which(distThu$dist == "ObsSelect"),]

#-----------------------------#
# Distribution of offspring-father distances Picardy
#-----------------------------#
dyadsH0=read.table("AssignationPic/outputs/dyadsH0.txt",header=TRUE)
dyadsObsSelect=read.table("AssignationPic/outputs/dyadsObsSelect.txt",header=TRUE)
icH0=read.table("AssignationPic/outputs/icH0.txt",header=TRUE)

##### Fitted distribution are estimated in the Dispersal kernel directory
# Empirical probability density function
paramPicKernel = read.table("Tables/paramPicKernel GAMETIC.txt", header = T)
z=seq(0,100,by=0.01)

# Making functions for distribution kernel
library(gamlss.dist)
library(MASS)
library(fitdistrplus)
library(gnorm)
expfit = function( z ){dexp(z, paramPicKernel$param1[1])}
normfit = function( z ){dnorm(z, paramPicKernel$param1[2], paramPicKernel$param2[2])}
weibullfit = function( z ){dweibull(z, paramPicKernel$param1[3], paramPicKernel$param2[3])}
gammafit = function( z ){dgamma(z, paramPicKernel$param1[4], paramPicKernel$param2[4])}
lnormfit = function( z ){dlnorm(z, paramPicKernel$param1[5], paramPicKernel$param2[5])}
expowerfit = function( z ){dgnorm(z, paramPicKernel$param1[6], paramPicKernel$param2[6], paramPicKernel$param3[6])}

library(RColorBrewer)
palette = brewer.pal(9, "Set1")

linesize = 0.7

# Dispersal Kernel (PDF)
# df = data.frame(d = dyadsObsSelect$distance[dyadsObsSelect$distance > 0])
df = data.frame(d = distPic$d[distPic$d > 0])

histPic = ggplot(data=df, aes(x=d)) +
  geom_histogram(aes(y = ..density..), fill = "white", color = 'black', bins = 20) +
  stat_function(fun = expfit,aes(color = "a"), linetype = "dashed", size = linesize) +
  stat_function(fun = normfit, aes(color = "b"), linetype = "dashed", size = linesize) +
  stat_function(fun = gammafit, aes(color = "d"), linetype = "dashed", size = linesize) +
  stat_function(fun = lnormfit, aes(color = "e"), linetype = "dashed", size = linesize) +
  stat_function(fun = expowerfit, aes(color = "f"), linetype = "dashed", size = linesize) +
  stat_function(fun = weibullfit, aes(color = "c"), linetype = "solid", size = linesize) +
  # geom_density(color="black", size = 1) +
  ylim(0,0.07) +
  xlim(0,100) +
  scale_colour_manual(name = 'Distribution', 
                      values =c('a'=palette[1],'b'=palette[2],'c'=palette[3],'d'=palette[4],
                                'e'=palette[5],'f'=palette[7]),
                      labels = c('Exponential','Gaussian','Weibull','Gamma','Log-normal','Exponential power')) +
  # scale_x_continuous(breaks=c(0,round(meanDispersal,digits=1),20,40,60),limits=c(0,60)) +
  # ggtitle(paste("Distribution of offspring-father distances\n(n=",nrow(df3),")"))+
  xlab("Geographic distance (km)") + ylab("Density") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=14),
        axis.text=element_text(size=14, colour="black"),
        legend.key = element_rect(fill = "white", colour = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(2,"line"),
        legend.text=element_text(size=10),
        legend.title=element_text(size=10, face = "bold"))
histPic

#-----------------------------#
# Distribution of offspring-father distances ThuAll
#-----------------------------#
dyadsH0=read.table("AssignationThu/outputs/dyadsH0.txt",header=TRUE)
dyadsObsSelect=read.table("AssignationThu/outputs/dyadsObsSelect.txt",header=TRUE)
icH0=read.table("AssignationThu/outputs/icH0.txt",header=TRUE)

##### Fitted distribution are estimated in the Dispersal kernel directory
# Empirical probability density function
paramThuKernel = read.table("Tables/ParamThuKernel GAMETIC.txt", header = T)
z=seq(0,100,by=0.01)

# Making functions for distribution kernel
expfit = function( z ){dexp(z, paramThuKernel$param1[1])}
normfit = function( z ){dnorm(z, paramThuKernel$param1[2], paramThuKernel$param2[2])}
weibullfit = function( z ){dweibull(z, paramThuKernel$param1[3], paramThuKernel$param2[3])}
gammafit = function( z ){dgamma(z, paramThuKernel$param1[4], paramThuKernel$param2[4])}
lnormfit = function( z ){dlnorm(z, paramThuKernel$param1[5], paramThuKernel$param2[5])}
expowerfit = function( z ){dgnorm(z, paramThuKernel$param1[6], paramThuKernel$param2[6], paramThuKernel$param3[6])}


# Dispersal Kernel (PDF)
# df = data.frame(d=dyadsObsSelect$distance[dyadsObsSelect$distance > 0])
df = data.frame(d=distThu$d[distThu$d > 0])

histThu = ggplot(data=df, aes(x=d)) +
  geom_histogram(aes(y = ..density..), fill = "white", color = 'black', bins = 20) +
  stat_function(fun = expfit,aes(color = "a"), linetype = "dashed", size = linesize) +
  stat_function(fun = normfit, aes(color = "b"), linetype = "dashed", size = linesize) +
  stat_function(fun = weibullfit, aes(color = "c"), linetype = "solid", size = linesize) +
  stat_function(fun = gammafit, aes(color = "d"), linetype = "dashed", size = linesize) +
  stat_function(fun = expowerfit, aes(color = "f"), linetype = "dashed", size = linesize) +
  stat_function(fun = lnormfit, aes(color = "e"), linetype = "dashed", size = linesize) +
  # geom_density(color="black", size = 1) +
  ylim(0,0.07) +
  xlim(0,100) +
  scale_colour_manual(name = 'Distribution', 
                      values =c('a'=palette[1],'b'=palette[2],'c'=palette[3],'d'=palette[4],
                                'e'=palette[5],'f'=palette[7]),
                      labels = c('Exponential','Gaussian','Weibull','Gamma','Log-normal','Exponential power')) +
  # scale_x_continuous(breaks=c(0,round(meanDispersal,digits=1),20,40,60),limits=c(0,60)) +
  # ggtitle(paste("Distribution of offspring-father distances\n(n=",nrow(df3),")"))+
  xlab("Geographic distance (km)") + ylab("Density") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=14),
        axis.text=element_text(size=14, colour="black"),
        legend.key = element_rect(fill = "white", colour = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(2,"line"),
        legend.text=element_text(size=10),
        legend.title=element_text(size=10, face = "bold"))

histThu

ggarrange(histPic,histThu,heights=1:1,align="v",nrow=1, ncol =2,labels="AUTO", common.legend = TRUE, legend = "right")

ggsave("Figures/FigS1.png",
       device="png",dpi=320,units="cm",width=26,height=10,
       create.dir = T)
ggsave("Figures/FigS1.jpeg",
       dpi=320,units="cm",width=26,height=10,
       create.dir = T)




#==========================================================#
# FIGURE S2. IBD ----
#==========================================================#

#==========================================================#
# Picardy ----
#==========================================================#

#----------------------------------------------------------#
# MATRIX OF GENETIC DISTANCES

allel=read.table("Data/Pic/uniqueGenotypesWithInfo.txt",h=T)
allel$idcol = as.factor(allel$idcol)

# Remove some colonies (too far away)
allel=allel[-which(allel$idcol=="M399" | allel$idcol=="M1079" | allel$idcol=="CXSGT"| allel$idcol=="M1975"  | allel$idcol=="M1979"),]
# allel=allel[allel$idcol != "M399",]

# Keep females only
allel=allel[-which(allel$sexe=="M"),]
# Remove juveniles
allel=allel[-which(allel$ageWhenFirstCaptur=="Juv"),]

# Remove individuals with uncomplete genotypes
allel=allel[!apply(allel[,2:17],1,function (x) 0 %in% x),]

colonies = as.character(allel$idcol)

genotypes = data.frame(rha101 = paste(allel[,2], allel[,3], sep = "/"),
                       rha109 = paste(allel[,4], allel[,5], sep = "/"),
                       rha4 = paste(allel[,6], allel[,7], sep = "/"),
                       rha7 = paste(allel[,8], allel[,9], sep = "/"),
                       rhc108 = paste(allel[,10], allel[,11], sep = "/"),
                       rhc3 = paste(allel[,12], allel[,13], sep = "/"),
                       rhd102 = paste(allel[,14], allel[,15], sep = "/"),
                       rhd103 = paste(allel[,16], allel[,17], sep = "/"))



matAlleles = df2genind(genotypes,
                       sep="/",
                       ploidy=2,
                       type="codom",
                       pop=colonies)
matAlleles
ploidy(matAlleles)

# Pairwise Fst between colonies
# matGen=adegenet::pairwise.fst(matAlleles,res.type="matrix")
# Since a change in adegenet, the function moved to hierfstat
matGen = genet.dist(matAlleles, method = "Nei87")
# matGen = pairwise.neifst(matAlleles)
matGen

# Genetic distance following Rousset (1997)
matGen=matGen/(1-matGen)
matGen

#-------------------------------#
# MATRIX OF GEOGRAPHIC DISTANCES

# Distance geographique a recalculer
# Calcul d'une matrice de distances
coordCol=read.table("Data/Pic/coordPic.txt",h=T)
coordCol=coordCol[-which(coordCol$Colonie=="M399" | coordCol$Colonie=="M1079" | coordCol$Colonie=="CXSGT"| coordCol$Colonie=="M1975"  | coordCol$Colonie=="M1979" | coordCol$Colonie=="P412"),]
# coordCol=coordCol[-which(coordCol$Colonie=="M399" | coordCol$Colonie=="P412"),]
# coordCol=coordCol[coordCol$Colonie != "M399",]

lvl=unique(colonies)
# COlonies in the same order as in the genetic matrix
coord=data.frame(Colonie=rep(NA,length(lvl)),
                 Long=rep(NA,length(lvl)),
                 Lat=rep(NA,length(lvl)))
coord$Colonie=lvl
for (i in 1:nrow(coord)) {
  coord$Long[i]=coordCol$Long[which(coordCol$Colonie==coord$Colonie[i])]
  coord$Lat[i]=coordCol$Lat[which(coordCol$Colonie==coord$Colonie[i])]
}

# Calculating the matrix
matGeo=matrix(nrow=nrow(coord),ncol=nrow(coord))
colnames(matGeo)=coord$Colonie
rownames(matGeo)=coord$Colonie

# create distance matrix
matGeo=distm(cbind(coord$Long,coord$Lat), fun=distCosine)
# Matrix in km
matGeo=matGeo/1000
# log distance as in Rousset (1997)
matGeo=log(matGeo)
# Correcting negative values due to log
matGeo[which(matGeo<0)]=0


#-------------------------------#
# MANTEL TEST

# "For this purpose, we first quantified the linear relationship between genetic distance FST/(1-FST)
# and the logarithm of Euclidean geographic distance (Rousset 1997), i.e. isolation-by-distance (IBD)
# and tested its significance with a Mantel test using 10000 permutations" (Jan 2017)

matGen[matGen == 0] = NA
mantel.rtest(as.dist(matGeo,diag=FALSE,upper=FALSE),
             as.dist(matGen,diag=FALSE,upper=FALSE),
             nrepet=1000)

# PLOT IBD PATTERN
# df=data.frame(Distance=as.vector(matGeo[lower.tri(matGeo,diag=FALSE)]),GeneDistance=as.vector(matGen[lower.tri(matGen,diag=FALSE)]))
df=data.frame(Distance=as.vector(matGeo[lower.tri(matGeo,diag=FALSE)]),GeneDistance=as.vector(matGen))
#df=data.frame(Distance=as.vector(Dgeo),GeneDistance=as.vector(Dgen))
plot(df)
mod=lm(GeneDistance~Distance,data=df)
abline(lm(GeneDistance~Distance,data=df))
mod$coefficients[1]
mod$coefficients[2]

#-------------------------------#
# GRAPH GGPLOT OF IBD PATTERN (13 colonies)

p1 = ggplot(data=df, aes(x=Distance, y=GeneDistance)) +
  geom_point(size=2)+
  geom_abline(intercept=mod$coefficients[1],slope=mod$coefficients[2],size=2) +
  xlab("Geographic distance (km)") + ylab("Genetic distance Fst/(1-Fst)") +
  scale_x_continuous(breaks=log(c(2,5,10,25)), labels=c(2,5,10,25)) +
  # scale_y_continuous(breaks=c(0, 0.006, 0.012), labels=c(0, 0.006, 0.012)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=14),
        axis.text=element_text(size=14, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(2,"line"),
        legend.text=element_text(size=10),
        legend.title=element_text(size=10, face = "bold"))
p1


#-------------------------------#
# Mean dispersal distance d'apres IBD pattern
#-------------------------------#

# Estimate of effective density
# 1000 sampled females ; 3382 sampled genotypes in 17 colonies, 2851 sampled gen in 13 colonies
# half of colonies were sampled (personal communication, Eric Petit)
# so 1000*2(half sampled)*2(for males)=4000 bats in the area (theoretical) or 5702/6764 (computed from data set)
# n effective = 10% * total = 400 or 570/675
# area of the sampling site = 500km2 (forêt) ou 2500km2 (conservatif, comprend tous les points de la zone)
# 
# d=c(400/500,400/2500) # Effective density (ind.km-2)
# p=c(0.005826133,0.01095223) # slope of IBD regression (for the subset of 13 colonies or 17 colonies)
mod$coefficients

d = 400/500
p = mod$coefficients[2]
# Formula the mean dispersal estimated from slope of IBD regression, given by Rousset (1997)
mean.dispersal = sqrt(1/(4*d*pi*p))
mean.dispersal


#==========================================================#
# Thuringia ----
#==========================================================#

#----------------------------------------------------------#
# MATRIX OF GENETIC DISTANCES

allel=read.table("Data/Thu/uniqueGenotypesWithInfo.txt",h=T)
allel$idcol = as.factor(allel$idcol)


allel = allel[allel$idcol != "Thu42",]

# enlever les males (femelles sont liees a la colonie beaucoup plus que les males)
allel=allel[-which(allel$sexe=="M"),]
# enlever les juveniles
allel=allel[-which(allel$ageWhenFirstCaptur=="Juv"),]
# Supprimer les individus avec un genotype incomplet
# i.e. individus avec un 0 parmi tous les alleles
#Enlever les individus avec missing data
allel=allel[!apply(allel[,2:17],1,function (x) 0 %in% x),]
# Construction de l'objet genind qui sera transmis pour la Fst
colonies = as.character(allel$idcol)

genotypes = data.frame(rha101 = paste(allel[,2], allel[,3], sep = "/"),
                       rha109 = paste(allel[,4], allel[,5], sep = "/"),
                       rha4 = paste(allel[,6], allel[,7], sep = "/"),
                       rha7 = paste(allel[,8], allel[,9], sep = "/"),
                       rhc108 = paste(allel[,10], allel[,11], sep = "/"),
                       rhc3 = paste(allel[,12], allel[,13], sep = "/"),
                       rhd102 = paste(allel[,14], allel[,15], sep = "/"),
                       rhd103 = paste(allel[,16], allel[,17], sep = "/"))



matAlleles = df2genind(genotypes,
                       sep="/",
                       ploidy=2,
                       type="codom",
                       pop=colonies)
matAlleles
ploidy(matAlleles)

matGen = genet.dist(matAlleles, method = "Nei87")


# Distance genetique selon Rousset (1997)
matGen=matGen/(1-matGen)

#-------------------------------#
# MATRICE DES DISTANCES GEOGRAPHIQUES

# Distance geographique a recalculer
# Calcul d'une matrice de distances
coordCol=read.table("Data/Thu/coordThu.txt",h=T)

coordCol = coordCol[coordCol$Colony != "Thu42",]

lvl=unique(colonies)
# Creation du fichier contenant les coordonnees dans le bon ordre
# !!!! Meme ordre que dans la matrice des distances genetiques matGen
coord=data.frame(Colonie=rep(NA,length(lvl)),
                 Long=rep(NA,length(lvl)),
                 Lat=rep(NA,length(lvl)))
coord$Colonie=lvl
for (i in 1:nrow(coord)) {
  coord$Long[i]=coordCol$Long[which(coordCol$Colony==coord$Colonie[i])]
  coord$Lat[i]=coordCol$Lat[which(coordCol$Colony==coord$Colonie[i])]
}

matGeo=matrix(nrow=nrow(coord),ncol=nrow(coord)) # Matrix de 17*17 -> 17 colonies prises en compte
colnames(matGeo)=coord$Colonie
rownames(matGeo)=coord$Colonie
# Calcul des distances paire a paire

# create distance matrix
matGeo=distm(cbind(coord$Long,coord$Lat), fun=distCosine)
# Matrix in km
matGeo=matGeo/1000
# Si on veut les distances en log selon Rousset (1997)
matGeo=log(matGeo)
# Correction des valeurs negatives apres le passage en log
matGeo[which(matGeo<0)]=0



#-------------------------------#
# MANTEL TEST

matGen[matGen == 0] = NA
mantel.rtest(as.dist(matGeo,diag=FALSE,upper=FALSE),
             as.dist(matGen,diag=FALSE,upper=FALSE),
             nrepet=1000)


# PLOT IBD PATTERN
# df=data.frame(Distance=as.vector(matGeo[lower.tri(matGeo,diag=FALSE)]),GeneDistance=as.vector(matGen[lower.tri(matGen,diag=FALSE)]))
df=data.frame(Distance=as.vector(matGeo[lower.tri(matGeo,diag=FALSE)]),GeneDistance=as.vector(matGen))

##### KEEPING ONLY DISTANCES INFERIOR TO 30KM (SHORT SCALE)
df = df[which(df$Distance < log(30)),]

#df=data.frame(Distance=as.vector(Dgeo),GeneDistance=as.vector(Dgen))
plot(df)
mod=lm(GeneDistance~Distance,data=df)
abline(lm(GeneDistance~Distance,data=df))
mod$coefficients[1]
mod$coefficients[2]

# save(df,file="Data/Pic/IBDPic.Rdata")
# save(mod,file="Data/Pic/modIBDPic.Rdata")

#-------------------------------#
# GRAPH GGPLOT OF IBD PATTERN

# export 1800x1300
p2 = ggplot(data=df, aes(x=Distance, y=GeneDistance)) +
  geom_point(size=2)+
  geom_abline(intercept=mod$coefficients[1],slope=mod$coefficients[2],size=2) +
  xlab("Geographic distance (km)") + ylab("") +
  scale_x_continuous(breaks=log(c(2,5,10,25)), labels=c(2,5,10,25)) +
  # scale_y_continuous(breaks=c(0, 0.02, 0.04, 0.06), labels=c(0, 0.02, 0.04, 0.06)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=14),
        axis.text=element_text(size=14, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(2,"line"),
        legend.text=element_text(size=10),
        legend.title=element_text(size=10, face = "bold"))
p2

#-------------------------------#
# Mean dispersal distance d'apres IBD pattern
#-------------------------------#

# Estimate of effective density
# 1000 sampled females ; 3382 sampled genotypes in 17 colonies, 2851 sampled gen in 13 colonies
# half of colonies were sampled (personal communication, Eric Petit)
# so 1000*2(half sampled)*2(for males)=4000 bats in the area (theoretical) or 5702/6764 (computed from data set)
# n effective = 10% * total = 400 or 570/675
# area of the sampling site = 500km2 (forêt) ou 2500km2 (conservatif, comprend tous les points de la zone)

mod$coefficients

d = 400/500
p = mod$coefficients[2]
# Formula the mean dispersal estimated from slope of IBD regression, given by Rousset (1997)
mean.dispersal=sqrt(1/(4*d*pi*p))
mean.dispersal




p = ggarrange(p1, p2, ncol = 2)
p

ggsave("Figures/FigS2.jpeg",
       dpi=320,units="cm",width=26,height=10,
       create.dir = T)



#==========================================================#
# END
#==========================================================#
