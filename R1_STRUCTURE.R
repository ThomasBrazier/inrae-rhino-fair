#==========================================================#
# LOADING ENVIRONMENT ----
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
# source("Sources/packages.R")
library(tidyverse)
library(ggplot2)

#=================================================#
# Import data ----
#=================================================#
dataPic = read.table("Data/Pic/uniqueGenotypesWithInfo.txt", header = TRUE,
                     stringsAsFactors = T)
dataThu = read.table("Data/Thu/uniqueGenotypesWithInfo.txt", header = TRUE,
                     stringsAsFactors = T)

resPic = read.table("Data/Pic/resultsWithInfo.txt", header = TRUE,
                    stringsAsFactors = T)
resThu = read.table("Data/Thu/resultsWithInfo.txt", header = TRUE,
                    stringsAsFactors = T)



# TODO Check the distribution of proba first/second
# TODO Describe natal dispersal more in details
# - Fix number of natal dispersers (with min-max range/interval estimated across 10 MCMC chains)
# - distribution of natal dispersal distances alone (Supp Fig)
# - output and distances of STRUCTURE (Supp Table)

# TODO Integrate uncertainty across proba of the best MCMC chain
# TODO Barplot of natal dispersal rate, for each rank of proba: first, second, third... up to six
# See how informative AND conservative (or not) we are
# TODO Dist of mean proba per rank: 1, 2, 3...

# DONE Integrate uncertainty with model averaging across 10 MCMC chains:
# - Average across MCMC chains, kernel with:
# Min distance
# Max distance
# Mean distance
# - Average of five best proba? NO


# TODO Check proba of assignment and prior error rates of known migrants:
# % of individuals born in one colony and sampled the year after in another. Birth colony already known for sampled fathers?
# The average natal dispersal rate vs the pop assignment-based natal dispersal rate


# TODO Prior error rate computation -> K-fold cross validation, estimate colony of individuals whom birth colony are already known
# DONE Consistency of proba of assignment across MCMC chains


#=================================================#
# Investigate the distributions of Proba of Assignment by STRUCTURE ----
#=================================================#



# Check if Father colony already known
# Age when first captured
# If Juvenile, then birth colony is known
table(resPic$ageFather)
table(resThu$ageFather)

# Only Adults, no known colony


#------------------------------------------------------------------------#
# PICARDY ----
#------------------------------------------------------------------------#
# Retrieve lines in the output file concerning only males involved in paternities
# idind = 2960 # Single  individual for testing
i = 2 # Select the run to estimate ancestry probabilities
individuals = unique(resPic$father) # List of individuals to retrieve (inferred fathers)

# Create the individuals infos
# Assign to pre- or post-dispersal
# pre- is for learning (adult females)
# post- is for assignment (adult males)
samplePic = dataPic[which(dataPic$ageWhenFirstCaptur == "Adult"),]
dfPic = data.frame(indivID = samplePic$idind,
                   PopID = as.numeric(samplePic$idcol),
                   PopFlag = as.numeric(samplePic$sexe == "F")
)
# Colonies are encoded in integers
# Translate from integers to names
colonynamesPic = levels(samplePic$idcol)[as.numeric(samplePic$idcol)]
# Add genotypes
dfPic = cbind(dfPic, samplePic[,2:17])



# Run a for loop on the list of individuals
# Re-init result file
individuals_inferred = dataPic$idind[which(dataPic$ageWhenFirstCaptur == "Adult" & dataPic$sexe == "M")]
# Individuals used for training (females)
individuals_known = dataPic$idind[which(dataPic$ageWhenFirstCaptur == "Adult" & dataPic$sexe == "F")]

res_Pic_inferred = read.table("STRUCTURE/output/Pic_R1/results_inferred.txt", header = F)
res_Pic_inferred = res_Pic_inferred[,-c(1, 3, 5)] # Remove unwanted columns
colnames(res_Pic_inferred) = c("idind", "idcol", levels(as.factor(colonynamesPic))) # Rename columns
res_Pic_inferred$idcol = (levels(as.factor(colonynamesPic))[res_Pic_inferred$idcol])
res_Pic_inferred = res_Pic_inferred[,1:19]

res_Pic_inferred$ancestry = NA
res_Pic_inferred$probancestry = NA
res_Pic_inferred$probancestry_second = NA

for (i in 1:nrow(res_Pic_inferred)) {
  res_Pic_inferred$ancestry[i] = as.character(colnames(res_Pic_inferred[,3:19])[which.max(res_Pic_inferred[i,3:19])])
  res_Pic_inferred$probancestry[i] = sort(as.numeric(res_Pic_inferred[i,3:19]), decreasing = T)[1]
  res_Pic_inferred$probancestry_second[i] = sort(as.numeric(res_Pic_inferred[i,3:19]), decreasing = T)[2]
}

(res_Pic_inferred$delta_proba = res_Pic_inferred$probancestry - res_Pic_inferred$probancestry_second)

res_Pic_inferred$nataldisp = (res_Pic_inferred$idcol != res_Pic_inferred$ancestry)
table(res_Pic_inferred$nataldisp)
# 50/50 chances of migrants/non-migrants
mean(res_Pic_inferred$probancestry)

# Are probabilities of ancestry for migrant or non-migrant the same?
mean(res_Pic_inferred$probancestry[res_Pic_inferred$nataldisp == TRUE])
mean(res_Pic_inferred$probancestry[res_Pic_inferred$nataldisp == FALSE])

mean(res_Pic_inferred$delta_proba[res_Pic_inferred$nataldisp == TRUE])
mean(res_Pic_inferred$delta_proba[res_Pic_inferred$nataldisp == FALSE])

hist(res_Pic_inferred$probancestry[res_Pic_inferred$nataldisp == TRUE], breaks = 20, xlim = c(0, 0.7))
hist(res_Pic_inferred$probancestry[res_Pic_inferred$nataldisp == FALSE], breaks = 20, xlim = c(0, 0.7))

hist(res_Pic_inferred$delta_proba[res_Pic_inferred$nataldisp == TRUE], breaks = 20, xlim = c(0, 0.7))
hist(res_Pic_inferred$delta_proba[res_Pic_inferred$nataldisp == FALSE], breaks = 20, xlim = c(0, 0.7))







PopAssignPic = read.table("STRUCTURE/output/Pic/results.txt", header = FALSE)
PopAssignPic = PopAssignPic[,-c(1, 3, 5)] # Remove unwanted columns
colnames(PopAssignPic) = c("idind", "idcol", levels(as.factor(colonynamesPic))) # Rename columns
PopAssignPic$idcol = (levels(as.factor(colonynamesPic))[PopAssignPic$idcol])
PopAssignPic = PopAssignPic[,1:19]

# Now which is the inferred ancestry and at which probability
PopAssignPic$ancestry = NA
PopAssignPic$probancestry = NA
PopAssignPic$probancestry_second = NA

for (i in 1:nrow(PopAssignPic)) {
  PopAssignPic$ancestry[i] = as.character(colnames(PopAssignPic[,3:19])[which.max(PopAssignPic[i,3:19])])
  PopAssignPic$probancestry[i] = sort(as.numeric(PopAssignPic[i,3:19]), decreasing = T)[1]
  PopAssignPic$probancestry_second[i] = sort(as.numeric(PopAssignPic[i,3:19]), decreasing = T)[2]
}

# Assess the validity of the inference
# Is the probability of the second probable ancestry much lower than the first one?
(PopAssignPic$delta_proba = PopAssignPic$probancestry - PopAssignPic$probancestry_second)

# Is there natal dispersal?
# i.e. ancestry and colony are not the same
# TRUE are migrants
# FALSE are non-migrants
PopAssignPic$nataldisp = (PopAssignPic$idcol != PopAssignPic$ancestry)
table(PopAssignPic$nataldisp)
# 50/50 chances of migrants/non-migrants
mean(PopAssignPic$probancestry)

# Are probabilities of ancestry for migrant or non-migrant the same?
mean(PopAssignPic$probancestry[PopAssignPic$nataldisp == TRUE])
mean(PopAssignPic$probancestry[PopAssignPic$nataldisp == FALSE])

mean(PopAssignPic$delta_proba[PopAssignPic$nataldisp == TRUE])
mean(PopAssignPic$delta_proba[PopAssignPic$nataldisp == FALSE])

hist(PopAssignPic$probancestry[PopAssignPic$nataldisp == TRUE], breaks = 20, xlim = c(0, 0.7))
hist(PopAssignPic$probancestry[PopAssignPic$nataldisp == FALSE], breaks = 20, xlim = c(0, 0.7))

hist(PopAssignPic$delta_proba[PopAssignPic$nataldisp == TRUE], breaks = 20, xlim = c(0, 0.7))
hist(PopAssignPic$delta_proba[PopAssignPic$nataldisp == FALSE], breaks = 20, xlim = c(0, 0.7))


# Test it with Mann-Whitney
wilcox.test(PopAssignPic$probancestry[PopAssignPic$nataldisp == TRUE], PopAssignPic$probancestry[PopAssignPic$nataldisp == FALSE])
# Slightly significant for individuals of interest

# Export results
dir.create("Tables")
# write.table(PopAssignPic, "Tables/MigrantFathersPic.txt", col.names = TRUE, row.names = FALSE,
#             quote = FALSE)
# PopAssignPic = read.table("Tables_supp/MigrantFathersPic.txt", header = TRUE)




# Yang et al. (2005): At least 15-20 markers to achieve a confident assignement. Huge loss of accuracy under 15 markers.
# If the most informative markers are selected, then we can reach >99% assignment with only 4-6 markers.
# Besides, STRUCTURE need a sufficient sample size to detect informations
# Broquet et al. 2009
# "The bias and precision of the results from Structure and BayesAss is
# nearly indistinguishable from our method (Fig. 8). Our
# method gave slightly more accurate results than Structure, while BayesAss returned the least-biased estimates
# (especially for philopatry). However, for BayesAss, the
# confidence intervals for the estimates with the largest error
# were often too small, which meant that they did not
# include the true value (this occurred in 17% of immigration
# estimates and 58% of philopatry estimates, as compared to
# 13% and 14%, respectively, with our approach)."
# STRUCTURE. Underestimation of philopatry.

# # Explore results for all individuals assigned (PopFlag == 0)
# i = 1 # Select the run to estimate ancestry probabilities
# individuals = unique(dfPic$indivID[dfPic$PopFlag == 0]) # List of individuals to retrieve (inferred fathers)
# 
# # Run a for loop on the list of individuals
# # Re-init result file
# system("rm STRUCTURE/output/Pic/resultsAll.txt")
# for (idind in individuals) {
#   # Copy the lines for individuals of interest in a result file
#   system(paste("grep ' ", idind, "    (' STRUCTURE/output/Pic/run", i, "_f >> STRUCTURE/output/Pic/resultsAll.txt", sep = ""))
# }
# # system("cat STRUCTURE/output/Pic/resultsAll.txt")
# 
# PopAssign = read.table("STRUCTURE/output/Pic/resultsAll.txt", header = FALSE)
# PopAssign = PopAssign[,-c(1, 3, 5)] # Remove unwanted columns
# colnames(PopAssign) = c("idind", "idcol", levels(as.factor(colonynamesPic))) # Rename columns
# PopAssign$idcol = (levels(as.factor(colonynamesPic))[PopAssign$idcol])
# PopAssign = PopAssign[,1:19]
# 
# # Now which is the inferred ancestry and at which probability
# PopAssign$ancestry = NA
# PopAssign$probancestry = NA
# for (i in 1:nrow(PopAssign)) {
#   PopAssign$ancestry[i] = as.character(colnames(PopAssign[,3:19])[which.max(PopAssign[i,3:19])])
#   PopAssign$probancestry[i] = max(PopAssign[i,3:19]) 
# }
# # Assess the validity of the inference
# # Is the probability of the second probable ancestry much lower than the first one?
# 
# # Is there natal dispersal?
# # i.e. ancestry and colony are the same
# PopAssign$nataldisp = (PopAssign$idcol != PopAssign$ancestry)
# table(PopAssign$nataldisp)
# # Proportion of migrants
# sum(PopAssign$nataldisp == TRUE)/nrow(PopAssign)
# # 68% of migrants seems a high percentage
# 
# # Are probabilities of ancestry for migrant or non-migrant the same?
# mean(PopAssign$probancestry[PopAssign$nataldisp == TRUE])
# mean(PopAssign$probancestry[PopAssign$nataldisp == FALSE])
# hist(PopAssign$probancestry[PopAssign$nataldisp == TRUE], breaks = 20, xlim = c(0, 0.7))
# hist(PopAssign$probancestry[PopAssign$nataldisp == FALSE], breaks = 20, xlim = c(0, 0.7))
# # Test it with Mann-Whitney
# wilcox.test(PopAssign$probancestry[PopAssign$nataldisp == TRUE], PopAssign$probancestry[PopAssign$nataldisp == FALSE])
# # Not significant for individuals of interest




#------------------------------------------------------------------------#
# THURINGIA ----
#------------------------------------------------------------------------#
# Retrieve lines in the output file concerning only males involved in paternities
# idind = 2960 # Single  individual for testing
i = 7 # Select the run to estimate ancestry probabilities
individuals = unique(resThu$father) # List of individuals to retrieve (inferred fathers)

sampleThu = dataThu[which(dataThu$ageWhenFirstCaptur == "Adult"),]
dfThu = data.frame(indivID = sampleThu$idind,
                   PopID = as.numeric(sampleThu$idcol),
                   PopFlag = as.numeric(sampleThu$sexe == "F")
)
# Colonies are encoded in integers
# Translate from integers to names
colonynamesThu = levels(sampleThu$idcol)[as.numeric(sampleThu$idcol)]
# Add genotypes
dfThu = cbind(dfThu, sampleThu[,2:17])

# Run a for loop on the list of individuals
# Re-init result file
# system("rm STRUCTURE/output/Thu/results.txt")
# for (idind in individuals) {
#   # Copy the lines for individuals of interest in a result file
#   system(paste("grep ' ", idind, " ' STRUCTURE/output/Thu/run", i, "_f >> STRUCTURE/output/Thu/results.txt", sep = ""))
# }
# system("cat STRUCTURE/output/Thu/results.txt")

PopAssignThu = read.table("STRUCTURE/output/Thu/results.txt", header = FALSE)
PopAssignThu = PopAssignThu[,-c(1, 3, 5)] # Remove unwanted columns
colnames(PopAssignThu) = c("idind", "idcol", levels(as.factor(colonynamesThu))) # Rename columns
PopAssignThu$idcol = (levels(as.factor(colonynamesThu))[PopAssignThu$idcol])
PopAssignThu = PopAssignThu[,1:22]

# Now which is the inferred ancestry and at which probability
PopAssignThu$ancestry = NA
PopAssignThu$probancestry = NA
PopAssignThu$probancestry_second = NA


for (i in 1:nrow(PopAssignThu)) {
  PopAssignThu$ancestry[i] = as.character(colnames(PopAssignThu[,3:22])[which.max(PopAssignThu[i,3:22])])
  PopAssignThu$probancestry[i] = sort(as.numeric(PopAssignThu[i,3:22]), decreasing = T)[1]
  PopAssignThu$probancestry_second[i] = sort(as.numeric(PopAssignThu[i,3:22]), decreasing = T)[2]
}

# Assess the validity of the inference
# Is the probability of the second probable ancestry much lower than the first one?
(PopAssignThu$delta_proba = PopAssignThu$probancestry - PopAssignThu$probancestry_second)

# Is there natal dispersal?
# i.e. ancestry and colony are not the same
PopAssignThu$nataldisp = (PopAssignThu$idcol != PopAssignThu$ancestry)
table(PopAssignThu$nataldisp)

# 7 migrants for 11 non-migrants
mean(PopAssignThu$probancestry)

# Are probabilities of ancestry for migrant or non-migrant the same?
mean(PopAssignThu$probancestry[PopAssignThu$nataldisp == TRUE])
mean(PopAssignThu$probancestry[PopAssignThu$nataldisp == FALSE])

hist(PopAssignThu$probancestry[PopAssignThu$nataldisp == TRUE], breaks = 20, xlim = c(0, 0.7))
hist(PopAssignThu$probancestry[PopAssignThu$nataldisp == FALSE], breaks = 20, xlim = c(0, 0.7))

# Test it with Mann-Whitney
wilcox.test(PopAssignThu$probancestry[PopAssignThu$nataldisp == TRUE], PopAssignThu$probancestry[PopAssignThu$nataldisp == FALSE])
# Not significant for individuals of interest
# Besides probabilities associated to migrnats and n on-migrants are high
# especially compared to those of Picardy
# Maybe sampling is better and we have more polymorphism/structure in genetic markers
# More power to detect migrants

# Export results
# write.table(PopAssignThu, "Tables/MigrantFathersThu.txt", col.names = TRUE, row.names = FALSE,
#             quote = FALSE)
# PopAssignThu = read.table("Tables_supp/MigrantFathersThu.txt", header = TRUE)



#=================================================#
# Integrate uncertainty with model averaging across 10 MCMC chains ----
#=================================================#
# TODO Integrate uncertainty with model averaging across 10 MCMC chains:
# - Average across MCMC chains, kernel with:
# Min distance
# Max distance
# Mean distance


#--------------------------------------------------#
# PICARDY ----
#--------------------------------------------------#
# Retrieve lines in the output file concerning only males involved in paternities
# idind = 2960 # Single  individual for testing
individuals = unique(resPic$father) # List of individuals to retrieve (inferred fathers)

# Create the individuals infos
# Assign to pre- or post-dispersal
# pre- is for learning (adult females)
# post- is for assignment (adult males)
samplePic = dataPic[which(dataPic$ageWhenFirstCaptur == "Adult"),]
dfPic = data.frame(indivID = samplePic$idind,
                   PopID = as.numeric(samplePic$idcol),
                   PopFlag = as.numeric(samplePic$sexe == "F")
)
# Colonies are encoded in integers
# Translate from integers to names
colonynamesPic = levels(samplePic$idcol)[as.numeric(samplePic$idcol)]
# Add genotypes
dfPic = cbind(dfPic, samplePic[,2:17])


# Format results of the 10 MCMC chains
# for i in {1..10};
# do
# echo "Iteration number: $i"
# cat STRUCTURE/output/Pic/run${i}_f | grep "(0.000," > STRUCTURE/output/Pic_R1/run${i}_res.txt
# done

# for i in {1..10};
# do
# echo "Iteration number: $i"
# cat STRUCTURE/output/Pic/run${i}_f | grep "| Pop" > STRUCTURE/output/Pic_R1/run${i}_training.txt
# done

rm(res_pic)
for (i in 1:10) {
  res = read.table(paste0("STRUCTURE/output/Pic_R1/run", i, "_res.txt"))
  res = res[,-c(1, 3, 5)] # Remove unwanted columns
  colnames(res) = c("idind", "idcol", levels(as.factor(colonynamesPic))) # Rename columns
  res$idcol = (levels(as.factor(colonynamesPic))[res$idcol])
  res = res[,1:19]
  
  res$run = i
  
  if (exists("res_pic")) {
    res_pic = bind_rows(res_pic, res)
  } else {
    res_pic = res
  }
}


res_pic$col_origin = apply(res_pic[,3:19], 1, which.max)
res_pic$proba_col_origin = apply(res_pic[,3:19], 1, function(x) {x[which.max(x)]})

res_pic$run = as.factor(res_pic$run)


individuals = unique(resPic$father)
res_pic_sampled = res_pic %>%
  filter(idind %in% individuals)

res_pic_sampled_summary = res_pic_sampled %>%
  group_by(idind) %>%
  summarise(mean_proba = mean(proba_col_origin),
            n = n())


# How many different colonies are found?
# How many times the best colony is found?
res_pic_summary = res_pic %>%
  group_by(idind) %>%
  summarise(mean_proba = mean(proba_col_origin),
            n_col_found = length(unique(col_origin)),
            n = n())

res_pic_summary %>%
  filter(idind %in% individuals)

# Only one colony found across ten MCMC chains
# High consistency across chains


# FIGURES
figs4a = ggplot(res_pic, aes(x = proba_col_origin, colour = run)) +
  geom_density() +
  geom_vline(data = res_pic_sampled_summary, aes(xintercept = mean_proba), linetype =  "dashed", alpha = 0.4) +
  geom_density(data = res_pic_summary, aes(x = mean_proba, colour = "Black"), colour = "Black") +
  xlab("Assignment probability") +
  ylab("Density") +
  labs(colour = "Run") +
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
        legend.key.height = unit(1.5,"line"),
        legend.key.width = unit(1.5,"line"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14, face = "bold"),
        legend.position = "none")
figs4a


dir.create("Figures_R1")
ggsave("Figures_R1/Probability_assignment_per_run_Pic.jpeg", width = 8, height = 5)

#--------------------------------------------------#
# THURINGIA ----
#--------------------------------------------------#
# Retrieve lines in the output file concerning only males involved in paternities
# idind = 2960 # Single  individual for testing
individuals = unique(resThu$father) # List of individuals to retrieve (inferred fathers)

# Create the individuals infos
# Assign to pre- or post-dispersal
# pre- is for learning (adult females)
# post- is for assignment (adult males)
sampleThu = dataThu[which(dataThu$ageWhenFirstCaptur == "Adult"),]
dfThu = data.frame(indivID = sampleThu$idind,
                   PopID = as.numeric(sampleThu$idcol),
                   PopFlag = as.numeric(sampleThu$sexe == "F")
)
# Colonies are encoded in integers
# Translate from integers to names
colonynamesThu = levels(sampleThu$idcol)[as.numeric(sampleThu$idcol)]
# Add genotypes
dfThu = cbind(dfThu, sampleThu[,2:17])


# Format results of the 10 MCMC chains
# for i in {1..10};
# do
# echo "Iteration number: $i"
# cat STRUCTURE/output/Thu/run${i}_f | grep "(0.000," > STRUCTURE/output/Thu_R1/run${i}_res.txt
# done

# for i in {1..10};
# do
# echo "Iteration number: $i"
# cat STRUCTURE/output/Thu/run${i}_f | grep "| Pop" > STRUCTURE/output/Thu_R1/run${i}_training.txt
# done

rm(res_thu)
for (i in 1:10) {
  res = read.table(paste0("STRUCTURE/output/Thu_R1/run", i, "_res.txt"))
  res = res[,-c(1, 3, 5)] # Remove unwanted columns
  colnames(res) = c("idind", "idcol", levels(as.factor(colonynamesThu))) # Rename columns
  res$idcol = (levels(as.factor(colonynamesThu))[res$idcol])
  res = res[,1:22]
  
  res$run = i
  
  if (exists("res_thu")) {
    res_thu = bind_rows(res_thu, res)
  } else {
    res_thu = res
  }
}


res_thu$col_origin = apply(res_thu[,3:22], 1, which.max)
res_thu$proba_col_origin = apply(res_thu[,3:22], 1, function(x) {x[which.max(x)]})

res_thu$run = as.factor(res_thu$run)


individuals = unique(resThu$father)
res_thu_sampled = res_thu %>%
  filter(idind %in% individuals)

res_thu_sampled_summary = res_thu_sampled %>%
  group_by(idind) %>%
  summarise(mean_proba = mean(proba_col_origin))


# How many different colonies are found?
# How many times the best colony is found?
res_thu_summary = res_thu %>%
  group_by(idind) %>%
  summarise(mean_proba = mean(proba_col_origin),
            n_col_found = length(unique(col_origin)))

res_thu_summary %>%
  filter(idind %in% individuals)

# Only one colony found across ten MCMC chains
# High consistency across chains


# FIGURES
figs4b = ggplot(res_thu, aes(x = proba_col_origin, colour = run)) +
  geom_density() +
  geom_vline(data = res_thu_sampled_summary, aes(xintercept = mean_proba),
             linetype =  "dashed",
             alpha = 0.4) +
  geom_density(data = res_thu_summary, aes(x = mean_proba, colour = "Black"),
               colour = "Black") +
  xlab("Assignment probability") +
  ylab("Density") +
  labs(colour = "Run") +
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
        legend.text=element_text(size=14),
        legend.title=element_text(size=14, face = "bold"),
        legend.position = "right")

figs4b

ggsave("Figures_R1/Probability_assignment_per_run_Thu.jpeg", width = 8, height = 5)


# Fig S4 ----
ggpubr::ggarrange(figs4a, figs4b,
                  ncol = 2,
                  widths = c(2, 2.5))

ggsave("Figures_R1/FigS4.jpeg", width = 15, height = 5)


#=================================================#
# Integrate uncertainty across proba of the best MCMC chain ----
#=================================================#
# TODO Integrate uncertainty across proba of the best MCMC chain
# TODO Barplot of natal dispersal rate, for each rank of proba: first, second, third... up to six
# See how informative AND conservative (or not) we are
# TODO Dist of mean proba per rank: 1, 2, 3...


#--------------------------------------------------#
# PICARDY ----
#--------------------------------------------------#
# Run 2 
i = 2
res = res_pic %>%
  filter(run == i)

individuals = unique(resPic$father) # List of individuals to retrieve (inferred fathers)
res_pic %>%
  filter(idind %in% individuals)

# Identify colonies and associated proba of the five best assignments
rm(res_five_best)
for (i in 1:5) {
  res_i = res
  res_i$proba_rank = i
  
  res_i$proba_col_origin = apply(res_i[,3:19], 1, function(x) {sort(as.numeric(x), decreasing = T)[i]})
  res_i$col_origin = apply(res_i[,3:19], 1, function(x) {names(x)[order(as.numeric(x), decreasing = T)[i]]})
  
  if (exists("res_five_best")) {
    res_five_best = bind_rows(res_five_best, res_i)
  } else {
    res_five_best = res_i
  }
}

res_five_best$proba_rank = as.factor(res_five_best$proba_rank)

res_five_best$natal_disperser = (res_five_best$idcol != res_five_best$col_origin)

# FIGURES
fig3a = ggplot(res_five_best, aes(x = proba_col_origin, colour = proba_rank)) +
  geom_vline(data = res_pic_sampled_summary, aes(xintercept = mean_proba),
             linetype =  "dashed",
             alpha = 0.4) +
  geom_density() +
  xlab("Assignment probability") +
  ylab("Density") +
  labs(colour = "Probability\nRank") +
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
        legend.title=element_text(size=10, face = "bold"),
        legend.position = "none")
  # geom_density(data = res_pic_summary, aes(x = mean_proba, colour = "Black"), colour = "Black") +
  # theme_bw()

fig3a

fig3a = ggplot(res_five_best, aes(x = proba_rank, y = proba_col_origin, fill = proba_rank)) +
  geom_violin(width = 1) +
  geom_boxplot(width=0.2, fill = "white", alpha=0.2) +
  geom_hline(data = res_pic_sampled_summary, aes(yintercept = mean_proba),
             linetype =  "dashed",
             alpha = 0.2) +
  xlab("Probability rank") +
  ylab("Assignment probability") +
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
        legend.title=element_text(size=10, face = "bold"),
        legend.position = "none")

fig3a

ggsave("Figures_R1/Probability_assignment_rank_best_run_Pic.jpeg", width = 8, height = 5)


# TODO Dist of mean proba per rank: 1, 2, 3...
# The decay in proba for sampled fathers, showing power to discriminate the best colony
figs3a = ggplot(res_five_best %>% filter(idind %in% individuals), aes(y = proba_col_origin, fill = proba_rank, x = as.factor(idind))) +
  geom_col(position = "dodge")+
  theme_bw()  +
  ylab("Assignment\nprobability") +
  xlab("Father ID") +
  labs(fill = "Probability\nrank") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=10, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=10),
        axis.title.y = element_text(color="black", size=10),
        axis.text=element_text(size=10, colour="black"),
        legend.key = element_rect(fill = "white", colour = "white", size = 1),
        legend.key.height = unit(1.5,"line"),
        legend.key.width = unit(1.5,"line"),
        legend.text=element_text(size=6),
        legend.title=element_text(size=6, face = "bold"),
        legend.position = "right")

figs3a
ggsave("Figures_R1/Probability_assignment_rank_per_individual_Pic.jpeg", width = 8, height = 5)


# Barplot of natal dispersal rate, for each rank of proba: first, second, third... up to six ----
# See how informative AND conservative (or not) we are
res_five_best_summary = res_five_best %>%
  group_by(proba_rank) %>%
  summarise(natal_dispersal_events = sum(natal_disperser),
            n = n(),
            natal_dispersal_rate = natal_dispersal_events / n)

knitr::kable(res_five_best_summary)

res_sampled_five_best_summary = res_five_best %>%
  filter(idind %in% individuals) %>%
  group_by(proba_rank) %>%
  summarise(natal_dispersal_events = sum(natal_disperser),
            n = n(),
            natal_dispersal_rate = natal_dispersal_events / n)



figs2a = ggplot(res_five_best_summary, aes(y = natal_dispersal_rate, x = proba_rank)) +
  geom_col(position = "dodge")+
  theme_bw() +
  geom_hline(aes(yintercept = 1 - (1 / length(unique(colonynamesPic))))) +
  geom_hline(aes(yintercept = res_sampled_five_best_summary$natal_dispersal_rate[1]), linetype = "dashed") +
  ylim(0, 1)  +
  xlab("Probability rank") +
  ylab("Natal dispersal rate") +
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
        legend.title=element_text(size=10, face = "bold"),
        legend.position = "none")

figs2a
ggsave("Figures_R1/Natal_dispersal_rate_per_rank_Pic.jpeg", width = 8, height = 5)



# Barplot of natal dispersal rate for bins of proba of assignment ----
res$col_inferred = NA
col_names = colnames(res)[3:19]

for (i in 1:nrow(res)) {
  res$col_inferred[i] = col_names[res$col_origin[i]]
}

res$natal_disperser = (res$idcol != res$col_inferred)

res$proba_bin = ifelse(res$proba_col_origin >= 0.5, "p > 0.5",
                       ifelse(res$proba_col_origin >= 0.25 & res$proba_col_origin < 0.5, "0.25 < p < 0.5", 
                              ifelse(res$proba_col_origin >= 0.1 & res$proba_col_origin < 0.25, "0.1 < p < 0.25", "p < 0.1")))



res_bin_proba_summary = res %>%
  group_by(proba_bin) %>%
  summarise(natal_dispersal_events = sum(natal_disperser),
            n = n(),
            natal_dispersal_rate = natal_dispersal_events / n)

res_bin_proba_summary$proba_bin = factor(res_bin_proba_summary$proba_bin,
                       levels = c("p < 0.1",
                                  "0.1 < p < 0.25",
                                  "0.25 < p < 0.5"))

figs2c = ggplot(res_bin_proba_summary, aes(y = natal_dispersal_rate, x = proba_bin)) +
  geom_col(position = "dodge")+
  theme_bw() +
  geom_hline(aes(yintercept = 1 - (1 / length(unique(colonynamesPic))))) +
  geom_hline(aes(yintercept = res_sampled_five_best_summary$natal_dispersal_rate[1]), linetype = "dashed") +
  ylim(0, 1)  +
  xlab("Probability bin") +
  ylab("Natal dispersal rate") +
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
        legend.title=element_text(size=10, face = "bold"),
        legend.position = "none")


figs2c
ggsave("Figures_R1/Natal_dispersal_rate_per_proba_bin_Pic.jpeg", width = 8, height = 5)




# Save sampled fathers for the dispersal kernel
df = res_five_best %>%
  filter(idind %in% individuals)

dir.create("R1_supp_analyses")
write_tsv(df, "R1_supp_analyses/Pic_fathers_natal_colony.tsv")


# Identify colonies and associated proba of assignment for all colonies ----
rm(res_all_probas)
for (i in 1:length(unique(colonynamesPic))) {
  res_i = res
  res_i$proba_rank = i
  
  res_i$proba_col_origin = apply(res_i[,3:19], 1, function(x) {sort(as.numeric(x), decreasing = T)[i]})
  res_i$col_origin = apply(res_i[,3:19], 1, function(x) {names(x)[order(as.numeric(x), decreasing = T)[i]]})
  
  if (exists("res_all_probas")) {
    res_all_probas = bind_rows(res_all_probas, res_i)
  } else {
    res_all_probas = res_i
  }
}

res_all_probas$proba_rank = as.factor(res_all_probas$proba_rank)

table(res_all_probas$proba_rank)

res_all_probas$natal_disperser = (res_all_probas$idcol != res_all_probas$col_origin)

write_tsv(res_all_probas, "R1_supp_analyses/Pic_fathers_natal_colony_all_probas.tsv")




#--------------------------------------------------#
# THURINGIA ----
#--------------------------------------------------#
# Run 7 
i = 7
res = res_thu %>%
  filter(run == i)

individuals = unique(resThu$father) # List of individuals to retrieve (inferred fathers)

res_thu %>%
  filter(idind %in% individuals)

# Identify colonies and associated proba of the five best assignments
rm(res_five_best)
for (i in 1:5) {
  res_i = res
  res_i$proba_rank = i
  
  res_i$proba_col_origin = apply(res_i[,3:22], 1, function(x) {sort(as.numeric(x), decreasing = T)[i]})
  res_i$col_origin = apply(res_i[,3:22], 1, function(x) {names(x)[order(as.numeric(x), decreasing = T)[i]]})
  
  if (exists("res_five_best")) {
    res_five_best = bind_rows(res_five_best, res_i)
  } else {
    res_five_best = res_i
  }
}

res_five_best$proba_rank = as.factor(res_five_best$proba_rank)

res_five_best$natal_disperser = (res_five_best$idcol != res_five_best$col_origin)

# FIGURES
fig3b = ggplot(res_five_best, aes(x = proba_col_origin, colour = proba_rank)) +
  geom_vline(data = res_thu_sampled_summary, aes(xintercept = mean_proba),
             linetype =  "dashed",
             alpha = 0.4) +
  geom_density() +
  xlab("Assignment probability") +
  ylab("Density") +
  labs(colour = "Probability\nRank") +
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
  # geom_density(data = res_pic_summary, aes(x = mean_proba, colour = "Black"), colour = "Black") +
  # theme_bw()


fig3b = ggplot(res_five_best, aes(x = proba_rank, y = proba_col_origin, fill = proba_rank)) +
  geom_violin(width = 1) +
  geom_boxplot(width=0.2, fill = "white", alpha=0.2) +
  geom_hline(data = res_thu_sampled_summary, aes(yintercept = mean_proba),
             linetype =  "dashed",
             alpha = 0.2) +
  xlab("Probability rank") +
  ylab("Assignment probability") +
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
        legend.title=element_text(size=10, face = "bold"),
        legend.position = "none")


fig3b

ggsave("Figures_R1/Probability_assignment_rank_best_run_Thu.jpeg", width = 16, height = 10)



# Fig 3 ----
ggpubr::ggarrange(fig3a, fig3b,
                  ncol = 2)

ggsave("Figures_R1/Fig3.jpeg", width = 15, height = 5)


# TODO Dist of mean proba per rank: 1, 2, 3...
# The decay in proba for sampled fathers, showing power to discriminate the best colony
figs3b = ggplot(res_five_best %>% filter(idind %in% individuals),
       aes(y = proba_col_origin, fill = proba_rank, x = as.factor(idind))) +
  geom_col(position = "dodge")  +
  xlab("Father ID") +
  ylab("Assignment\nprobability") +
  labs(fill = "Probability\nrank") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=10, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=10),
        axis.title.y = element_text(color="black", size=10),
        axis.text=element_text(size=10, colour="black"),
        legend.key = element_rect(fill = "white", colour = "white", size = 1),
        legend.key.height = unit(1.5,"line"),
        legend.key.width = unit(1.5,"line"),
        legend.text=element_text(size=6),
        legend.title=element_text(size=6, face = "bold"),
        legend.position = "right")

figs3b
ggsave("Figures_R1/Probability_assignment_rank_per_individual_Thu.jpeg", width = 8, height = 5)


# Fig S3 ----
ggpubr::ggarrange(figs3a, figs3b, nrow = 2)

ggsave("Figures_R1/FigS3.jpeg", width = 10, height = 5)



# TODO Barplot of natal dispersal rate, for each rank of proba: first, second, third... up to six
# See how informative AND conservative (or not) we are
res_five_best_summary = res_five_best %>%
  group_by(proba_rank) %>%
  summarise(natal_dispersal_events = sum(natal_disperser),
            n = dplyr::n(),
            natal_dispersal_rate = natal_dispersal_events / n)

knitr::kable(res_five_best_summary)

res_sampled_five_best_summary = res_five_best %>%
  filter(idind %in% individuals) %>%
  group_by(proba_rank) %>%
  summarise(natal_dispersal_events = sum(natal_disperser),
            n = dplyr::n(),
            natal_dispersal_rate = natal_dispersal_events / n)



figs2b = ggplot(res_five_best_summary, aes(y = natal_dispersal_rate, x = proba_rank)) +
  geom_col(position = "dodge")+
  theme_bw() +
  geom_hline(aes(yintercept = 1 - (1 / length(unique(colonynamesThu))))) +
  geom_hline(aes(yintercept = res_sampled_five_best_summary$natal_dispersal_rate[1]), linetype = "dashed") +
  ylim(0, 1)  +
  xlab("Probability rank") +
  ylab("Natal dispersal rate") +
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
        legend.title=element_text(size=10, face = "bold"),
        legend.position = "none")

figs2b
ggsave("Figures_R1/Natal_dispersal_rate_per_rank_Thu.jpeg", width = 8, height = 5)


# Barplot of natal dispersal rate for bins of proba of assignment ----
res$col_inferred = NA
col_names = colnames(res)[3:22]

for (i in 1:nrow(res)) {
  res$col_inferred[i] = col_names[res$col_origin[i]]
}

res$natal_disperser = (res$idcol != res$col_inferred)

res$proba_bin = ifelse(res$proba_col_origin >= 0.5, "p > 0.5",
                       ifelse(res$proba_col_origin >= 0.25 & res$proba_col_origin < 0.5, "0.25 < p < 0.5", 
                              ifelse(res$proba_col_origin >= 0.1 & res$proba_col_origin < 0.25, "0.1 < p < 0.25", "p < 0.1")))



res_bin_proba_summary = res %>%
  group_by(proba_bin) %>%
  summarise(natal_dispersal_events = sum(natal_disperser),
            n = n(),
            natal_dispersal_rate = natal_dispersal_events / n)

res_bin_proba_summary$proba_bin = factor(res_bin_proba_summary$proba_bin,
                                         levels = c("p < 0.1",
                                                    "0.1 < p < 0.25",
                                                    "0.25 < p < 0.5",
                                                    "p > 0.5"))

figs2d = ggplot(res_bin_proba_summary, aes(y = natal_dispersal_rate, x = proba_bin)) +
  geom_col(position = "dodge")+
  theme_bw() +
  geom_hline(aes(yintercept = 1 - (1 / length(unique(colonynamesPic))))) +
  geom_hline(aes(yintercept = res_sampled_five_best_summary$natal_dispersal_rate[1]), linetype = "dashed") +
  ylim(0, 1)  +
  xlab("Probability rank") +
  ylab("Natal dispersal rate") +
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
        legend.title=element_text(size=10, face = "bold"),
        legend.position = "none")

figs2d
ggsave("Figures_R1/Natal_dispersal_rate_per_proba_bin_Thu.jpeg", width = 8, height = 5)


# Fig S2 ----
ggpubr::ggarrange(figs2a, figs2b, figs2c, figs2d, ncol = 2, nrow = 2)

ggsave("Figures_R1/FigS2.jpeg", width = 16, height = 8)



# Save sampled fathers for the dispersal kernel
df = res_five_best %>%
  filter(idind %in% individuals)

dir.create("R1_supp_analyses")
write_tsv(df, "R1_supp_analyses/Thu_fathers_natal_colony.tsv")


# Identify colonies and associated proba of assignment for all colonies ----
rm(res_all_probas)
for (i in 1:length(unique(colonynamesThu))) {
  res_i = res
  res_i$proba_rank = i
  
  res_i$proba_col_origin = apply(res_i[,3:22], 1, function(x) {sort(as.numeric(x), decreasing = T)[i]})
  res_i$col_origin = apply(res_i[,3:22], 1, function(x) {names(x)[order(as.numeric(x), decreasing = T)[i]]})
  
  if (exists("res_all_probas")) {
    res_all_probas = bind_rows(res_all_probas, res_i)
  } else {
    res_all_probas = res_i
  }
}

res_all_probas$proba_rank = as.factor(res_all_probas$proba_rank)

table(res_all_probas$proba_rank)

res_all_probas$natal_disperser = (res_all_probas$idcol != res_all_probas$col_origin)

write_tsv(res_all_probas, "R1_supp_analyses/Thu_fathers_natal_colony_all_probas.tsv")




#=================================================#
# Prior error rate based on training individuals ----
#=================================================#


#--------------------------------------------------#
# PICARDY ----
#--------------------------------------------------#






#--------------------------------------------------#
# THURINGIA ----
#--------------------------------------------------#



#=================================================#
# End -----
#=================================================#
