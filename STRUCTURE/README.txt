INRA - ANR Pseudorasbora
Scott McCainrs - Thomas Brazier

README of STRUCTURE analysis

Thinning (default): 100

#===========================
# Directories

#---------------------------
# Worldwide populations
1a_AllPops: 2,112 loci, burnin 100,000 + 100,000 iterations
1b_AllPops: 2,112 loci, burnin 100,000 + 100,000 iterations, POPFLAG + PFROMPOPFLAGONLY
1c_AllPops: 2,112 loci, burnin 100,000 + 100,000 iterations, POPFLAG + PFROMPOPFLAGONLY, demes in Asia instead of sampled sites

#---------------------------
# Native populations
2a_Native: 2,112 loci, burnin 20,000 + 100,000 iterations
2b_Native: 2,112 loci, burnin 100,000 + 100,000 iterations
2b1_Native: 2,112 loci, burnin 100,000 + 100,000 iterations, subset of 13 populations in native area that showed admixture
2b_Native_OC: 2,112 loci, burnin 100,000 + 100,000 iterations --> only convergent (OC) retained
2c_Native: 3,000 loci, burnin 100,000 + 100,000 iterations

#---------------------------
# Invasive populations
3a_Invasive: 2,112 loci, burnin 100,000 + 100,000 iterations, thinning of 1
3a1_Invasive: 2,112 loci, burnin 100,000 + 100,000 iteration, subset of Western Europe
3b_Invasive: 3,000 loci, burnin 100,000 + 100,000 iterations
