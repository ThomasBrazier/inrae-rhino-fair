# FAIR code and data to reproduce the results of The geometry of gametic dispersal in a flying mammal



## Contact


Thomas Brazier
brazier.thomas@gmail.com
thomas.brazier@univ-rennes.fr



## Open and FAIR code and data policy

As of March 2025, PCI Ecology and PCI Evol Biol are implementing data and code curation for new submissions. Submissions are expected to fill all the criteria listed below to be recommended:

1. Is the data accessible, interoperable, and associated with appropriate metadata and licenses?

2. Is the data as described by the manuscript?

3. Is the code accessible, interoperable, and associated with appropriate metadata and licenses?

4. Does the code match the methods?

5. Does the code run without error?

6. Does running the code with the provided data produce the same results as written in the manuscript?


## Data and files

```
Data/Pic/coordPic.txt
Data/Pic/distancesDist.txt
Data/Pic/distancesDistPic_gametic.txt
Data/Pic/dyadsObs.txt
Data/Pic/dyadsObsSelect.txt
Data/Pic/resultsWithInfo.txt
Data/Pic/uniqueGenotypesWithInfo.txt
```

```
Data/Thu/coordThu.txt
Data/Thu/distancesDist.txt
Data/Thu/distancesDistThu_gametic.txt
Data/Thu/dyadsObsSelect.txt
Data/Thu/resultsWithInfo.txt
Data/Thu/uniqueGenotypesWithInfo.txt
```


`Analyses.zip`

`Data.zip`



## Reproduce analyses

I assume that you run these scripts on a Linux system. Before running each script, set your working directory to the script location.

Because some parts of the workflow are interactive, necessitate to copy and paste files manually, or launch commands manually, but also they are computationally intensive, I tag them as *DO NOT RE-RUN UNLESS NECESSARY*.

In addition, some scripts are experimental, and their outputs are not presented in the final manuscript, hence they are tagged as *EXPERIMENTAL*.




The script `data_for_manuscript.R` is the central script to reproduce all the results presented in the paper. In addition, figures are produced by the script `Figures_preprint.R`.



The scripts to reproduce analyses must be run in this order:

(1) IBD

The script `IBD.R` analyses Isolation by Distance and produces `Fig. S2`. The intercept and slope of the IBD regression line are also calculated and a Mantel test is performed.

Basic summary statistics of genetic diversity are computed in `GeneticDiversity.R` *EXPERIMENTAL* and not presented in the final manuscript.


(2) COLONY

The COLONY analysis is performed in `AssignationPic` and `AssignationThu` for France and Germany, respectively.

The complete workflow is in these directories:


* `AssignationThu/01_prepadata/` prepare the dataset *DO NOT RE-RUN UNLESS NECESSARY*
* `AssignationThu/01_prepadata/01a_ExtractData.R` Get the data and format the dataset *DO NOT RE-RUN UNLESS NECESSARY*
* `AssignationThu/02_colony/` where to run COLONY *DO NOT RE-RUN UNLESS NECESSARY*
* `AssignationThu/02_colony/02_constructColony.R` Make the input files, config files and shell command files to launch COLONY *DO NOT RE-RUN UNLESS NECESSARY*
* `AssignationThu/03_simulations/` Where to run COLONY on simulations *DO NOT RE-RUN UNLESS NECESSARY*
* `AssignationThu/03_simulations/03_simulations.R` Prepare input for COLONY simulation module *DO NOT RE-RUN UNLESS NECESSARY*
* `AssignationThu/04_data/` The final data produced
* `AssignationThu/AnalyseSimulations.R` to produce the COLONY empirical output dataset for downstream analyses
* `AssignationThu/AnalyseColony.R` to produce the COLONY simulation output dataset for downstream analyses




For Thuringia (Thu) output files are stored in `AssignationThu/outputs` and the scripts to analyse and reformat the results are `AssignationThu/AnalyseColony.R` and `AssignationThu/AnalyseSimulations.R`, for the empirical data and the simulations, respectively. In addition,`AssignationThu/fonctionsColony.R` and `AssignationThu/fonctionsSim.R` are helper functions. The same organization applies to the Picardy (Pic) dataset.



!!! IMPORTANT !!!



NOTE The COLONY workflow is quite complex because it involved many manual steps and is also intensive in computation, hence I recommend not to run these scripts unless you wish to perform a new complete analysis. In addition, many comments are in french (sorry!). All the outputs and results files are already stored in the appropriate directories to run downstream analyses.



(3) MasterBayes

`MasterBayes\MasterBayes.R` is a standalone script to run the MasterBayes analysis.
All input files are already in the directory.


NOTE that `MasterBayes` does not seem to be installable anymore (removed from CRAN). All the commands, formulas and hyperparameters necessary top re-run MasterBayes are commmented and tabulated. Moreover, the models inferred with MasterBayes are saved in `MasterBayes\` and can be loaded and interpreted in `MasterBayes\MasterBayes.R`.


(4) STRUCTURE

Everything needed to run the STRUCTURE analysis is in `STRUCTURE Dispersal.R`  *DO NOT RE-RUN UNLESS NECESSARY*. Inputs, outputs and the STRUCTURE software and config files are stored in `./STRUCTURE`.


NOTE that STRUCTURE is computationally expensive. Run it on a cluster if necessary, or redo analyses from the outputs already in `./STRUCTURE`.


(5) Dispersal distances

In these scripts, I analyse the outputs of COLONY and STRUCTURE to compute mating, natal and gametic dispersal distances per offspring.

These scripts are exploratory and are not presented in the manuscript.


* `Dispersal_Kernel_Gametic.R` *EXPERIMENTAL*
* `Dispersal_Kernel.R` *EXPERIMENTAL*
* `Gametic dispersal.R` where I added the natal origin of the father to recompute the distribution of gametic dispersal distances
* `Natal dispersal.R` *EXPERIMENTAL*



(6) FitDistR


```
fitDistR_model_comparison_Gametic reparameterized.R
fitDistR_model_comparison_Gametic.R
fitDistR_model_comparison.R
```

Outputs:

`fitDistR_model_comparison.R` output the estimated parameters of the mating dispersal kernel. The two following tables are necessary for `Fig. 2`.

```
Tables/ParamPic.txt # The best fit Kernel on the Mating Dispersal distances in Picardy
Tables/ParamThu.txt # The best fit Kernel on the Mating Dispersal distances in Thuringia
```

(7) Figures and Tables

Run `Figures_preprint.R` to generate the figures presented in the paper.

```
Fig 1
Fig 2
Fig 3
Fig S1
```

The following script produces the results presented in the paper: `data_for_manuscript.R`

The output generated by each analysis and presented in the manuscript are stored in `Tables\`.

```
Table 1
Table S1
```



## Reproduce Tables and Figures

Figures are stored in `Figures\` and tables in `Tables\`.

The script to reproduce all the figures in the main text is `Figures_preprint.R`.