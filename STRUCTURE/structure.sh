#!/bin/sh
# file: structure.sh

cd ./STRUCTURE/Pic

# Set programs as executables
chmod +x ./structure

# Launch batch runs for all values of K and every iterations
./structure -K 1 -o Results/K1run1 & ./structure -K 1 -o Results/K1run2 & ./structure -K 1 -o Results/K1run3
./structure -K 2 -o Results/K2run1 & ./structure -K 2 -o Results/K2run2 & ./structure -K 2 -o Results/K2run3
./structure -K 3 -o Results/K3run1 & ./structure -K 3 -o Results/K3run2 & ./structure -K 3 -o Results/K3run3
./structure -K 4 -o Results/K4run1 & ./structure -K 4 -o Results/K4run2 & ./structure -K 4 -o Results/K4run3
./structure -K 5 -o Results/K5run1 & ./structure -K 5 -o Results/K5run2 & ./structure -K 5 -o Results/K5run3
./structure -K 6 -o Results/K6run1 & ./structure -K 6 -o Results/K6run2 & ./structure -K 6 -o Results/K6run3
./structure -K 7 -o Results/K7run1 & ./structure -K 7 -o Results/K7run2 & ./structure -K 7 -o Results/K7run3
./structure -K 8 -o Results/K8run1 & ./structure -K 8 -o Results/K8run2 & ./structure -K 8 -o Results/K8run3
./structure -K 9 -o Results/K9run1 & ./structure -K 9 -o Results/K9run2 & ./structure -K 9 -o Results/K9run3
./structure -K 10 -o Results/K10run1 & ./structure -K 10 -o Results/K10run2 & ./structure -K 10 -o Results/K10run3

# STRUCTURE Harvester on results to:
# apply Evanno method (2005) to choose the best K
# prepare files for CLUMPP
python structureHarvester.py --dir=./Results/ --out=./Harvester/ --evanno --clumpp

# CLUMPP on files formated by STRUCTURE Harvester
# Launch a CLUMPP task on each K (except K=1, no segmentation possible on a single cluster)
chmod +x ./CLUMPP/CLUMPP
for (( i = 1; i < 11; i++ )); do
  ./CLUMPP/CLUMPP ./CLUMPP/paramfile -i './Harvester/K'$i'.indfile' -p './Harvester/K'$i'.popfile' -o './CLUMPP/K'$i'.outfile' -j './CLUMPP/K'$i'.miscfile' -k $i -c 484 -r 3
done


# Distruct to display nice graphs of STRUCTURE, with labels
# Automatically produce bar plots for each K value
chmod +x ./Distruct/distruct
for (( i = 0; i < 11; i++ )); do
  ./Distruct/distruct -K $i -M 18 -N 484 -p './CLUMPP/K'$i'.outfile' -i './Harvester/K'$i'.indfile' -o './Distruct/DistructPlotK'$i'.ps'
done
