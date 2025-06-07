##########################################################################
###!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!###
#========================================================================#
#			EXTRACTING RESULTS FROM STRUCTURE
#========================================================================#
###!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!###
##########################################################################

#========================================================================#
# Loop to Extract Posterior Estimates & Write to Dataframe
#------------------------------------------------------------------------#
#!#!#! N.B. the names of standard output files will vary with each Structure run
	# e.g. outputK1run1.txt
		# 'output' will be standard
		# 'K' & 'run' will be semi-standard
			# integers will vary based on the number of 'K' defined in parameter set & the number of replicates (i.e. runs)

#------------------------------------------------------------------------#
# Alpha & posterior estimates of P(D) (sampling chain) -- full chain with burnin
#------------------------------------------------------------------------#

#------------------------------------------------------------------------#
# Alpha & posterior estimates of P(D) (sampling chain) -- full chain with burnin
#------------------------------------------------------------------------#
# This part get values of Alpha, Ln Likelihood of P(Data) for each MCMC iteration (burnin + whole length of MCMC), for all replicates (runs) of a given K and write it into the file trace_alphaK${k}
# !!! One sampling chain per cluster
# Supplementary data: Values of F, D and r

#------------------------------------------------------------------------#
# Init K value and number of replicates (number of runs)
# With K in {k1-k2}
# Values are taken in arguments
k1=17
k2=17
run=1

for k in $(seq $k1 $k2);
do
  (echo Run; (grep -A 1 'BURNIN completed' outputK${k}run${run}.txt | grep -v 'BURNIN completed')) | paste - - > sampling_chain_K${k} # Init the file where to save the trace of alpha and the sampling chain, with names of columns
  # Reformating space characters
  sed -i s/'Ln Like'/'Ln.Like'/'' sampling_chain_K${k}
  sed -i s/'Est Ln P(D)'/'Est.Ln.P'/'' sampling_chain_K${k}
  # print the number of fields
  fields="$(awk '{print NF}' sampling_chain_K${k})"

  for r in $(seq $run)
  do
    (grep -A 11 'Rep#:' outputK${k}run${r}.txt | grep -v 'Rep#:' \
  	| grep '[^[:blank:]]' | grep '[[:blank:]$]') | awk 'BEGIN {FS=OFS="\t"; r='$r'};  {print r, (seq $fields)}' >> sampling_chain_K${k}
  done
  # Reformatting iterations: removing ':' character
  sed -i s/':'//'g' sampling_chain_K${k}
  # and ',' characters
  sed -i s/','//'g' sampling_chain_K${k}
  # Reformating BURNIN iterations: those with '--' became '--    --'
  sed -i s/'--'/'--  --'/'' sampling_chain_K${k}

  sed -i 's/\s/#/g' sampling_chain_K${k} # visualise how many spaces
  # Not elegant, but I needed to move on...
  sed -i 's/--##--##/NA NA/' sampling_chain_K${k}
  sed -i 's/#A/ A/g' sampling_chain_K${k}
  sed -i 's/###0./ 0./' sampling_chain_K${k}
  sed -i 's/##0./ 0./' sampling_chain_K${k}
  sed -i 's/########/ /' sampling_chain_K${k}
  sed -i 's/######/ /' sampling_chain_K${k}
  sed -i 's/####/ /g' sampling_chain_K${k}
  sed -i 's/Like##/Like /' sampling_chain_K${k}
  sed -i 's/###/ /g' sampling_chain_K${k}
  sed -i 's/##/ /g' sampling_chain_K${k}
  sed -i 's/#/ /g' sampling_chain_K${k}

  sed -i 's/\s\s/ /g' sampling_chain_K${k}

  head sampling_chain_K${k}
done
