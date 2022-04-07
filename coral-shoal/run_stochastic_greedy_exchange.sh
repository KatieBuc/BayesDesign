#!/bin/bash
K=1
#utilityname=D.posterior.precision
utilityname=KL.divergence

# data seed (controls the sample from the true model)
for (( j = 1; j <= $K; j++ ))
do

	# job submit script name
	subfile='stochasticgreedy_'${utilityname:0:1}'_'$j'.sh'
	errorfile='stochasticgreedy_'${utilityname:0:1}'_'$j'_error.txt'
	outputfile='stochasticgreedy_'${utilityname:0:1}'_'$j'_output.txt'
	
	#write submit script to file 
	cat > scripts/$subfile << EOF
#!/bin/bash -l
#
#PBS -N CORAL-GREEDY-${utilityname:0:1}-$j-MEDIAN
#PBS -l ncpus=10
#PBS -l cputype=6140
#PBS -l mem=150gb
#PBS -l walltime=2:00:00
#PBS -e scripts/$errorfile
#PBS -o scripts/$outputfile

cd \$PBS_O_WORKDIR
module load gdal/3.2.1-foss-2020b
module load r/4.0.3-foss-2020b
R -e "j <- $j; source('./stochastic_greedy_exchange.R');"
EOF
	# submit the script to the HPC queue    
	qsub scripts/$subfile

done
