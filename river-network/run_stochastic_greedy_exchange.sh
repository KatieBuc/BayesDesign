#!/bin/bash
K=5
nmontedraws=500
nchoose=5
utilityname=KL.divergence
postrname=laplace
priors=frequentist
filename=stochastic_greedy_exchange

# data seed (controls the sample from the true model)
for (( j = 1; j <= $K; j++ ))
do
	# job submit script name
	subfile=''$filename'_'${utilityname:0:1}'_'$j'_'$nchoose'.sh'
	errorfile=''$filename'_'${utilityname:0:1}'_'$j'_'$nchoose'_error.txt'
	outputfile=''$filename'_'${utilityname:0:1}'_'$j'_'$nchoose'_output.txt'
	
	#write submit script to file 
	cat > scripts/$subfile << EOF
#!/bin/bash -l
#!/bin/bash -l
#
#PBS -N $filename-${utilityname:0:1}-$j-$nchoose
#PBS -l ncpus=10
#PBS -l cputype=6140
#PBS -l mem=150gb
#PBS -l walltime=80:00:00
#PBS -e scripts/$errorfile
#PBS -o scripts/$outputfile

cd \$PBS_O_WORKDIR
module load gdal/3.2.1-foss-2020b
module load r/4.0.3-foss-2020b
R -e "j <- $j; priors = '$priors'; n.monte.draws <- $nmontedraws; n.choose <- $nchoose; utility.function.name <- '$utilityname'; postr.approx.function.name <- '$postrname'; source('./$filename.R');"
EOF
	# submit the script to the HPC queue    
	qsub scripts/$subfile 

done
