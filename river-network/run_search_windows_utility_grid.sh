#!/bin/bash

# data seed (controls the sample from the true model)
for pt1 in 29073 29074 29075 29049 28603 28604 28605 29039
do

	# job submit script name
	subfile='windows_'$pt1'.sh'
	#errorfile='windows_'${utilityname:0:1}'_'$j'_error.txt'
	#outputfile='windows_'${utilityname:0:1}'_'$j'_output.txt'
	
	#write submit script to file 
	cat > scripts/$subfile << EOF
#!/bin/bash -l
#
#PBS -N windows-clearwater-n2-$pt1
#PBS -l ncpus=10
#PBS -l cputype=6140
#PBS -l mem=10gb
#PBS -l walltime=45:00:00
#PBS -e scripts/windows-e.txt
#PBS -o scripts/windows-o.txt

cd \$PBS_O_WORKDIR
module load gdal/3.2.1-foss-2020b
module load r/4.0.3-foss-2020b
R -e "j <- 1; pt.1= '$pt1';source('./search_windows_utility_grid.R');"
EOF
	# submit the script to the HPC queue    
	#qsub -e scripts/$errorfile -o scripts/$outputfile scripts/$subfile
	qsub scripts/$subfile

done
