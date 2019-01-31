#!/bin/bash -l
##SBATCH --account=lcbc
##SBATCH --mail-type=END,FAIL,TIME_LIMIT                   # send mail when job is done for whatever reason
#SBATCH --mail-user=baudin.pablo@gmail.com
#SBATCH --error=%x.stde                                # Redirect std error file to <job-name>_<jobid>.stde
#SBATCH --output=%x.stdo                               # Redirect std output file to <job-name>_<jobid>.stdo
#SBATCH --constraint=E5v4

# JOBNAME

# PROG ENVIRONMENT

# WORK 

# SCRATCH


# PRINT PROG INFO
echo 'Executable: '$EXE
np="$(grep -c ^processor /proc/cpuinfo)"
echo 'Number of available cores: '$np

echo "Running job: " $job
echo "From working directory: " $wrk
echo "In scratch directory: " $scr
echo ""

# START SCRATCH FROM SCRATCH
rm -r $scr
mkdir -p $scr


# CP FILES TO SCR

# RUN PROG

# SAVE FILES FROM SCR


# send very elaborate email to myself
hst="$(hostname)"
tim="$(date)"
outlines="$(tail $out/${job}.out)"
tail $out/${job}.out

# put multiple line string in variable
read -r -d '' text <<- EOM
	Dear Dr. Baudin,
	
	Your calculation with job title: $job 
   running on: $hst 
   terminated on: $tim.

	Final outpul lines:
	
	$outlines
	
	I wish you a very pleasant day.
	
	Yours faithfully,
	
	Dr. Baudin ;)
EOM


# SEND EMAIL



