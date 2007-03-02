#! /bin/bash 
#####################################################
# Trivial Runner Script 
# Suitable for running on scalar machines or places
# where no special invocation is needed for running
# parallel code
####################################################

retstring=`cqsub -e BGL_APP_L1_SWOA=1 -m vn -n 2 -c 4 -t 10:00 $*`
jobid=`echo $retstring | awk '{ print $4}'`
echo Job Id is $jobid

qstat_string="NOT_DONE"
while [ "X${qstat_string}X" != "XX" ];
do
	qstat_string=`cqstat | grep $jobid`
	job_state=`echo $qstat_string | awk '{ print $5}'`
	if [ "X${job_state}X" != "XX" ];
        then 
#	  echo Job state is $job_state. Sleeping for 10 seconds
	  sleep 10;
        fi
done
echo Job Finished
