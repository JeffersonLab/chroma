#!/bin/sh -xev
#
# $Id: run_complete_job.sh,v 3.0 2006-04-03 04:58:40 edwards Exp $
#


###############################################################
#                                                             #
#    Bourne shell script for submitting a parallel MPICH job  #
#    to the PBS queue using the qsub command.                 #
#                                                             #
###############################################################

#     Remarks: A line beginning with # is a comment.
#	       A line beginning with #PBS is a PBS directive.
#              PBS directives must come first; any directives
#                 after the first executable statement are ignored.
#
   
##########################
#                        #
#   The PBS directives   #
#                        #
##########################

#          Set the name of the job (up to 15 characters, 
#          no blank spaces, start with alphanumeric character)

#PBS -N __JOBNAME__

#          Specify the number of nodes requested and the
#          number of processors per node. 

#PBS -l nodes=1:ppn=1

#          By default, the standard output and error streams are sent
#          to files in the current working directory with names:
#              job_name.osequence_number  <-  output stream
#              job_name.esequence_number  <-  error stream
#          where job_name is the name of the job and sequence_number 
#          is the job number assigned when the job is submitted.
#          Use the directives below to change the files to which the
#          standard output and error streams are sent.

#    #PBS -o stdout_file
#    #PBS -e stderr_file

#          The directive below directs that the standard output and
#          error streams are to be merged, intermixed, as standard
#          output. 

#PBS -j oe

#          Specify the maximum cpu and wall clock time. The wall
#          clock time should take possible queue waiting time into
#          account.  Format:   hhhh:mm:ss   hours:minutes:seconds
#          Be sure to specify a reasonable value here.
#          If the job does not finish by the time reached,
#          the job is terminated.

#PBS -l     cput=12:00:00
#PBS -l walltime=12:00:00

#          Specify the queue.  The CMU cluster currently has two queues:
#          "fermion" and "boson".  Jobs submitted to the "fermion" queue
#          will run in cpu-dedicated mode; if all cpu's assigned to the
#          queue are occupied with a job, then new jobs are queued and will
#          not run until a cpu is freed up.  You should take this waiting
#          time into account when setting "walltime".  Jobs submitted to
#          the "boson" queue will run in cpu-shared mode.  Such jobs will
#          generally begin execution immediately.

#PBS -q fermion_gold

#          Specify the maximum amount of physical memory required.
#          kb for kilobytes, mb for megabytes, gb for gigabytes.
#          Take some care in setting this value.  Setting it too large
#          can result in your job waiting in the queue for sufficient
#          resources to become available.

#PBS -l mem=512mb

#          PBS can send informative email messages to you about the
#          status of your job.  Specify a string which consists of
#          either the single character "n" (no mail), or one or more
#          of the characters "a" (send mail when job is aborted),
#          "b" (send mail when job begins), and "e" (send mail when
#          job terminates).  The default is "a" if not specified.
#          You should also specify the email address to which the
#          message should be send via the -M option.

#  #PBS -m abe
#  #PBS -m ae

#PBS -M user_email_address

#          Declare the time after which the job is eligible for execution.
#          If you wish the job to be immediately eligible for execution,
#          comment out this directive.  If you wish to run at some time in 
#          future, the date-time argument format is
#                      [DD]hhmm
#          If the day DD is not specified, it will default to today if the
#          time hhmm is in the future, otherwise, it defaults to tomorrow.
#          If the day DD is specified as in the future, it defaults to the
#          current month, otherwise, it defaults to next month.

# #PBS -a 2215  commented out

#          Specify the priority for the job.  The priority argument must be
#          an integer between -1024 and +1023 inclusive.  The default is 0.

#  #PBS -p 0

#          Define the interval at which the job will be checkpointed,
#          if checkpointing is desired, in terms of an integer number
#          of minutes of CPU time.

#  #PBS -c c=2

##########################################
#                                        #
#   Output some useful job information.  #
#                                        #
##########################################

#
#  Some job-dependent run parameters
CONFIG_START=__CONFIG_START__
CONFIG_SKIP=__CONFIG_SKIP__
NCONFIG=__NCONFIG__

#
#  Physics sizes and ensemble parameters
LSIZE=__LSIZE__
TSIZE=__TSIZE__
CFGLABEL=__CFGLABEL__
GAUGEROOT=__GAUGEROOT__
#
# Fermion parameters
NU=__NU__
XI_0=__XI_0__
MASS=__MASS__

_RSDCG=1.0e-7
_MAXCG=1000

_APE_FACT=2.5
_APE_NUM=5
_GAUS_RAD_SRC=3.6
_GAUS_ITR_SRC=32
#
#  Where we store the data
_ARCHROOT=/raid3/dgr/

#
# Filenames
_UD_PROP_ROOT=ud_prop280_test
_S_PROP_ROOT=s_prop280_test


NCPU=`wc -l < $PBS_NODEFILE`
echo ------------------------------------------------------
echo ' This job is allocated on '${NCPU}' cpu(s)'
echo 'Job is running on node(s): '
cat $PBS_NODEFILE
echo ------------------------------------------------------
echo PBS: qsub is running on $PBS_O_HOST
echo PBS: originating queue is $PBS_O_QUEUE
echo PBS: executing queue is $PBS_QUEUE
echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: execution mode is $PBS_ENVIRONMENT
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: node file is $PBS_NODEFILE
echo PBS: current home directory is $PBS_O_HOME
echo PBS: PATH = $PBS_O_PATH
echo ------------------------------------------------------

##############################################################
#                                                            #
#   The prologue script automatically makes a directory      #
#   on the local disks for you.  The name of this directory  #
#   depends on the job id, but you need only refer to it     #
#   using ${WORKDIR}.  In the mpirun command below, the      #
#   machines are specified using the copy of the             #
#   $PBS_NODEFILE which is placed on the /scratch partition  #
#   of the local disk by the prologue script.  You could     #
#   instead use                                              #
#            -machinefile $PBS_NODEFILE                      #
#                                                            #
##############################################################

SERVER=$PBS_O_HOST
WORKDIR=/scratch/PBS_$PBS_JOBID
SCP=/usr/bin/scp
SSH=/usr/bin/ssh
MACHINES=${WORKDIR}/NODEFILE
#LAUNCH="$HOME/mpich-1.2.5/ch_p4/bin/mpirun -np $NCPU -machinefile $PBS_NODEFILE "
LAUNCH=

######################################################################
#                                                                    #
#   To minimize communications traffic, it is best for your job      #
#   to work with files on the local disk of the compute node.        #
#   Hence, one needs to transfer files from your permanent home      #
#   directory tree to the directory ${WORKDIR} automatically         #
#   created by PBS on the local disk before program execution,       #
#   and to transfer any important output files from the local        #
#   disk back to the permanent home directory tree after program     #
#   execution is completed.                                          #
#                                                                    #
#   There are essentially two ways to achieve this: (1) to use the   #
#   PBS stagein and stageout utilities, or (2) to manually copy the  #
#   files by commands in this script.  The stagein and stageout      #
#   features of OpenPBS are somewhat awkward, especially since       #
#   wildcards and macros in the file lists cannot be used.  This     #
#   method also has some timing issues.  Hence, we ask you to use    #
#   the second method, and to use secure copy (scp) to do the file   #
#   transfers to avoid NSF bottlenecks.                              #
#                                                                    #
######################################################################

#####################################################
#                                                   #
#    Specify the permanent directory(ies) on the    #
#    server host.  Note that when the job begins    #
#    execution, the current working directory at    #
#    the time the qsub command was issued becomes   #
#    the current working directory of the job.      #
#                                                   #
#####################################################


PERMDIR=${HOME}/work
#PROGDIR=$HOME/bin/qdp++-1.9.3
#PROGDIR=$HOME/qcd/src/chroma/scalar/mainprogs/main
PROGDIR=$HOME/bin/scalar

PROGLIST="make_source propagator sink_smearing qqq_w"


SERVPERMDIR=${PBS_O_HOST}:${PERMDIR}

echo server is $SERVER
echo workdir is $WORKDIR
echo permdir is $PERMDIR
echo servpermdir is $SERVPERMDIR
echo ------------------------------------------------------
echo 'Job is running on node(s): '
cat $PBS_NODEFILE
echo ------------------------------------------------------
echo ' '
echo ' '

###############################################################
#                                                             #
#    Transfer files from server to local disks.               #
#                                                             #
###############################################################

stagein()
{
 if [ -r $MACHINES ] ; then
    machines=$(sort $MACHINES | uniq )
 else
    machines=$(hostname)
 fi

 for machine in $machines ; do

    echo ' '
    echo Transferring files from server to compute node ${machine}
    echo Writing files in node directory ${WORKDIR}

    for names in $PROGLIST
    do
      ${SCP} ${SERVER}:${PROGDIR}/${names} ${machine}:${WORKDIR}
    done

    ${SCP} ${SERVER}:/raid3/baryons/simulation_input/First_Run/Nucleons/qqq_props ${WORKDIR}

    echo Files in node work directory are as follows:
    ${SSH} ${machine} ls -l ${WORKDIR}
    
 done
}

############################################################
#                                                          #
#    Execute the run.  Do not run in the background.       #
#                                                          #
############################################################

runprogram()
{
 cd ${WORKDIR}
}

###########################################################
#                                                         #
#   Copy necessary files back to permanent directory.     #
#                                                         #
###########################################################

stageout()
{
 echo ' '
 echo Transferring files from compute nodes to server
 echo Writing files in permanent directory  ${PERMDIR}
 cd ${WORKDIR}

 #${SCP} output_files  ${SERVPERMDIR}

 echo Final files in permanent data directory:
 ${SSH} ${SERVER} "cd ${PERMDIR}; ls -l"
 }

#####################################################################
#                                                                   #
#  The "qdel" command is used to kill a running job.  It first      #
#  sends a SIGTERM signal, then after a delay (specified by the     #
#  "kill_delay" queue attribute (set to 60 seconds), unless         #
#  overridden by the -W option of "qdel"), it sends a SIGKILL       #
#  signal which eradicates the job.  During the time between the    #
#  SIGTERM and SIGKILL signals, the "cleanup" function below is     #
#  run. You should include in this function commands to copy files  #
#  from the local disk back to your home directory.  Note: if you   #
#  need to transfer very large files which make take longer than    #
#  60 seconds, be sure to use the -W option of qdel.                #
#                                                                   #
#####################################################################

early()
{
 echo ' '
 echo ' ############ WARNING:  EARLY TERMINATION #############'
 echo ' '
 }

trap 'early; stageout' 2 9 15


##################################################
#                                                #
#   Staging in, running the job, and staging out #
#   were specified above as functions.  Now      #
#   call these functions to perform the actual   #
#   file transfers and program execution.        #
#                                                #
##################################################

stagein
cd ${WORKDIR}

#
#  We now loop over the configurations, copying data to the 
ctr=0
cfg=${CONFIG_START}
while [ ${ctr} -lt ${NCONFIG} ]
do

#
#  Copy over the gauge configuration

szingauge=${GAUGEROOT}${cfg}

${SCP} ${SERVER}:${_ARCHROOT}/${CFGLABEL}/${LSIZE}^3x${TSIZE}/cfgs/${szingauge} .

#
#  Now create the make source template file
#

cat > src_template <<EOF
<?xml version="1.0"?>

<make_source>
<annotation>
; $Id: run_complete_job.sh,v 3.0 2006-04-03 04:58:40 edwards Exp $
;
; MAKE_SOURCE input file.
;
; This program is the input file for a  make_source  test run on Wilson-type
; propagators
;
</annotation>

<Param>
 <version>5</version>
 <wave_state>S_WAVE</wave_state>
 <source_type>SHELL_SOURCE</source_type>
 <j_decay>3</j_decay>
 <direction>2</direction>
 <t_source>0 0 0 0</t_source>

 <ShellSource>
   <SourceSmearingParam>
     <wvf_kind>GAUGE_INV_GAUSSIAN</wvf_kind>
     <wvf_param>${_GAUS_RAD_SRC}</wvf_param>
     <wvfIntPar>${_GAUS_ITR_SRC}</wvfIntPar>
   </SourceSmearingParam>
   <laplace_power>0</laplace_power>
   <link_smear_fact>${_APE_FACT}</link_smear_fact>
   <link_smear_num>${_APE_NUM}</link_smear_num>
   <disp_length>_DISP_LENGTH</disp_length>
   <disp_dir>_DISP_DIR</disp_dir>
 </ShellSource>
 <nrow>${LSIZE} ${LSIZE} ${LSIZE} ${TSIZE}</nrow>
</Param>
<Cfg>
 <cfg_type>SZIN</cfg_type>
 <cfg_file>${szingauge}</cfg_file>
</Cfg>
<Prop>
 <source_file>_SOURCE_NAME</source_file>
 <source_volfmt>MULTIFILE</source_volfmt>
</Prop>
</make_source>
EOF

echo "Source template file is"
echo "***************************************************"
cat src_template
echo "***************************************************"

#
#  Now create the propagator template file
#

cat > prop_template <<EOF
<?xml version="1.0"?>

<propagator>
<annotation>
; $Id: run_complete_job.sh,v 3.0 2006-04-03 04:58:40 edwards Exp $
;
; PROPAGATOR input file.
;
; This program is the input file for a propagator (spectroscopy)
; test of the preconditioned Clover operator.
;
; This program is the input file for a propagator test run on Wilson-type
; propagators
</annotation>

<Param>
 <version>6</version>
 <FermTypeP>WILSON</FermTypeP>
 <nonRelProp>false</nonRelProp>
 <FermionAction>
  <FermAct>WILSON</FermAct>
  <Mass>${MASS}</Mass>
  <AnisoParam>
   <anisoP>true</anisoP>
   <t_dir>3</t_dir>
   <xi_0>${XI_0}</xi_0>
   <nu>${NU}</nu>
  </AnisoParam>
 </FermionAction>
 <InvertParam>
   <invType>CG_INVERTER</invType>
   <RsdCG>${_RSDCG}</RsdCG>
   <MaxCG>${_MAXCG}</MaxCG>
 </InvertParam>
 <nrow>${LSIZE} ${LSIZE} ${LSIZE} ${TSIZE}</nrow>
 <boundary>1 1 1 -1</boundary>
 <t_source>0 0 0 0</t_source>
</Param>
<Cfg>
 <cfg_type>SZIN</cfg_type>
 <cfg_file>${szingauge}</cfg_file>
</Cfg>
<Prop>
 <source_file>_SOURCE_NAME</source_file>
 <prop_file>_PROP_NAME</prop_file>
 <prop_volfmt>SINGLEFILE</prop_volfmt>
</Prop>
</propagator>
EOF

echo "Propagator template file is"
echo "***************************************************"
cat prop_template
echo "***************************************************"

#
#  Now create the smearing template file
#

cat > snk_template <<EOF
<?xml version="1.0"?>

<sink_smearing>
<annotation>
; 
;
; SINK_SMEARING input file
;
; This program is the input file for a sink_smearing
;
</annotation>

<Param>
 <version>4</version>
 <wave_state>S_WAVE</wave_state>
 <sink_type>SHELL_SINK</sink_type>
 <direction>2</direction>

 <ShellSink>
   <SinkSmearingParam>
     <wvf_kind>GAUGE_INV_GAUSSIAN</wvf_kind>
     <wvf_param>${_GAUS_RAD_SRC}</wvf_param>
     <wvfIntPar>${_GAUS_ITR_SRC}</wvfIntPar>
   </SinkSmearingParam>
   <laplace_power>0</laplace_power>
   <link_smear_fact>${_APE_FACT}</link_smear_fact>
   <link_smear_num>${_APE_NUM}</link_smear_num>
   <disp_length>_DISP_LENGTH</disp_length>
   <disp_dir>_DISP_DIR</disp_dir>
 </ShellSink>

 <nrow>${LSIZE} ${LSIZE} ${LSIZE} ${TSIZE}</nrow>
</Param>
<Cfg>
 <cfg_type>SZIN</cfg_type>
 <cfg_file>${szingauge}</cfg_file>
</Cfg>
<Prop>
 <prop_file>_PROP_NAME_IN</prop_file>
 <smeared_prop_file>_PROP_NAME_OUT</smeared_prop_file>
 <smeared_prop_volfmt>MULTIFILE</smeared_prop_volfmt>
</Prop>
</sink_smearing>
EOF

echo "Sink template file is"
echo "***************************************************"
cat snk_template
echo "***************************************************"

#
#  Now create the qqq template file
#

cat > qqq_template <<EOF
<?xml version="1.0"?>

<qqq_w>
<annotation>
; : qqq_w_v4.xml.ini,v 1.1 2004/06/16 14:51:35 dgr Exp $
;
; QQQ_W input file.
;
; This program is the input file for binding up bunch of propagators
; into a generalized propagator
;
</annotation>

<Param>
 <version>3</version>
 <Dirac_basis>true</Dirac_basis>
 <nrow>${LSIZE} ${LSIZE} ${LSIZE} ${TSIZE}</nrow>
</Param>
<Cfg>
 <cfg_type>SZIN</cfg_type>
 <cfg_file>${szingauge}</cfg_file>
</Cfg>
<Prop>
 <prop_file>
<elem>_PROPAGATOR1</elem>
<elem>_PROPAGATOR2</elem>
<elem>_PROPAGATOR3</elem>
</prop_file>
</Prop>
<Barcomp>
 <qqq_file>_QQQNAME</qqq_file>
</Barcomp>
</qqq_w>
EOF

echo "QQQ template file is"
echo "***************************************************"
cat qqq_template
echo "***************************************************"

cat > disp_table <<EOF
1
2
3
EOF
echo "Displacement table is "
echo "***************************************************"
cat disp_table
echo "***************************************************"



#
#  We now create all the DATA files we need for the propagators

$HOME/bin/parse_qs.pl qqq_props disp_table src_template prop_template prop_template snk_template ${_UD_PROP_ROOT} $cfg

#
#  We now create all the DATA files we need for the qqq files

$HOME/bin/parse_qqq.pl qqq_props disp_table qqq_template ${_UD_PROP_ROOT} $cfg

#
#  We now generate the sources and run the propagators

echo "Generated all the input files"
echo
echo "Files in working directory are"
ls -l
echo

prop_ctr=0

#  First check the file exists...
while [ -f src_ini.$prop_ctr ]
do

mv src_ini.$prop_ctr DATA
echo "Making source number $prop_ctr at "`date`
echo "DATA file is "
echo "*********************************************************"
cat DATA
echo "*********************************************************"
echo

${LAUNCH} ./make_source
echo "..Done at "`date`

${SCP} XMLDAT ${SERVER}:${_ARCHROOT}/${CFGLABEL}/${LSIZE}^3x${TSIZE}/xml/prop_${_UD_PROP_ROOT}_src_xmldat.${prop_ctr}.cfg${cfg}
rm XMLDAT

mv prop_ini.${prop_ctr} DATA

echo "Making ud propagator number $prop_ctr at "`date`
echo "DATA file is "
echo "*********************************************************"
cat DATA
echo "*********************************************************"
echo
${LAUNCH} ./propagator
echo "..Done at "`date`

${SCP} XMLDAT ${SERVER}:${_ARCHROOT}/${CFGLABEL}/${LSIZE}^3x${TSIZE}/xml/prop_${_UD_PROP_ROOT}_xmldat.${prop_ctr}.cfg${cfg}
rm XMLDAT

echo "Files in working directory are "
ls -l 
echo "*******************************************************************"

echo "Finished source and prop number $prop_ctr for config ${cfg}"
echo

prop_ctr=`expr $prop_ctr + 1 `


done

echo Copying local sink propagator to archive at `date`
${SCP} *PS*.cfg${cfg} ${SERVER}:${_ARCHROOT}/${CFGLABEL}/${LSIZE}^3x${TSIZE}/solver/.

echo

echo "Files on workdir are"
ls -l

# Now compute the sink smearing
#
prop_ctr=0

#  First check the file exists...
while [ -f sink_ini.$prop_ctr ]
do

mv sink_ini.${prop_ctr} DATA
echo "Smearing number $prop_ctr at "`date`
echo "DATA file is "
echo "*********************************************************"
cat DATA
echo "*********************************************************"
echo

${LAUNCH} ./sink_smearing

${SCP} XMLDAT ${SERVER}:${_ARCHROOT}/${CFGLABEL}/${LSIZE}^3x${TSIZE}/xml/sink_${_UD_PROP_ROOT}_xmldat.${prop_ctr}.cfg${cfg}
rm XMLDAT

echo "Finished sink smearing $prop_ctr for config ${cfg}"
echo

prop_ctr=`expr $prop_ctr + 1 `

done

echo Files on disk are
ls -l

#
#  Now we need to compute the generalised propagators

qqq_ctr=0

#  First check the file exists...
while [ -f qqq_ini.${qqq_ctr} ]
do

mv qqq_ini.${qqq_ctr} DATA

echo "Running qqq at "`date`
echo "DATA file is "
echo "*********************************************************"
cat DATA
echo "*********************************************************"
echo
${LAUNCH} ./qqq_w

${SCP} XMLDAT ${SERVER}:${_ARCHROOT}/${CFGLABEL}/${LSIZE}^3x${TSIZE}/xml/qqq_${_UD_PROP_ROOT}_${_S_PROP_ROOT}.${qqq_ctr}.cfg${cfg}
rm XMLDAT

echo "Finished qqq number ${qqq_ctr} for config ${cfg}"
echo

qqq_ctr=`expr $qqq_ctr + 1 `

done

echo "Moving qqq files to the appropriate places"

${SCP} qqq_*.cfg${cfg} ${SERVER}:${_ARCHROOT}/${CFGLABEL}/${LSIZE}^3x${TSIZE}/qqq/.

rm qqq_*.cfg${cfg}

ctr=`expr ${ctr} + 1 `
cfg=`expr ${cfg} + ${CONFIG_SKIP} `

done

#stageout 

###############################################################
#                                                             #
#   The epilogue script automatically deletes the directory   #
#   created on the local disk (including all files contained  #
#   therein.                                                  #
#                                                             #
###############################################################

