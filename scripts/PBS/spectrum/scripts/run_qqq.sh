#!/bin/sh -xev
#
#  This script computes the qqq spectrum as a single job on the JLAB 
#  (or other) cluster

#
#  We begin with various PBS directives
#PBS -N __JOBNAME__
#PBS -l nodes=__NODES____ATTRIBUTE__

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
S_MASS=__S_MASS__
UD_MASS=__UD_MASS__

_RSDCG=1.0e-7
_MAXCG=1000

_APE_FACT=2.5
_APE_NUM=0
_GAUS_RAD_SRC=__GAUSS_RAD__
_GAUS_ITR_SRC=__GAUSS_ITR__
#
#  Where we store the data

#CMU
_ARCHROOT=/raid3/dgr/
_SERVER=qcdsvr:

#_ARCHROOT=/home/dgr/qcd/data/
#JLAB
#_ARCHROOT=/cache/LHPC/
#_SERVER=hpcdata3:
#
# Filenames
PROP_ROOT=prop
QQQ_PROPS=__QQQ_PROPS__

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

#
#  Here we have the launch command for the MYRINET cluster
#LAUNCH="/usr/local/mpich-gm-1.2.5/bin/mpirun -machinefile $PBS_NODEFILE"
#PROGDIR=$HOME/qcd/src/chroma/parscalar-gm/mainprogs/main
#SCP=rcp

#  gigE test cluster
#LAUNCH="/usr/local/qmp/bin/QMP_run.gige --qmp-l /etc/qmp/2x2x2-bigMem/3d_list_l --qmp-f /etc/qmp/2x2x2-bigMem/3d_conf_l"
#PROGDIR=$HOME/qcd/src/chroma/parscalar-gigE/mainprogs/main
#SCP=rcp

#  CMU cluster
LAUNCH=
PROGDIR=$HOME/qcd/src/chroma/scalar/mainprogs/main
SCP=/usr/bin/scp

#
#  Now the machine-dependent directory cluster
WORKDIR=/scratch/PBS_$PBS_JOBID

#
#  Begin by listing the nodes on which the job is running
echo
echo 'Job is running on nodes:'
echo ------------------------------------------------------
cat $PBS_NODEFILE
echo ------------------------------------------------------

#
#  Create the work directories on the scratch nodes

#pbs_create $PBS_NODEFILE $WORKDIR

cd $WORKDIR

echo "Working in $WORKDIR"
echo

#
#
#
#  We now loop over the configurations, copying data to the 
ctr=0
cfg=${CONFIG_START}
while [ ${ctr} -lt ${NCONFIG} ]
do

#
#  Copy over the gauge configuration

szingauge=${GAUGEROOT}${cfg}
${SCP} ${_SERVER}${_ARCHROOT}/${CFGLABEL}/cfgs/${szingauge} .

#
#  Now create the make source template file
#

cat > src_template <<EOF
<elem>
<Name>MAKE_SOURCE</Name>
<Frequency>1</Frequency>
<Param>
 <version>6</version>
 <Source>
 <version>1</version>
 <SourceType>SHELL_SOURCE</SourceType>
 <j_decay>3</j_decay>
 <t_srce>0 0 0 0</t_srce>

   <SmearingParam>
     <wvf_kind>GAUGE_INV_GAUSSIAN</wvf_kind>
     <wvf_param>${_GAUS_RAD_SRC}</wvf_param>
     <wvfIntPar>${_GAUS_ITR_SRC}</wvfIntPar>
     <no_smear_dir>3</no_smear_dir>
   </SmearingParam>

   <disp_length>_DISP_LENGTH</disp_length>
   <disp_dir>_DISP_DIR</disp_dir>
  
   <LinkSmearing>
      <LinkSmearingType>APE_SMEAR</LinkSmearingType>
      <link_smear_fact>${_APE_FACT}</link_smear_fact>
      <link_smear_num>${_APE_NUM}</link_smear_num>
      <no_smear_dir>3</no_smear_dir>
   </LinkSmearing>

 </Source>
 <nrow>${LSIZE} ${LSIZE} ${LSIZE} ${TSIZE}</nrow>
</Param>
<NamedObject>
 <source_id>_SOURCE_NAME</source_id>
</NamedObject>
</elem>
EOF

echo "Source template file is"
echo "***************************************************"
cat src_template
echo "***************************************************"

cat > prop_template <<EOF
<elem>
<Name>PROPAGATOR</Name>
<Frequency>1</Frequency>
<Param>
 <version>8</version>
 <FermTypeP>WILSON</FermTypeP>
 <nonRelProp>false</nonRelProp>
 <obsvP>false</obsvP>
 <FermionAction>
  <FermAct>WILSON</FermAct>
  <Mass>_MASS</Mass>
  <AnisoParam>
   <anisoP>true</anisoP>
   <t_dir>3</t_dir>
   <xi_0>${XI_0}</xi_0>
   <nu>${NU}</nu>
  </AnisoParam>
  <FermionBC>
  <FermBC>SIMPLE_FERMBC</FermBC>
  <boundary>1 1 1 -1</boundary>
  </FermionBC>
 </FermionAction>
 <InvertParam>
   <invType>CG_INVERTER</invType>
   <RsdCG>${_RSDCG}</RsdCG>
   <MaxCG>${_MAXCG}</MaxCG>
 </InvertParam>
 <nrow>${LSIZE} ${LSIZE} ${LSIZE} ${TSIZE}</nrow>
 <t_source>0 0 0 0</t_source>
</Param>
<NamedObject>
<source_id>_SOURCE_NAME</source_id>
<prop_id>_PROP_NAME</prop_id>
</NamedObject>
</elem>
EOF
echo "Propagator template file is"
echo "***************************************************"
cat prop_template
echo "***************************************************"


#
#  Now create the smearing template file
#

cat > snk_template <<EOF
<elem>
  <Name>SINK_SMEAR</Name>
  <Frequency>1</Frequency>
<Param>
 <version>5</version>
 <Sink>
 <version>1</version>
 <SinkType>SHELL_SINK</SinkType>
 <j_decay>3</j_decay>
   <SmearingParam>
     <wvf_kind>GAUGE_INV_GAUSSIAN</wvf_kind>
     <wvf_param>${_GAUS_RAD_SRC}</wvf_param>
     <wvfIntPar>${_GAUS_ITR_SRC}</wvfIntPar>
     <no_smear_dir>3</no_smear_dir>
   </SmearingParam>
     <disp_length>_DISP_LENGTH</disp_length>
     <disp_dir>_DISP_DIR</disp_dir>
 <LinkSmearing>
 <LinkSmearingType>APE_SMEAR</LinkSmearingType>
<link_smear_fact>${_APE_FACT}</link_smear_fact>
<link_smear_num>${_APE_NUM}</link_smear_num>
<no_smear_dir>3</no_smear_dir>
</LinkSmearing>
 </Sink>
 <nrow>${LSIZE} ${LSIZE} ${LSIZE} ${TSIZE}</nrow>
</Param>
<NamedObject>
 <prop_id>_PROP_NAME_IN</prop_id>
 <smeared_prop_id>_PROP_NAME_OUT</smeared_prop_id>
</NamedObject>
</elem>
<elem>
<Name>QIO_WRITE_NAMED_OBJECT</Name>
<Frequency>1</Frequency>
<NamedObject>
<object_id>_PROP_NAME_OUT</object_id>
<object_type>LatticePropagator</object_type>
</NamedObject>
<File>
<file_name>./_PROP_NAME_OUT</file_name>
<file_volfmt>MULTIFILE</file_volfmt>
</File>
</elem>
EOF
echo "Sink template file is"
echo "***************************************************"
cat snk_template
echo "***************************************************"

#
#  Now create the qqq template file
#

cat > qqq_template <<EOF
<elem>
<Name>QQQ</Name>
<Frequency>1</Frequency>
<Param>
 <version>4</version>
 <Dirac_basis>true</Dirac_basis>
 <nrow>${LSIZE} ${LSIZE} ${LSIZE} ${TSIZE}</nrow>
</Param>
<NamedObject>
<prop_ids>
<elem>_PROPAGATOR1</elem>
<elem>_PROPAGATOR2</elem>
<elem>_PROPAGATOR3</elem>
</prop_ids>
<qqq_file>_QQQNAME</qqq_file>
</NamedObject>
</elem>
EOF
echo "QQQ template file is"
echo "***************************************************"
cat qqq_template
echo "***************************************************"

#
# Create the cfg_template
cat > cfg_template <<EOF
</InlineMeasurements>
<nrow>${LSIZE} ${LSIZE} ${LSIZE} ${TSIZE}</nrow>
</Param>
<Cfg>
  <cfg_type>SZIN</cfg_type>
<cfg_file>${szingauge}</cfg_file>
</Cfg>
EOF
echo "CFG template file is"
echo "***************************************************"
cat cfg_template
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
#  Now create the input file for chroma

/usr/bin/perl $HOME/bin/parse_baryons.pl ${QQQ_PROPS} disp_table src_template prop_template ${UD_MASS} ${S_MASS} snk_template ${PROP_ROOT} $cfg qqq_template cfg_template

#
#
echo "Data input file is "
echo
echo --------------------------------------------------------------------
cat DATA
echo --------------------------------------------------------------------


#
#  Run the code
echo "Beginning execution of chroma at "`date`
echo
${LAUNCH} $PROGDIR/chroma
echo "Finished at "`date`
echo

echo "Files in working directory are "
echo -------------------------------------------------------------------
ls -l
echo -------------------------------------------------------------------
echo

#
#  Now copy the files back to the archive directory
#

${SCP} XMLDAT ${_SERVER}${_ARCHROOT}/${CFGLABEL}/xml/qqq_${PROP_ROOT}_xmldat.cfg${cfg}
${SCP} qqq_*.cfg${cfg} ${_SERVER}${_ARCHROOT}/${CFGLABEL}/qqq/.


rm qqq_*.cfg${cfg}

ctr=`expr ${ctr} + 1 `
cfg=`expr ${cfg} + ${CONFIG_SKIP} `

done




