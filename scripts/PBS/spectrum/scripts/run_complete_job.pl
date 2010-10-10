#!/usr/bin/perl
#
#  Perl version of CMU job script (created by David) for JLAB clusters:
#  Construction of generalized propagator using displacement operators
#  and using Colin's source-operator tables .... Subhasish, 02/18/05
#
#  ###########################################################

#  read in command line arguments
#
die "Usage: $0 <Start Cfg> <Skip Cfg> <N Cfgs> <pbs_name>\n" unless $#ARGV eq 3;

#  stdin assignments
#
$_START_CFG=$ARGV[0];
$_SKIP_CFG=$ARGV[1];
$_NCFG=$ARGV[2];
$_PBS_NAME=$ARGV[3];

#  site and node parameters
#  JLAB-4G (128), JLAB-3G (256) or JLAB-2M (64)
#  for test-cluster use: JLAB-3G and 8 nodes w/ keyword MESH
#
$_SITE="JLAB-3G";
site:	{
    (${_SITE} eq "JLAB-2M") && do{
	$_NODES=64; 
	$_WALLTH=12; $_CPUTH=${_NODES} * ${_WALLTH};
	last site;
    };
    (${_SITE} eq "JLAB-3G") && do{
	$_PLANE=2; $_NODES=256; 
	$_WALLTH=12; $_CPUTH=${_NODES} * ${_WALLTH};
	last site;
    };
    (${_SITE} eq "JLAB-4G") && do{
	$_PANEL=01; $_NODES=128;	#  panel can be 01, 23 or 45
	$_WALLTH=12; $_CPUTH=${_NODES} * ${_WALLTH};
	last site;
    };
    die "Unsupported site!\n";
	}

#  list of chroma executables
#
$_MAKE_SRC=make_source;
$_PROPAGATOR=propagator;
$_SINK_SMR=sink_smearing;
$_SPECTRUM=spectrum_w;
$_QQQ=qqq_w;

#  link smearing types, executable and parameters
#  APE, HYP or STOUT
#
$_SMEAR=APE;
smear:	{
     (${_SMEAR} eq APE) && do{
	$_LINK_SMR=apesmear;
        $_SMEAR_FACT=2.5; $_SMEAR_NUM=5;
        last smear;
     };
     (${_SMEAR} eq HYP) && do{
	$_LINK_SMR=hypsmear3d;
        $_ALPHA1=0.72; $_ALPHA2=0.6; $_ALPHA3=0.3; $_SMEAR_NUM=2;
        last smear;
     };
     (${_SMEAR} eq STOUT) && do{
	$_LINK_SMR=stoutsmear;
        $_SMEAR_FACT=0.2; $_SMEAR_NUM=2;
        last smear;
     };
     die "Unsupported smearing type!\n";
	}

#  rest of the simulation parameters
#
$_LSIZE=8;
$_TSIZE=48;
$_ANISO=true;

$_GAUS_RAD_SRC=3.6;
$_GAUS_ITR_SRC=32;
$_GAUS_RAD_SNK=3.6;
$_GAUS_ITR_SNK=32;

$_LAP_POW_SRC=0;
$_LAP_POW_SNK=0;

$_SRC_TYPE=SHELL_SOURCE;
$_SNK_TYPE=SHELL_SINK;
$_SRC_WAVE_TYPE=S_WAVE;
$_SNK_WAVE_TYPE=S_WAVE;
$_SRC_NDIR=0;
$_SNK_NDIR=0;

$_RSDCG=1.0e-7;
$_MAXITER=10000;

mass:	{
    (${_LSIZE} eq 4) && do{
        $_MASS=0.11;
  	last mass;
    };
    (${_LSIZE} eq 8) && do{
        $_MASS=-0.373; $_XI_0=2.464; $_NU=0.9700;
	last mass;
    };
    ((${_LSIZE} eq 16)||(${_LSIZE} eq 24)) && do{
        $_MASS=-0.313; $_XI_0=2.464; $_NU=0.8980;
	last mass;
    };
    die "Unsupported lattice size!\n";
	}

# labeling stuff
#
$_UD_PROP_ROOT=ud_propagator;
$_S_PROP_ROOT=s_propagator;

#  ###########################################################

#  setting up header for pbs job
#
open(OUT,">${_PBS_NAME}.sh");
print OUT <<EOF;
#!/bin/sh
#PBS -V
#PBS -r n
#PBS -l file=1000mb
#PBS -l nodes=${_NODES}
#PBS -l walltime=${_WALLTH}:00:00
#PBS -l cput=${_CPUTH}:00:00
#PBS -N ${_PBS_NAME}
#PBS -S /bin/bash
set -xev
echo -------------------------------
cat \$PBS_NODEFILE
echo -------------------------------
set -xev
EOF
close(OUT);

#  creating rest of the script
#
open(OUT,">>${_PBS_NAME}.sh");
print OUT <<EOF1;

SITE=${_SITE}
RCP=rcp

#  machine specific instructions
#
case \${SITE} in
JLAB-2M)
MACHINE_TYPE=parscalar-gm
MPIRUN=/usr/local/mpich-gm-1.2.5/bin/mpirun
;;
JLAB-3G)
MACHINE_TYPE=parscalar-gige
MPIRUN=/usr/local/qmp/bin/QMP_run.sh
if [ ${_NODES} -eq 8 ]
then
  QMP_CONF_FILE=/etc/qmp/2x2x2/3d_conf_l
  QMP_LIST_FILE=/etc/qmp/2x2x2/3d_list_l
elif [ ${_NODES} -eq 256 ]
then
  QMP_CONF_FILE=/etc/qmp/4x4x16/3d_conf_l
  QMP_LIST_FILE=/etc/qmp/4x4x16/3d_list_l
elif [ ${_NODES} -eq 64 ]
then
  QMP_CONF_FILE=/etc/qmp/4x8x8/p${_PLANE}/2d_conf_l
  QMP_LIST_FILE=/etc/qmp/4x8x8/p${_PLANE}/2d_list_l
else
echo " ${_NODES} not supported "
exit 1
fi
;;
JLAB-4G)
MACHINE_TYPE=parscalar-gige
MPIRUN=/usr/local/qmp/bin/QMP_run.sh
if [ ${_NODES} -eq 128 ]
then
  QMP_CONF_FILE=/etc/qmp/4x4x8-${_PANEL}/3d_conf_l
  QMP_LIST_FILE=/etc/qmp/4x4x8-${_PANEL}/3d_list_l
elif [ ${_NODES} -eq 384 ]
then
  QMP_CONF_FILE=/etc/qmp/6x8x8/3d_conf_l
  QMP_LIST_FILE=/etc/qmp/6x8x8/3d_list_l
else
echo " ${_NODES} not supported "
exit 1
fi
;;
esac
echo "Using executable type \${MACHINE_TYPE} at ${_SITE}"

#  servers & root directories
#
WORKROOT=/scratch
RESULTROOT=//cache/users/sbasak
FILESERVER=qcdi01:
MAINSERVER=qcdi01
BINROOT=\${HOME}/bin
LIBROOT=\${HOME}/lib
PBSROOT=\${HOME}/pbs
ARCHROOT=/cache/LHPC/NF0

case ${_LSIZE} in
4)
_CFGLABEL=iso/${_LSIZE}^3x${_TSIZE}_wl
_GAUGEROOT=test_purgaug.cfg;;
8)
_CFGLABEL=aniso/5p8_${_LSIZE}_${_TSIZE}_xi3p0_wl
_GAUGEROOT=wqa_5p8_${_LSIZE}_${_TSIZE}_xi3p0.cfg;;
16)
_CFGLABEL=aniso/6p1_${_LSIZE}_${_TSIZE}_xi3p0_wl
_GAUGEROOT=wqa_6p1_${_LSIZE}_${_TSIZE}_x3p0.cfg;;
24)
_CFGLABEL=aniso/6p1_${_LSIZE}_${_TSIZE}_xi3p0_wl
_GAUGEROOT=wlq_6p1_${_LSIZE}_${_TSIZE}_x3p0.cfg;;
esac
GAUGEDIR=\${ARCHROOT}/\${_CFGLABEL}/cfgs
PROGROOT=\${HOME}/chroma/\${MACHINE_TYPE}/mainprogs/main
XMLDIR=\${RESULTROOT}/xml_${_LSIZE}_${_TSIZE}
QQQDIR=\${RESULTROOT}/qqq_${_LSIZE}_${_TSIZE}
SPCDIR=\${RESULTROOT}/spc_${_LSIZE}_${_TSIZE}

#  ###########################################################

#  machine specific files
#
case \${SITE} in
JLAB-2M)
PBSPROGS="pbs_cleanup pbs_create pbs_stagein pbs_stageout pbs_replicate pbs_replicate_pwd"
PROGLIST="${_LINK_SMR} ${_MAKE_SRC} ${_PROPAGATOR} ${_SINK_SMR} ${_SPECTRUM} ${_QQQ}"
;;
esac

#  Creating working directories
#
WORKDIR=\${WORKROOT}
cd \${WORKDIR}

#  machine specific initialization
#
case \${SITE} in
JLAB-2M)
#  Copying over the PBS utilities
for name in \${PBSPROGS}
do
echo "Copying over ... " \${name}
\${RCP} \${FILESERVER}\${PBSROOT}/\${name} .
done
conf_file=conf.\${PBS_JOBID}
for name in \${PROGLIST}
do
echo "Copying over and replicating program "  \${name}
\${RCP} \${FILESERVER}\${PROGROOT}/\${name} .
./pbs_replicate \${PBS_NODEFILE} \${WORKDIR} \${name}
done
echo "Finished setting up \${SITE}"
;;
esac

\${RCP} \${FILESERVER}\${LIBROOT}/uud_props .

#  ###########################################################

#  ##################################
#  begin the loop over configurations
#
ctr=0
cfg=${_START_CFG}
while [ \${ctr} -lt ${_NCFG} ]
do

echo
echo "*******************************"
echo "Begining configuration  \${cfg}"
echo "*******************************"

szingauge=\${_GAUGEROOT}\${cfg}

#  test gauge config file exists or not and act accordingly
#
GG_RESULT=\$(rsh \${MAINSERVER} ls \${GAUGEDIR}/\${szingauge})
if [ \${GG_RESULT} ]; then

echo -n "Copying over \${szingauge} ..."
\${RCP} \${FILESERVER}\${GAUGEDIR}/\${szingauge} .
echo "files are "
ls -lrt

#  ########################################
#  creating the link smearing template file
#
case ${_SMEAR} in
APE)
cat > lsmr_template <<EOF
<?xml version="1.0"?>

<apesmear>
<annotation>
;
; Ape smearing input
;
</annotation>

<Param>
  <version>2</version>
  <nrow>${_LSIZE} ${_LSIZE} ${_LSIZE} ${_TSIZE}</nrow>
  <j_decay>3</j_decay>
  <link_smear_fact>${_SMEAR_FACT}</link_smear_fact>
  <link_smear_num>${_SMEAR_NUM}</link_smear_num>
</Param>
<Cfg>
  <cfg_type>SZIN</cfg_type>
  <cfg_file>\${szingauge}</cfg_file>
</Cfg>
<Ape>
  <volfmt>MULTIFILE</volfmt>
  <ape_type>SZINQIO</ape_type>
  <ape_file>smr_\${szingauge}</ape_file>
</Ape>

</apesmear>
EOF
;;
HYP)
cat > lsmr_template <<EOF
<?xml version="1.0"?>

<hypsmear3d>
<annotation>
;
; Hyp smearing input
;
</annotation>

<Param>
  <version>2</version>
  <alpha1>${_ALPHA1}</alpha1>
  <alpha2>${_ALPHA2}</alpha2>
  <alpha3>${_ALPHA3}</alpha3>
  <link_smear_num>${_SMEAR_NUM}</link_smear_num>
  <nrow>${_LSIZE} ${_LSIZE} ${_LSIZE} ${_TSIZE}</nrow>
  <j_decay>3</j_decay>
</Param>
<Cfg>
  <cfg_type>SZIN</cfg_type>
  <cfg_file>\${szingauge}</cfg_file>
</Cfg>
<Hyp>
  <volfmt>MULTIFILE</volfmt>
  <hyp_type>SZINQIO</hyp_type>
  <hyp_file>smr_\${szingauge}</hyp_file>
</Hyp>

</hypsmear3d>
EOF
;;
STOUT)
cat > lsmr_template <<EOF
<?xml version="1.0"?>

<stoutsmear>
<annotation>
;
; Stout smearing input
;
</annotation>

<Param>
  <version>2</version>
  <nrow>${_LSIZE} ${_LSIZE} ${_LSIZE} ${_TSIZE}</nrow>
  <j_decay>3</j_decay>
  <link_smear_fact>${_SMEAR_FACT}</link_smear_fact>
  <link_smear_num>${_SMEAR_NUM}</link_smear_num>
</Param>
<Cfg>
  <cfg_type>SZIN</cfg_type>
  <cfg_file>\${szingauge}</cfg_file>
</Cfg>
<Stout>
  <volfmt>MULTIFILE</volfmt>
  <stout_type>SZINQIO</stout_type>
  <stout_file>smr_\${szingauge}</stout_file>
</Stout>

</stoutsmear>
EOF
;;
esac

echo "LINK template file is ..."
echo "************************************************"
cat lsmr_template
echo "************************************************"

#  ######################################
#  creating the make_source template file
#
cat > src_template <<EOF
<?xml version="1.0"?>

<make_source>
<annotation>
;
; Make_source input
;
</annotation>

<Param>
 <version>5</version>
 <wave_state>${_SRC_WAVE_TYPE}</wave_state>
 <source_type>${_SRC_TYPE}</source_type>
 <j_decay>3</j_decay>
 <direction>${_SRC_NDIR}</direction>
 <t_source>0 0 0 0</t_source>

 <ShellSource>
   <SourceSmearingParam>
     <wvf_kind>GAUGE_INV_GAUSSIAN</wvf_kind>
     <wvf_param>${_GAUS_RAD_SRC}</wvf_param>
     <wvfIntPar>${_GAUS_ITR_SRC}</wvfIntPar>
   </SourceSmearingParam>
   <laplace_power>${_LAP_POW_SRC}</laplace_power>
   <link_smear_fact>0.0</link_smear_fact>
   <link_smear_num>0</link_smear_num>
   <disp_length>_DISP_LENGTH</disp_length>
   <disp_dir>_DISP_DIR</disp_dir>
 </ShellSource>

 <nrow>${_LSIZE} ${_LSIZE} ${_LSIZE} ${_TSIZE}</nrow>
</Param>
<Cfg>
 <cfg_type>SZIN</cfg_type>
 <cfg_file>\${szingauge}</cfg_file>
</Cfg>
<Prop>
 <source_file>_SOURCE_NAME</source_file>
 <source_volfmt>MULTIFILE</source_volfmt>
</Prop>
</make_source>
EOF

echo "SOURCE template file is ..."
echo "************************************************"
cat src_template
echo "************************************************"

#  #################################
#  creating propagator template file
#
cat > prop_template <<EOF
<?xml version="1.0"?>

<propagator>
<annotation>
;
;  Propagator input
;
</annotation>

<Param>
 <version>6</version>
 <FermTypeP>WILSON</FermTypeP>
 <nonRelProp>false</nonRelProp>
 <FermionAction>
  <FermAct>WILSON</FermAct>
  <Mass>${_MASS}</Mass>
  <AnisoParam>
    <anisoP>${_ANISO}</anisoP>
    <t_dir>3</t_dir>
    <xi_0>${_XI_0}</xi_0>
    <nu>${_NU}</nu>
  </AnisoParam>
 </FermionAction>
 <InvertParam>
   <invType>CG_INVERTER</invType>
   <RsdCG>${_RSDCG}</RsdCG>
   <MaxCG>${_MAXITER}</MaxCG>
 </InvertParam>
 <nrow>${_LSIZE} ${_LSIZE} ${_LSIZE} ${_TSIZE}</nrow>
 <boundary>1 1 1 -1</boundary>
 <t_source>0 0 0 0</t_source>
</Param>
<Cfg>
 <cfg_type>SZIN</cfg_type>
 <cfg_file>\${szingauge}</cfg_file>
</Cfg>
<Prop>
 <source_file>_SOURCE_NAME</source_file>
 <prop_file>_PROP_NAME</prop_file>
 <prop_volfmt>SINGLEFILE</prop_volfmt>
</Prop>
</propagator>
EOF

echo "PROPAGATOR template file is ..."
echo "************************************************"
cat prop_template
echo "************************************************"

#  ####################################
#  creating sink_smearing template file
#
cat > snk_template <<EOF
<?xml version="1.0"?>

<sink_smearing>
<annotation>
;
; Sink_smearing input
;
</annotation>

<Param>
 <version>4</version>
 <wave_state>${_SNK_WAVE_TYPE}</wave_state>
 <sink_type>${_SNK_TYPE}</sink_type>
 <direction>${_SNK_NDIR}</direction>

 <ShellSink>
   <SinkSmearingParam>
     <wvf_kind>GAUGE_INV_GAUSSIAN</wvf_kind>
     <wvf_param>${_GAUS_RAD_SNK}</wvf_param>
     <wvfIntPar>${_GAUS_ITR_SNK}</wvfIntPar>
   </SinkSmearingParam>
   <laplace_power>${_LAP_POW_SNK}</laplace_power>
   <link_smear_fact>0.0</link_smear_fact>
   <link_smear_num>0</link_smear_num>
   <disp_length>_DISP_LENGTH</disp_length>
   <disp_dir>_DISP_DIR</disp_dir>
 </ShellSink>

 <nrow>${_LSIZE} ${_LSIZE} ${_LSIZE} ${_TSIZE}</nrow>
</Param>
<Cfg>
 <cfg_type>SZIN</cfg_type>
 <cfg_file>\${szingauge}</cfg_file>
</Cfg>
<Prop>
 <prop_file>_PROP_NAME_IN</prop_file>
 <smeared_prop_file>_PROP_NAME_OUT</smeared_prop_file>
 <smeared_prop_volfmt>MULTIFILE</smeared_prop_volfmt>
</Prop>
</sink_smearing>
EOF

echo "SINK template file is ..."
echo "************************************************"
cat snk_template
echo "************************************************"

#  ############################
#  creating qqq_w template file
#
cat > qqq_template <<EOF
<?xml version="1.0"?>

<qqq_w>
<annotation>
;
; QQQ_W input
;
</annotation>

<Param>
 <version>3</version>
 <Dirac_basis>true</Dirac_basis>
 <nrow>${_LSIZE} ${_LSIZE} ${_LSIZE} ${_TSIZE}</nrow>
</Param>
<Cfg>
 <cfg_type>SZIN</cfg_type>
 <cfg_file>\${szingauge}</cfg_file>
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

echo "QQQ template file is ..."
echo "************************************************"
cat qqq_template
echo "************************************************"

cat > disp_table <<EOF
1
2
3
EOF

echo "Displacement table is ..."
echo "************************************************"
cat disp_table
echo "************************************************"

#  ###################################################
#  creating all DATA files we need for the propagators
#
\${BINROOT}/parse_qs.pl uud_props disp_table src_template prop_template prop_template snk_template ${_UD_PROP_ROOT} \${cfg}

#  #############################################
#  creating all DATA files we need for the qqq_w
#
\${BINROOT}/parse_qqq.pl uud_props disp_table qqq_template ${_UD_PROP_ROOT} \${cfg}

#  generate the sources and run the propagators
#
echo "Generated all the input files"
echo
echo "Files in working directory are"
ls -lrt
echo

#  #################################################  #
#                                                     #
#  Link smearing execution needs to be included here  #
#                                                     #
#  #################################################  #

prop_ctr=0
#  first check the file exists ...
#
while [ -f src_ini.\${prop_ctr} ]
do

#  #####################
#  executing make_source
#
mv src_ini.\${prop_ctr} DATA
echo "Making source number \${prop_ctr} at " \`date\`
echo "DATA file is "
echo "*********************************************************"
cat DATA
echo "*********************************************************"
echo

sleep 40
echo "Running ${_MAKE_SRC} at " \`date\`
case \${SITE} in
JLAB-2M)
(\${MPIRUN} -np ${_NODES} -machinefile \${PBS_NODEFILE} ${_MAKE_SRC}) >& RESULT_mksrc
;;
JLAB-3G | JLAB-4G)
(\${MPIRUN} --qmp-f \${QMP_CONF_FILE} --qmp-l \${QMP_LIST_FILE} \${PROGROOT}/${_MAKE_SRC}) >& RESULT_mksrc
;;
esac
echo "... done at " \`date\`
mv XMLDAT prop_src_xmldat.\${prop_ctr}.cfg\${cfg}

#  ####################
#  executing propagator
#
mv prop_ini.\${prop_ctr} DATA
echo "Making ud propagator number \${prop_ctr} at " \`date\`
echo "DATA file is "
echo "*********************************************************"
cat DATA
echo "*********************************************************"
echo

sleep 40
echo "Running ${_PROPAGATOR} at " \`date\`
case \${SITE} in
JLAB-2M)
(\${MPIRUN} -np ${_NODES} -machinefile \${PBS_NODEFILE} ${_PROPAGATOR}) >& RESULT_prop
;;
JLAB-3G | JLAB-4G)
(\${MPIRUN} --qmp-f \${QMP_CONF_FILE} --qmp-l \${QMP_LIST_FILE} \${PROGROOT}/${_PROPAGATOR}) >& RESULT_prop
;;
esac
echo "... done at " \`date\`
echo "files are "
ls -lrt
mv XMLDAT prop_xmldat.\${prop_ctr}.cfg\${cfg}

echo "Finished source and prop number \${prop_ctr} for config \${cfg}"
echo

prop_ctr=\`expr \${prop_ctr} + 1 \`
done    #  end of make_source and propagator

prop_ctr=0
#  first check the file exists ...
#
while [ -f sink_ini.\${prop_ctr} ]
do

#  ###########################
#  computing the sink smearing
#
mv sink_ini.\${prop_ctr} DATA
echo "Smearing number \${prop_ctr} at " \`date\`
echo "DATA file is "
echo "*********************************************************"
cat DATA
echo "*********************************************************"
echo

sleep 40
echo "Running ${_SINK_SMR} at " \`date\`
case \${SITE} in
JLAB-2M)
(\${MPIRUN} -np ${_NODES} -machinefile \${PBS_NODEFILE} ${_SINK_SMR}) >& RESULT_snksmr
;;
JLAB-3G | JLAB-4G)
(\${MPIRUN} --qmp-f \${QMP_CONF_FILE} --qmp-l \${QMP_LIST_FILE} \${PROGROOT}/${_SINK_SMR}) >& RESULT_snksmr
;;
esac
echo "... done at " \`date\`
echo "files are "
ls -lrt
mv XMLDAT sink_xmldat.\${prop_ctr}.cfg\${cfg}

echo "Finished qqq number \${qqq_ctr} for config \${cfg}"
echo

prop_ctr=\`expr \${prop_ctr} + 1 \`
done    #  end of sink smearing

qqq_ctr=0
#  First check the file exists...
while [ -f qqq_ini.\${qqq_ctr} ]
do

#  barcomps / generalized propagators construction
#
mv qqq_ini.\${qqq_ctr} DATA
echo "Running qqq at " \`date\`
echo "DATA file is "
echo "*********************************************************"
cat DATA
echo "*********************************************************"
echo

sleep 40
echo "Running ${_QQQ} at " \`date\`
case \${SITE} in
JLAB-2M)
(\${MPIRUN} -np ${_NODES} -machinefile \${PBS_NODEFILE} ${_QQQ}) >& RESULT_qqq
;;
JLAB-3G | JLAB-4G)
(\${MPIRUN} --qmp-f \${QMP_CONF_FILE} --qmp-l \${QMP_LIST_FILE} \${PROGROOT}/${_QQQ}) >& RESULT_qqq
;;
esac
echo "... done at " \`date\`
echo "files are "
ls -lrt
mv XMLDAT qqq_\${_UD_PROP_ROOT}_\${_S_PROP_ROOT}.\${qqq_ctr}.cfg\${cfg}

echo "Finished qqq number \${qqq_ctr} for config \${cfg}"
echo

qqq_ctr=\`expr \${qqq_ctr} + 1 \`
done    #  end of qqq_w

# moving files to the appropriate places
#

sleep 40
if [ \${SITE} == JLAB-3G ]; then
if [ \$HOSTNAME == qcd3g000 ]; then
\${RCP} *xmldat*  \${FILESERVER}\${XMLDIR}/.
\${RCP} *PS*.cfg\${cfg}  \${FILESERVER}\${SPCDIR}/.
\${RCP} qqq_*.cfg\${cfg} \${FILESERVER}\${QQQDIR}/.
fi
else
\${RCP} *xmldat*  \${FILESERVER}\${XMLDIR}/.
\${RCP} *PS*.cfg\${cfg}  \${FILESERVER}\${SPCDIR}/.
\${RCP} qqq_*.cfg\${cfg} \${FILESERVER}\${QQQDIR}/.
fi

rm -f *PS*
rm -f qqq_*.cfg*
rm -f RESULT*
rm -f *xmldat*

else
echo "Yikes! You have no \${FILESERVER}\${GAUGEDIR}/\${szingauge}!!"
fi

ctr=\`expr \${ctr} + 1 \`
cfg=\`expr \${cfg} + ${_SKIP_CFG} \`
done    #  end of configuration loop

EOF1

system "chmod u+x ${_PBS_NAME}.sh";
print "Job script is ${_PBS_NAME}.sh\n"


