#!/usr/bin/perl
#
# $id:$
#
#  This creates a job script given a template job script
#
#  There are various parameters that are set below.  The template must be
#  written so that these parameters are passed to the template script using
#  the appropriate substitions as written below.
#


die "Usage: create_pbs <pbs template> <first config> <config_skip> <nconfig> <pbs job>" unless $#ARGV eq 4;

#
#  Check template file exists

die "File $ARGV[0] does not exist\n" unless -f $ARGV[0];

$pbs_template = $ARGV[0];
$first_cfg = $ARGV[1];
$skip = $ARGV[2];
$ncfg = $ARGV[3];
$pbs_name = $ARGV[4];

$full_name = $pbs_name."_".${first_cfg}."_".$skip."_".$ncfg;

###################################################################
#
#  Here we have various parameters that change between runs
#
###################################################################

#
#  Job stuff
#$nodes=2;			# Number of nodes
#$attributes=":MYRINET";		# Any attributes

$nodes=8;			# Number of nodes
$attributes=":BIGMEM";		# Any attributes



#
#  Lattice sizes and ensembles
$lsize=12;			# Spatial lattice size
$tsize=48;			# Temporal lattice size
$cfglabel="NF0/aniso/6p1_12_48_xi3p0_wl/"; # Ensemble directory
$gaugeroot="wlq_6p1_12_48_xi3p0.cfg"; # Root name of configs
$qqq_props="/home/dgr/simulation_inputs/qqq_nucleon"; # input file

#
#  Kappa values etc
$mass=-0.28;			# mass
$s_mass=-0.28;
$ud_mass=-0.28;
$nu=0.915;			# Fermion anistropy renormaliation
$xi_0=2.464;			# Bare anisotropy

#
#  Smearing values
$gauss_rad=3.6;
$gauss_itr=32;

#
#  Open the template file
open(TEMPLATE, "< $pbs_template");
open(PBS_OUT, "> $full_name.sh");

#
#  Now set the parameters in the output file
while (<TEMPLATE>){
    s/__NODES__/$nodes/;
    s/__ATTRIBUTE__/$attributes/;
    s/__CONFIG_START__/$first_cfg/;
    s/__NCONFIG__/$ncfg/;
    s/__CONFIG_SKIP__/$skip/;
    s/__JOBNAME__/$full_name/;
    s/__LSIZE__/$lsize/;
    s/__TSIZE__/$tsize/;
    s/__GAUGEROOT__/$gaugeroot/;
    s/__CFGLABEL__/$cfglabel/;
    s/__MASS__/$mass/;
    s/__UD_MASS__/$ud_mass/;
    s/__S_MASS__/$s_mass/;
    s/__XI_0__/$xi_0/;
    s/__NU__/$nu/;
    s/__GAUSS_RAD__/$gauss_rad/;
    s/__GAUSS_ITR__/$gauss_itr/;
    s/__QQQ_PROPS__/$qqq_props/;
    print PBS_OUT $_;
}
close(TEMPLATE);
close(PBS_OUT);

#
# Submit the job
#system("qsub $full_name.sh");
