#!/bin/sh

# submits a single job to the CMU cluster with the given parameters

if [ $# -ne 5 ]; then
    echo "Usage: $0 <Gauss rad> <Gauss itr> <Stout rho> <Stout itr> <queue>"
    exit 1
fi

# specify location of this script and create_pbs-template.pl
topDir="/home/username/research/scripts"  # CHANGE THIS

# specify the location of run_qqq_spectrum.sh and the output files
runDir="/home/username/research/runs"  # CHANGE THIS

# specify the place to put the qqq propagator files
outputdir_root="/raid3/username/NFO/aniso/6p1_12_48_xi3p0_wl/12^3x48"  # CHANGE THIS

template="${topDir}/create_pbs.pl-template"
run_script="${runDir}/create_pbs.pl"

gauss_rad_param=$1
gauss_itr_param=$2
stout_rho_param=$3
stout_itr_param=$4
queue_name=$5
proproot_param="QG"`echo $gauss_rad_param | tr . p`"_${gauss_itr_param}_LS"`echo $stout_rho_param | tr . p`"_${stout_itr_param}"
qqq_dir="${outputdir_root}/qqq/qqq-${proproot_param}"

if [ ! -d $qqq_dir ]; then
    echo "Creating $qqq_dir"
    mkdir $qqq_dir
fi

sed -e "s/__QUEUE_NAME__/$queue_name/g; \
        s/__GAUSS_RAD_PARAM__/$gauss_rad_param/g; \
        s/__GAUSS_ITR_PARAM__/$gauss_itr_param/g; \
        s/__STOUT_RHO_PARAM__/$stout_rho_param/g; \
        s/__STOUT_ITR_PARAM__/$stout_itr_param/g; \
        s/__PROPROOT_PARAM__/$proproot_param/g" < $template > $run_script

firstConfig=10
lastConfig=59

configNum=$firstConfig
while [ $configNum -le $lastConfig ]; do 
    ssh -x enrico "cd ${runDir}; ./create_pbs.pl run_qqq_spectrum.sh $configNum 1 5 QLsmr"
    configNum=`expr $configNum + 5`
done

exit 0
