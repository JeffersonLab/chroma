#!/bin/bash
# $Id: build.sh,v 1.3 2005-12-20 19:32:36 edwards Exp $
#
#  Original author: Zbigniew Sroczynski
#  See README_buildtest.sh for more information.
#

cd `dirname $0`

checkdir=$PWD

logdir=$checkdir/logs

modules="qmp qdp++ chroma"


# Send mail to these addresses
 
mailto="edwards@jlab.org bjoo@jlab.org"

# Send mail if something goes wrong

failmail(){
    cat $logfile | mail -s "buildtest: $module $build $1 failed" $mailto &> /dev/null
}

# Send mail if something goes right

successmail(){
    mail -s "buildtest: $module $build $1 succeeded" $mailto &> /dev/null <<EOF
$logfile
EOF
}

# Identify the subdirectories for the different types of build

builddirs(){
    ls $moduledir/*/build/configure.sh | sed -e s,$moduledir/,, -e s,/build/configure.sh,,
}


perform_action(){

    logfile=$logdir/$module/$build/${action_name}_$now
    uname -a > $logfile
    echo '------------------------' >> $logfile

    $action >> $logfile 2>&1
    local status=$?

    if [ $status -ne 0 ]
    then
    	failmail $action_name 
    else
    	successmail $action_name 
    fi 

    gzip $logfile
    return $status 
}

for module in $modules
do

    moduledir=$checkdir/$module

    # Check out new code

    srcdir=$moduledir/$module

    cd $srcdir

    now=`date +"%F_%H.%M_%Z"`
    if ! [ -d $logdir/${module} ] 
    then	
	mkdir -p $logdir/${module}
    fi

    action_name="update"
    # Assume that the module is already checked out and the .cvspass file 
    # already exists	
    action="cvs -d :pserver:anonymous@cvs.jlab.org:/group/lattice/cvsroot update -d -C"
    build=""
    perform_action
    [ $? -eq 0 ] || continue
    	

    # Build the code

    for build in `builddirs`
    do

        # Configure

        builddir=$moduledir/$build/build
        cd $builddir

        if ! [ -d $logdir/$module/$build ] 
	then	
	    mkdir -p $logdir/$module/$build
        fi

	action_name="configure"
	action="sh configure.sh"
	perform_action
	[ $? -eq 0 ] || continue

	# Build 

	action_name="build"
	action="gmake -k"
	perform_action
	[ $? -eq 0 ] || continue

	# Check

	action_name="check"
	action="gmake -k check"
	perform_action
	[ $? -eq 0 ] || continue
	
	# XCheck

        if [ $module == "chroma" ]
        then
	        action_name="xcheck"
                action="gmake -k xcheck"
        	perform_action
        fi 
	
	# Install

	action_name="install"
	action="gmake install"
	perform_action
	[ $? -eq 0 ] || continue

	# Link
	
	[ ! -d ../link ] && mkdir ../link
        cd ../link

	action_name="link"
	action="perl $checkdir/link $module"
	perform_action
	[ $? -eq 0 ] || continue

	cd ..
	rm -r link

    done	# loop over builds
		
done   # loop over modules


# Tidy up

for module in $modules
do
    moduledir=$checkdir/$module
    for build in `builddirs`
    do
       cd $moduledir/$build/build
       gmake -k uninstall  distclean &> /dev/null	
    done
done

# Remove old logfiles

find $logdir -type f -mtime +14 -exec rm '{}' ';'
