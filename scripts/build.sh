#!/bin/bash
# $Id: build.sh,v 3.0 2006-04-03 04:59:22 edwards Exp $
#
#  Original author: Zbigniew Sroczynski
#  See README_buildtest.sh for more information.
#

cd `dirname $0`

checkdir=$PWD

logdir=$checkdir/logs

modules="bagel_qdp qmp qdp++ chroma adat"

failcnt=0

# Send mail to these addresses
# In case of failure mail everyone on this list  
failmailto="edwards@jlab.org bjoo@jlab.org mcneile@jlab.org"

# In case of success send to archive list only
successmailto=""

# Send mail if something goes wrong

failmail(){
    cat $logfile | mail -s "buildtest: $module $build $1 failed" $failmailto &> /dev/null
}

# Send mail if something goes right
# NOTE: had to comment this out since the list is null. Should have a
# check instead.

successmail(){
#    mail -s "buildtest: $module $build $1 succeeded" $successmailto &> /dev/null <<EOF
#$logfile
#EOF
echo "do not send out successes" > /dev/null
}

# Send mail to indicate the tests are all done -- send to everyone
finishmail(){
    mail -s "buildtest: All Tests Complete: $failcnt errors" $failmailto &> /dev/null <<EOF
All tests done: $failcnt errors
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
        let failcnt++
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
	    # Skip regressions on IB - dont know how to do the running
	    # Need to generalise system to be able to launch mpi jobs
	    # As regressions
	    if [ $build != "parscalar-ib" ]
            then
	        action_name="xcheck"
                action="gmake -k xcheck"
        	perform_action
            fi
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
finishmail
