$Id: README_buildtest.sh,v 3.0 2006-04-03 04:59:22 edwards Exp $
Original author: Zbigniew Sroczynski

Check the build of the modules QMP, QDP++ and Chroma.

The script takes each specified module and performs the following operations:

1. Update the source from the CVS repository (cvs update)
2. Build the library (gmake)
3. Compile and link the test programs (gmake check)
4. Install the library (gmake install)
5. Compile an link the test programs against the installed module

Each of  modules, has its own directory.
In each of these is 

- a subdirectory which is the checked-out CVS module
e.g. buildtest/qmp/qmp/

- a subdirectory for each type of build
e.g. buildtest/qmp/mpich/

To work automatically, the script uses anonymous pserver access to the CVS 
server.  This works if the user has previously performed a cvs login and the
~/.cvspass file still exists. 
The CVS module should be checked out, because the script uses cvs update.

In the build-type directories is a directory called build, e.g.
buildtest/qmp/mpich/build/
This contains a script, configure.sh, which defines the type of
build. The library is configured and built in this directory.

For each build, the built library is installed in its own subdirectory. These
installation directories are defined and created by configure.sh
At the moment they are in a directory called install/
next to the build/ directory, e.g.  buildtest/qmp/mpich/install/

The results of each operation are logged to a file in subdirectories of the 
logs directory corresponding to the module and build, 
e.g. buildtest/logs/qmp/mpich/
These directories are created by the script if need be.

When an operation succeeds an empty email is sent to a list of addresses.
If any operation fails, the log file is emailed.


The builds to be done are identified by looking for the configure.sh files 
in subdirectories. This behaviour is encapsulated in the function builddirs() 
so you can change it if you want.

The location of the installation directories is set in the configure.sh 
scripts; it has to be known there for the --prefix argument. 
This is set to be an install/ directory with build subdirectories, e.g.
buildtest/qmp/mpich/install/
These directories are created by configure.sh if need be. The location of the 
installation directories is arbitrary, but remember that installation location
of some libraries needs to be known to configure others, e.g. chroma depends 
on qdp++. This is all encapsulated in the configure.sh files if you want to 
change the set-up.
