#!/bin/tcsh

set what = $1

(cd parscalar-gigE ; do_config ; make ; cd mainprogs/main; make $what)
(cd parscalar-gigE-double ; do_config ; make ; cd mainprogs/main; make $what)
