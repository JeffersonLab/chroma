#!/bin/tcsh

set what = $1

(cd scalar-intel ; do_config ; make ; cd mainprogs/main; make $what)
(cd parscalar-gm-intel ; do_config ; make ; cd mainprogs/main; make $what)
