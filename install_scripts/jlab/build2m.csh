#!/bin/tcsh

set what = $1

(cd scalar ; do_config ; make ; cd mainprogs/main; make $what)
(cd parscalar-gm ; do_config ; make ; cd mainprogs/main; make $what)
(cd scalar-double ; do_config ; make ; cd mainprogs/main; make $what)
(cd parscalar-gm-double ; do_config ; make ; cd mainprogs/main; make $what)
