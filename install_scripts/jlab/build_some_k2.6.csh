#!/bin/tcsh

set what = $1

(cd parscalar-gigE-k2.6 ; do_config ; make ; make $what)
#(cd parscalar-gm ; do_config ; make ; make $what)
(cd scalar ; do_config ; make ; make $what)
(cd parscalar-single ; do_config ; make ; make $what)
