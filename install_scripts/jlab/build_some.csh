#!/bin/tcsh

set what = $1

(cd parscalar-gigE ; do_config ; make ; make $what)
(cd parscalar-gm ; do_config ; make ; make $what)
(cd scalar ; do_config ; make ; make $what)

