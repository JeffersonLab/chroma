#!/bin/tcsh

set what = $1

./build2m.csh $what >&! output.2m
./build3g.csh $what >&! output.3g
