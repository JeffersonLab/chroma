#! /bin/bash
cp  propagator_comp_zolo5_c0.xml.ini DATA
./propagator_comp
cp XMLDAT XMLDAT_c0

cp propagator_comp_zolo5_c1.xml.ini DATA
./propagator_comp
cp XMLDAT XMLDAT_c1

cp propagator_comp_zolo5_c2.xml.ini DATA
./propagator_comp
cp XMLDAT XMLDAT_c2

# Re use same data file
./collect_propcomp
