#! /bin/bash
cp  multi_propagator_comp_zolo4_c0.xml.ini DATA
./multi_propagator_comp
cp XMLDAT XMLDAT_c0

cp  multi_propagator_comp_zolo4_c1.xml.ini DATA
./multi_propagator_comp
cp XMLDAT XMLDAT_c1

cp multi_propagator_comp_zolo4_c2.xml.ini DATA
./multi_propagator_comp
cp XMLDAT XMLDAT_c2

# Re use same data file
./collect_multi_propcomp
