#!/bin/tcsh

set builddir = ../../scalar/mainprogs/main
set gauge_type = SZINQIO
set gauge_cfg = ../gfix.cfg1

#
# Forward source
#
/bin/rm -f DATA
set version = 5

cat << **EOF** >! DATA
<?xml version="1.0"?>

<make_source>
<annotation>
;
; MAKE_SOURCE input file.
;
; This program is the input file for a  make_source  test run on Wilson-type
; propagators
;
</annotation>

<Param>
 <version>${version}</version>
 <wave_state>S_WAVE</wave_state>
 <source_type>POINT_SOURCE</source_type>
 <j_decay>3</j_decay>
 <t_source>0 0 0 1</t_source>

 <nrow>4 4 4 8</nrow>
</Param>
<Cfg>
 <cfg_type>${gauge_type}</cfg_type>
 <cfg_file>${gauge_cfg}</cfg_file>
</Cfg>
<Prop>
 <source_file>make_source_pt</source_file>
 <source_volfmt>SINGLEFILE</source_volfmt>
</Prop>
</make_source>
**EOF**

# Run the make_source program
$builddir/make_source
/bin/mv XMLDAT make_source_pt_v${version}.xml


#
# Forward propagators
#
/bin/rm -f DATA
set version = 6

cat << **EOF** >! DATA
<?xml version="1.0"?>

<propagator>
<annotation>
; PROPAGATOR input file.
;
; This program is the input file for a propagator (spectroscopy)
; test of the preconditioned Wilson operator.
;
; Use this input file after seqsource
</annotation>

<Param>
 <version>${version}</version>
 <FermTypeP>WILSON</FermTypeP>
 <nonRelProp>false</nonRelProp>
 <FermionAction>
  <FermAct>WILSON</FermAct>
  <Kappa>0.11</Kappa>
 </FermionAction>
 <InvertParam>
   <invType>CG_INVERTER</invType>
   <RsdCG>1.0e-10</RsdCG>
   <MaxCG>1000</MaxCG>
 </InvertParam>
 <nrow>4 4 4 8</nrow>
 <boundary>1 1 1 0</boundary>
</Param>
<Cfg>
 <cfg_type>${gauge_type}</cfg_type>
 <cfg_file>${gauge_cfg}</cfg_file>
</Cfg>
<Prop>
 <source_file>make_source_pt</source_file>
 <prop_file>propagator_pt</prop_file>
 <prop_volfmt>SINGLEFILE</prop_volfmt>
</Prop>
</propagator>
**EOF**

# Run the propagator program
$builddir/propagator
/bin/mv XMLDAT prop_pt_v${version}.xml


#
# Backward source
#
/bin/rm -f DATA
set version = 5

cat << **EOF** >! DATA
<?xml version="1.0"?>

<make_source>
<annotation>
;
; MAKE_SOURCE input file.
;
; This program is the input file for a  make_source  test run on Wilson-type
; propagators
;
</annotation>

<Param>
 <version>${version}</version>
 <wave_state>S_WAVE</wave_state>
 <source_type>WALL_SOURCE</source_type>
 <j_decay>3</j_decay>
 <t_source>0 0 0 6</t_source>

 <nrow>4 4 4 8</nrow>
</Param>
<Cfg>
 <cfg_type>${gauge_type}</cfg_type>
 <cfg_file>${gauge_cfg}</cfg_file>
</Cfg>
<Prop>
 <source_file>make_source_wl</source_file>
 <source_volfmt>SINGLEFILE</source_volfmt>
</Prop>
</make_source>
**EOF**

# Run the make_source program
$builddir/make_source
/bin/mv XMLDAT make_source_wl_v${version}.xml


#
# Forward propagators
#
/bin/rm -f DATA
set version = 6

cat << **EOF** >! DATA
<?xml version="1.0"?>

<propagator>
<annotation>
; PROPAGATOR input file.
;
; This program is the input file for a propagator (spectroscopy)
; test of the preconditioned Wilson operator.
;
; Use this input file after seqsource
</annotation>

<Param>
 <version>${version}</version>
 <FermTypeP>WILSON</FermTypeP>
 <nonRelProp>false</nonRelProp>
 <FermionAction>
  <FermAct>WILSON</FermAct>
  <Kappa>0.11</Kappa>
 </FermionAction>
 <InvertParam>
   <invType>CG_INVERTER</invType>
   <RsdCG>1.0e-10</RsdCG>
   <MaxCG>1000</MaxCG>
 </InvertParam>
 <nrow>4 4 4 8</nrow>
 <boundary>1 1 1 0</boundary>
</Param>
<Cfg>
 <cfg_type>${gauge_type}</cfg_type>
 <cfg_file>${gauge_cfg}</cfg_file>
</Cfg>
<Prop>
 <source_file>make_source_wl</source_file>
 <prop_file>propagator_wl</prop_file>
 <prop_volfmt>SINGLEFILE</prop_volfmt>
</Prop>
</propagator>
**EOF**

# Run the propagator program
$builddir/propagator
/bin/mv XMLDAT prop_wl_v${version}.xml


#
# Forward spectrum
#
/bin/rm -f DATA
set version = 10

cat << **EOF** >! DATA
<?xml version="1.0"?>
<spectrum_w>
<annotation>
; SPECTRUM_W input file.
;
; This program is the input file for a spectrum_w test run on Wilson-type
; propagators
</annotation>

<Param>
 <version>${version}</version>
 <Pt_snk>true</Pt_snk>
 <Sl_snk>false</Sl_snk>
 <Wl_snk>true</Wl_snk>
 <MesonP>true</MesonP>
 <CurrentP>true</CurrentP>
 <BaryonP>true</BaryonP>
 <time_rev>false</time_rev>
 <mom2_max>3</mom2_max>
 <avg_equiv_mom>true</avg_equiv_mom>
 <wvf_kind>GAUGE_INV_GAUSSIAN</wvf_kind>
 <wvf_param>0</wvf_param>
 <wvfIntPar>0</wvfIntPar>
 <nrow>4 4 4 8</nrow>
</Param>
<Cfg>
 <cfg_type>${gauge_type}</cfg_type>
 <cfg_file>${gauge_cfg}</cfg_file>
</Cfg>
<Prop>
 <prop_files>
   <elem>propagator_pt</elem>
 </prop_files>
</Prop>
</spectrum_w>
**EOF**

# Run the spectrum program
$builddir/spectrum_w
/bin/mv XMLDAT spectrum_pt_v${version}.xml


#
# Backward spectrum
#
/bin/rm -f DATA
set version = 10

cat << **EOF** >! DATA
<?xml version="1.0"?>
<spectrum_w>
<annotation>
; SPECTRUM_W input file.
;
; This program is the input file for a spectrum_w test run on Wilson-type
; propagators
</annotation>

<Param>
 <version>${version}</version>
 <Pt_snk>true</Pt_snk>
 <Sl_snk>false</Sl_snk>
 <Wl_snk>true</Wl_snk>
 <MesonP>true</MesonP>
 <CurrentP>true</CurrentP>
 <BaryonP>true</BaryonP>
 <time_rev>false</time_rev>
 <mom2_max>3</mom2_max>
 <avg_equiv_mom>true</avg_equiv_mom>
 <wvf_kind>GAUGE_INV_GAUSSIAN</wvf_kind>
 <wvf_param>0</wvf_param>
 <wvfIntPar>0</wvfIntPar>
 <nrow>4 4 4 8</nrow>
</Param>
<Cfg>
 <cfg_type>${gauge_type}</cfg_type>
 <cfg_file>${gauge_cfg}</cfg_file>
</Cfg>
<Prop>
 <prop_files>
   <elem>propagator_wl</elem>
 </prop_files>
</Prop>
</spectrum_w>
**EOF**

# Run the spectrum program
$builddir/spectrum_w
/bin/mv XMLDAT spectrum_wl_v${version}.xml


#
# Wallformfac
#
/bin/rm -f DATA
set version = 1

cat << **EOF** >! DATA
<?xml version="1.0"?>

<WallFormFac>
<annotation>
; WALLFORMFAC input file.
;
; This program is the input file for a baryon 3pt func. calculation
; test run on Wilson-type propagators
; 
</annotation>

<Param>
 <version>${version}</version>
 <formfac_type>0 1</formfac_type>
 <mom2_max>3</mom2_max>
 <nrow>4 4 4 8</nrow>
</Param>
<Cfg>
 <cfg_type>${gauge_type}</cfg_type>
 <cfg_file>${gauge_cfg}</cfg_file>
</Cfg>
<Prop>
 <forwprop_file>propagator_pt</forwprop_file>
 <backprop_file>propagator_wl</backprop_file>
</Prop>
</WallFormFac>
**EOF**

# Run the wallformfac program
$builddir/wallformfac
/bin/mv XMLDAT wallformfac_v${version}.xml

