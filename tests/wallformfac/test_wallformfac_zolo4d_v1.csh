#!/bin/tcsh

set builddir = ../../build/double/mainprogs/main
#set gauge_type = UNIT
set gauge_type = SZINQIO
set gauge_cfg = ../gfix.cfg1
set inv_type = REL_GMRESR_CG_INVERTER
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
 <t_source>0 0 0 0</t_source>

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
  <FermAct>ZOLOTAREV_4D</FermAct>
  <Mass>0.04</Mass>
  <RatPolyDeg>17</RatPolyDeg>
  <RatPolyDegPrecond>5</RatPolyDegPrecond>
  <InnerSolve>
    <MaxCG>1000</MaxCG>
    <RsdCG>1.0e-7</RsdCG>
    <ReorthFreq>20</ReorthFreq>
  </InnerSolve> 
  <AuxFermAct>
    <FermAct>WILSON</FermAct>
    <Mass>-1.4</Mass>
  </AuxFermAct>
  <StateInfo>
    <!-- The ApproxMin and ApproxMax are ignored if NWilsVec &gt 0 -->
    <ApproxMin>0.66</ApproxMin>
    <ApproxMax>6.5</ApproxMax>
    <NWilsVec>0</NWilsVec>
    <Eig>
        <eigen_file_stem>../t_ritz/g5_wils_eigen</eigen_file_stem>
        <eigen_volfmt>SINGLEFILE</eigen_volfmt>
    </Eig>
  </StateInfo>
 </FermionAction>

 <InvertParam>
   <invType>${inv_type}</invType>
   <RsdCG>1.0e-6</RsdCG>
   <RsdCGPrec>0.1</RsdCGPrec>
   <MaxCG>1000</MaxCG>
   <MaxCGPrec>20</MaxCGPrec>
 </InvertParam>
 <nrow>4 4 4 8</nrow>
 <boundary>1 1 1 -1</boundary>
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
 <t_source>0 0 0 4</t_source>

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
  <FermAct>ZOLOTAREV_4D</FermAct>
  <Mass>0.04</Mass>
  <RatPolyDeg>17</RatPolyDeg>
  <InnerSolve>
    <MaxCG>1000</MaxCG>
    <RsdCG>1.0e-7</RsdCG>
    <ReorthFreq>20</ReorthFreq>
  </InnerSolve>
  <AuxFermAct>
    <FermAct>WILSON</FermAct>
    <Mass>-1.4</Mass>
  </AuxFermAct>
  <StateInfo>
    <!-- The ApproxMin and ApproxMax are ignored if NWilsVec &gt 0 -->
    <ApproxMin>0.66</ApproxMin>
    <ApproxMax>6.5</ApproxMax>
    <NWilsVec>0</NWilsVec>
    <Eig>
        <eigen_file_stem>../t_ritz/g5_wils_eigen</eigen_file_stem>
        <eigen_volfmt>SINGLEFILE</eigen_volfmt>
    </Eig>
  </StateInfo>
 </FermionAction>
 <InvertParam>
   <invType>${inv_type}</invType>
   <RsdCG>1.0e-6</RsdCG>
   <RsdCGPrec>0.1</RsdCGPrec>
   <MaxCG>1000</MaxCG>
   <MaxCGPrec>12</MaxCGPrec>
 </InvertParam>
 <nrow>4 4 4 8</nrow>
 <boundary>1 1 1 -1</boundary>
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
 <formfac_type>0 1 2</formfac_type>
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

