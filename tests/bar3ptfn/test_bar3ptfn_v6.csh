#!/bin/tcsh

set builddir = ../../scalar/mainprogs/main
set gauge_cfg = ../gfix.cfg1

#
# Forward source
#
/bin/rm -f DATA

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
 <version>5</version>
 <wave_state>S_WAVE</wave_state>
 <source_type>POINT_SOURCE</source_type>
 <j_decay>3</j_decay>
 <t_source>0 0 0 1</t_source>

 <nrow>4 4 4 8</nrow>
</Param>
<Cfg>
 <cfg_type>SZINQIO</cfg_type>
 <cfg_file>${gauge_cfg}</cfg_file>
</Cfg>
<Prop>
 <source_file>make_source_0</source_file>
 <source_volfmt>SINGLEFILE</source_volfmt>
</Prop>
</make_source>
**EOF**

# Run the make_source program
$builddir/make_source
/bin/mv XMLDAT make_source.xml


#
# Forward propagators
#
/bin/rm -f DATA

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
 <version>6</version>
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
 <cfg_type>SZINQIO</cfg_type>
 <cfg_file>${gauge_cfg}</cfg_file>
</Cfg>
<Prop>
 <source_file>make_source_0</source_file>
 <prop_file>propagator_0</prop_file>
 <prop_volfmt>SINGLEFILE</prop_volfmt>
</Prop>
</propagator>
**EOF**

# Run the propagator program
$builddir/propagator
/bin/mv XMLDAT prop.xml


#
# Sequential propagators
#
foreach i (NUCL_U_UNPOL NUCL_D_UNPOL NUCL_U_POL NUCL_D_POL PION)

/bin/rm -f DATA

cat << **EOF** >! DATA
<?xml version="1.0"?>

<seqsource>
<annotation>
; SEQSOURCE input file.
;
; This program is the input file for a seqsource
; test of the preconditioned Wilson operator.
</annotation>

<Param>
 <version>1</version>
 <seq_src>$i</seq_src>
 <t_sink>6</t_sink>
 <sink_mom>0 0 0</sink_mom>
 <nrow>4 4 4 8</nrow>
</Param>

<PropSink>
 <version>4</version>
 <wave_state>S_WAVE</wave_state>
 <sink_type>WALL_SINK</sink_type>

 <nrow>4 4 4 8</nrow>
</PropSink>

<Cfg>
 <cfg_type>SZINQIO</cfg_type>
 <cfg_file>${gauge_cfg}</cfg_file>
</Cfg>
<Prop>
 <prop_file>propagator_0</prop_file>
 <seqsource_file>seqsource_0_$i</seqsource_file>
 <seqsource_volfmt>SINGLEFILE</seqsource_volfmt>
</Prop>
</seqsource>
**EOF**

# Run the seqsource program
$builddir/seqsource
/bin/mv XMLDAT seqsource.$i.xml



/bin/rm -f DATA

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
 <version>6</version>
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
 <cfg_type>SZINQIO</cfg_type>
 <cfg_file>${gauge_cfg}</cfg_file>
</Cfg>
<Prop>
 <source_file>seqsource_0_$i</source_file>
 <prop_file>seqprop_0_$i</prop_file>
 <prop_volfmt>SINGLEFILE</prop_volfmt>
</Prop>
</propagator>
**EOF**

# Run the propagator program
$builddir/propagator
/bin/mv XMLDAT prop_seqsource_v6.$i.xml

end


#
# Bar3ptfn
#
/bin/rm -f DATA

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
 <version>10</version>
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
 <cfg_type>SZINQIO</cfg_type>
 <cfg_file>${gauge_cfg}</cfg_file>
</Cfg>
<Prop>
 <prop_files>
   <elem>propagator_0</elem>
 </prop_files>
</Prop>
</spectrum_w>
**EOF**

# Run the spectrum program
$builddir/spectrum_w
/bin/mv XMLDAT spectrum.xml


#
# Bar3ptfn
#
/bin/rm -f DATA

cat << **EOF** >! DATA
<?xml version="1.0"?>

<bar3ptfn>
<annotation>
;
; BAR3PTFN input file.
;
; This program is the input file for a baryon 3pt func. calculation
; test run on Wilson-type propagators
;
; Run the   test_prop_seqsource_vXXX.csh  script in  tests/seqsource
</annotation>

<Param>
 <version>6</version>
 <nrow>4 4 4 8</nrow>
 <mom2_max>3</mom2_max>
</Param>
<Cfg>
 <cfg_type>SZINQIO</cfg_type>
 <cfg_file>${gauge_cfg}</cfg_file>
</Cfg>
<Prop>
 <prop_file>propagator_0</prop_file>
 <seqprop_files>
   <elem>seqprop_0_NUCL_U_UNPOL</elem>
   <elem>seqprop_0_NUCL_D_UNPOL</elem>
   <elem>seqprop_0_NUCL_U_POL</elem>
   <elem>seqprop_0_NUCL_D_POL</elem>
   <elem>seqprop_0_PION</elem>
 </seqprop_files>
</Prop>
</bar3ptfn>
**EOF**

# Run the bar3ptfn program
$builddir/bar3ptfn
/bin/mv XMLDAT bar3ptfn.xml

