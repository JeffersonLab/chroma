#!/bin/tcsh

set builddir = .

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
 <t_sink>5</t_sink>
 <sink_mom>0 0 0</sink_mom>
 <nrow>4 4 4 8</nrow>
</Param>

<PropSink>
 <version>4</version>
 <wave_state>S_WAVE</wave_state>
 <sink_type>SHELL_SINK</sink_type>

 <ShellSink>
   <SinkSmearingParam>
     <wvf_kind>GAUGE_INV_GAUSSIAN</wvf_kind>
     <wvf_param>2.0</wvf_param>
     <wvfIntPar>5</wvfIntPar>
   </SinkSmearingParam>
   <laplace_power>0</laplace_power>
   <link_smear_fact>2.5</link_smear_fact>
   <link_smear_num>0</link_smear_num>
   <disp_length>0</disp_length>
   <disp_dir>0</disp_dir>
 </ShellSink>

 <nrow>4 4 4 8</nrow>
</PropSink>

<Cfg>
 <cfg_type>SZIN</cfg_type>
 <cfg_file>../test_purgaug.cfg1</cfg_file>
</Cfg>
<Prop>
 <prop_file>../propagator_comp/propagator_0</prop_file>
 <seqsource_file>seqsource_0_$i</seqsource_file>
 <seqsource_volfmt>SINGLEFILE</seqsource_volfmt>
</Prop>
</seqsource>
**EOF**

# Run the seqsource program
$builddir/seqsource
/bin/mv XMLDAT seqsource_v1.$i.xml



/bin/rm -f DATA

cat << **EOF** >! DATA
<?xml version="1.0"?>

<propagatorComp>
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
  <AnisoParam>
    <anisoP>false</anisoP>
    <t_dir>3</t_dir>
    <xi_0>1</xi_0>
    <nu>1</nu>
  </AnisoParam>
 </FermionAction>
 <InvertParam>
   <invType>CG_INVERTER</invType>
   <RsdCG>1.0e-10</RsdCG>
   <MaxCG>1000</MaxCG>
 </InvertParam>
 <nrow>4 4 4 8</nrow>
 <boundary>1 1 1 -1</boundary>
 <t_source>0 0 0 0</t_source>
</Param>
<Cfg>
 <cfg_type>SZIN</cfg_type>
 <cfg_file>../test_purgaug.cfg1</cfg_file>
</Cfg>
<Prop>
 <source_file>../seqsource/seqsource_0_$i</source_file>
 <prop_file>../seqsource/seqprop_0_$i</prop_file>
 <prop_volfmt>SINGLEFILE</prop_volfmt>
</Prop>
<Components>
    <elem><spin>0</spin><color>0</color></elem>
    <elem><spin>1</spin><color>0</color></elem>
    <elem><spin>2</spin><color>0</color></elem>
    <elem><spin>3</spin><color>0</color></elem>
</Components>
</propagatorComp>
**EOF**

# Run the propagator program
./propagator_comp

/bin/rm -f DATA

cat << **EOF** >! DATA
<?xml version="1.0"?>

<propagatorComp>
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
  <AnisoParam>
    <anisoP>false</anisoP>
    <t_dir>3</t_dir>
    <xi_0>1</xi_0>
    <nu>1</nu>
  </AnisoParam>
 </FermionAction>
 <InvertParam>
   <invType>CG_INVERTER</invType>
   <RsdCG>1.0e-10</RsdCG>
   <MaxCG>1000</MaxCG>
 </InvertParam>
 <nrow>4 4 4 8</nrow>
 <boundary>1 1 1 -1</boundary>
 <t_source>0 0 0 0</t_source>
</Param>
<Cfg>
 <cfg_type>SZIN</cfg_type>
 <cfg_file>../test_purgaug.cfg1</cfg_file>
</Cfg>
<Prop>
 <source_file>../seqsource/seqsource_0_$i</source_file>
 <prop_file>../seqsource/seqprop_0_$i</prop_file>
 <prop_volfmt>SINGLEFILE</prop_volfmt>
</Prop>
<Components>
    <elem><spin>0</spin><color>1</color></elem>
    <elem><spin>1</spin><color>1</color></elem>
    <elem><spin>2</spin><color>1</color></elem>
    <elem><spin>3</spin><color>1</color></elem>
</Components>
</propagatorComp>
**EOF**

# Run the propagator program
./propagator_comp

/bin/rm -f DATA

cat << **EOF** >! DATA
<?xml version="1.0"?>

<propagatorComp>
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
  <AnisoParam>
    <anisoP>false</anisoP>
    <t_dir>3</t_dir>
    <xi_0>1</xi_0>
    <nu>1</nu>
  </AnisoParam>
 </FermionAction>
 <InvertParam>
   <invType>CG_INVERTER</invType>
   <RsdCG>1.0e-10</RsdCG>
   <MaxCG>1000</MaxCG>
 </InvertParam>
 <nrow>4 4 4 8</nrow>
 <boundary>1 1 1 -1</boundary>
 <t_source>0 0 0 0</t_source>
</Param>
<Cfg>
 <cfg_type>SZIN</cfg_type>
 <cfg_file>../test_purgaug.cfg1</cfg_file>
</Cfg>
<Prop>
 <source_file>../seqsource/seqsource_0_$i</source_file>
 <prop_file>../seqsource/seqprop_0_$i</prop_file>
 <prop_volfmt>SINGLEFILE</prop_volfmt>
</Prop>
<Components>
    <elem><spin>0</spin><color>2</color></elem>
    <elem><spin>1</spin><color>2</color></elem>
    <elem><spin>2</spin><color>2</color></elem>
    <elem><spin>3</spin><color>2</color></elem>
</Components>
</propagatorComp>
**EOF**

# Run the propagator program
./propagator_comp


./collect_propcomp

mv XMLDAT prop_comp_seqsource_comp_v6.$i.xml

end

