#!/bin/tcsh

set gauge_type = WEAK_FIELD
set gauge_cfg = "dummy"

set anisoP = "false"
set Kappa = 0.11
set xi_0 = 1
set nu = 1
set bc = "1 1 1 1"

set N = 4

set nrow = (2 2 2 4)
#####set nrow = (4 4 4 8)

set RsdCG = 1.0e-10
#set RsdCG = 1.0e-1
set MaxCG = 2000

#
# Build chroma input
#
#/bin/rm -f  DATA

cat << **EOF** 
<?xml version="1.0"?>
<chroma>
<annotation>
;
; 
;
</annotation>
<Param> 
  <InlineMeasurements>

**EOF**

foreach n (1 2 3)

if ($n == 1) then
  set rnd = (471 1694 3965 563)
endif

if ($n == 2) then
  set rnd = (714 1573 3042 517)
endif

if ($n == 3) then
  set rnd = (2781 1884 3019 716)
endif

set t = 0
while ($t < $nrow[4])
cat << **EOF** 
    <elem>
      <annotation>
         Diluted complex Z(2) = Z(4) random source.
      </annotation>
      <Name>MAKE_SOURCE_FERM</Name>
      <Frequency>1</Frequency>
      <Param>
        <version>6</version>
        <Source>
          <SourceType>RAND_DILUTE_ZN_SOURCE</SourceType>
          <version>1</version>
          <N>${N}</N>
          <j_decay>3</j_decay>
          <t_source>${t}</t_source>
          <ran_seed>
            <Seed>	
              <elem>$rnd[1]</elem>
              <elem>$rnd[2]</elem>
              <elem>$rnd[3]</elem>
              <elem>$rnd[4]</elem>
            </Seed>
          </ran_seed>

          <spatial_mask_size>1 1 1</spatial_mask_size>
          <spatial_mask>
             <elem>0 0 0</elem>
          </spatial_mask>

          <color_mask>0 1 2</color_mask>
          <spin_mask>0 1 2 3</spin_mask>
        </Source>
      </Param>
      <NamedObject>
        <gauge_id>default_gauge_field</gauge_id>
        <source_id>zN_source</source_id>
      </NamedObject>
    </elem>

    <elem>
      <Name>PROPAGATOR_FERM</Name>
      <Frequency>1</Frequency>
      <Param>
        <version>10</version>
        <quarkSpinType>FULL</quarkSpinType>
        <obsvP>false</obsvP>
        <numRetries>1</numRetries>
        <FermionAction>
         <FermAct>WILSON</FermAct>
         <Kappa>0.11</Kappa>
         <AnisoParam>
           <anisoP>${anisoP}</anisoP>
           <t_dir>3</t_dir>
           <xi_0>${xi_0}</xi_0>
           <nu>${nu}</nu>
         </AnisoParam>
         <FermionBC>
           <FermBC>SIMPLE_FERMBC</FermBC>
           <boundary>${bc}</boundary>
         </FermionBC>
        </FermionAction>
        <InvertParam>
          <invType>CG_INVERTER</invType>
          <RsdCG>${RsdCG}</RsdCG>
          <MaxCG>${MaxCG}</MaxCG>
        </InvertParam>
      </Param>
      <NamedObject>
        <gauge_id>default_gauge_field</gauge_id>
        <source_id>zN_source</source_id>
        <prop_id>zN_prop</prop_id>
      </NamedObject>
    </elem>

    <elem>
      <Name>ERASE_NAMED_OBJECT</Name>
      <Frequency>1</Frequency>
      <NamedObject>
        <object_id>zN_source</object_id>
      </NamedObject>
    </elem>

    <elem>
      <Name>QIO_WRITE_ERASE_NAMED_OBJECT</Name>
      <Frequency>1</Frequency>
      <NamedObject>
        <object_id>zN_prop</object_id>
        <object_type>LatticeFermion</object_type>
      </NamedObject>
      <File>
        <file_name>./zN_prop_q${n}_t${t}</file_name>
        <file_volfmt>MULTIFILE</file_volfmt>
      </File>
    </elem>
**EOF**

@ t ++

end  # while(t < Lt)

end  # foreach n

# More stuff
cat << **EOF** 

    <elem>
      <annotation>
      ; STOCH_BARYON input file.
      </annotation>
      
      <Name>STOCH_BARYON</Name>
      <Frequency>1</Frequency>
      <Param> 
        <version>2</version>
        <mom2_max>0</mom2_max>

        <BaryonOperator>
          <version>2</version>
          <BaryonOperatorType>NUCLEON</BaryonOperatorType>

          <!-- BaryonOperatorType>GROUP_BARYON</BaryonOperatorType -->
          <!-- operator_coeff_file>./Single_Site</operator_coeff_file -->
          <!-- displacement_length>1</displacement_length -->

          <SourceQuarkSmearing>
            <wvf_kind>NONE</wvf_kind>
          </SourceQuarkSmearing>

          <SinkQuarkSmearing>
            <wvf_kind>NONE</wvf_kind>

            <!-- wvf_kind>GAUGE_INV_GAUSSIAN</wvf_kind -->
            <wvf_param>2.0</wvf_param>
            <wvfIntPar>5</wvfIntPar>
            <no_smear_dir>3</no_smear_dir>
          </SinkQuarkSmearing>

          <LinkSmearing>
            <LinkSmearingType>NONE</LinkSmearingType>

            <!-- LinkSmearingType>STOUT_SMEAR</LinkSmearingType -->
            <link_smear_fact>0.1625</link_smear_fact>
            <link_smear_num>4</link_smear_num>
            <no_smear_dir>3</no_smear_dir>
          </LinkSmearing>
        </BaryonOperator>
      </Param>
      <NamedObject>
        <gauge_id>default_gauge_field</gauge_id>
        <Prop>
          <operator_file>baryon_operator.dat</operator_file>
          <operator>
            <elem>
              <soln_files>
**EOF**

set t = 0
while ($t < $nrow[4])
cat << **EOF** 
                <elem>./zN_prop_q1_t${t}</elem>
**EOF**

@ t ++

end  # while(t < Lt)

cat << **EOF** 
              </soln_files>
            </elem>
            <elem>
              <soln_files>
**EOF**

set t = 0
while ($t < $nrow[4])
cat << **EOF** 
                <elem>./zN_prop_q2_t${t}</elem>
**EOF**

@ t ++

end  # while(t < Lt)

cat << **EOF** 
              </soln_files>
            </elem>
            <elem>
              <soln_files>
**EOF**

set t = 0
while ($t < $nrow[4])
cat << **EOF** 
                <elem>./zN_prop_q3_t${t}</elem>
**EOF**

@ t ++

end  # while(t < Lt)

cat << **EOF** 
              </soln_files>
            </elem>
          </operator>
        </Prop>
      </NamedObject>
    </elem>

  </InlineMeasurements>
  <nrow>${nrow}</nrow>
</Param>
<Cfg>
  <cfg_type>${gauge_type}</cfg_type>
  <cfg_file>${gauge_cfg}</cfg_file>
</Cfg>
</chroma>
**EOF**



