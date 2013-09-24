#!/bin/tcsh

set gauge_type = WEAK_FIELD
set gauge_cfg = "dummy"

set anisoP = "false"
set Kappa = 0.11
set xi_0 = 1
set nu = 1
set bc = "1 1 1 1"

set N = 4

set nrow = (4 4 4 8)

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

foreach n (1 2)

if ($n == 1) then
  set rnd = (201 213 215 217)
endif

if ($n == 2) then
  set rnd = (113 115 117 119)
endif

set t = 0
while ($t < $nrow[4])
cat << **EOF** 
    <elem>
      <annotation>
         Diluted complex Z(2) = Z(4) random source. Even sites.
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

          <spatial_mask_size>2 2 2</spatial_mask_size>
          <spatial_mask>
             <elem>0 0 0</elem>
             <elem>1 1 0</elem>
             <elem>1 0 1</elem>
             <elem>0 1 1</elem>
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
        <file_name>./zN_prop_q${n}_e_t${t}</file_name>
        <file_volfmt>MULTIFILE</file_volfmt>
      </File>
    </elem>


    <elem>
      <annotation>
         Diluted complex Z(2) = Z(4) random source. Odd sites.
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

          <spatial_mask_size>2 2 2</spatial_mask_size>
          <spatial_mask>
             <elem>1 0 0</elem>
             <elem>0 1 0</elem>
             <elem>0 0 1</elem>
             <elem>1 1 1</elem>
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
        <file_name>./zN_prop_q${n}_o_t${t}</file_name>
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
      ; STOCH_MESON input file.
      </annotation>
      
      <Name>STOCH_MESON</Name>
      <Frequency>1</Frequency>
      <Param> 
        <version>1</version>
        <mom2_max>0</mom2_max>
      </Param>
      <SourceSmearing>
        <version>6</version>
        <Source>
          <version>2</version>
          <SourceType>POINT_SOURCE</SourceType>
          <j_decay>3</j_decay>

          <Displacement>
            <version>1</version>
            <DisplacementType>NONE</DisplacementType>
          </Displacement>
        </Source>
      </SourceSmearing>

      <SinkSmearing>
        <version>5</version>
        <Sink>
          <version>2</version>
          <SinkType>POINT_SINK</SinkType>
          <j_decay>3</j_decay>

          <Displacement>
            <version>1</version>
            <DisplacementType>NONE</DisplacementType>
          </Displacement>
        </Sink>
      </SinkSmearing>

      <annotation>
      <SinkSmearing>
        <version>5</version>
              <Sink>
                <version>2</version>
                <SinkType>SHELL_SINK</SinkType>
                <j_decay>3</j_decay>

                <SmearingParam>
                  <wvf_kind>GAUGE_INV_GAUSSIAN</wvf_kind>
                  <wvf_param>2.0</wvf_param>
                  <wvfIntPar>5</wvfIntPar>
                  <no_smear_dir>3</no_smear_dir>
                </SmearingParam>

                <Displacement>
                  <version>1</version>
                  <DisplacementType>NONE</DisplacementType>
                </Displacement>

                <LinkSmearing>
                  <LinkSmearingType>APE_SMEAR</LinkSmearingType>
                  <link_smear_fact>2.5</link_smear_fact>
                  <link_smear_num>0</link_smear_num>
                  <no_smear_dir>3</no_smear_dir>
                </LinkSmearing>
              </Sink>
      </SinkSmearing>
      </annotation>

      <NamedObject>
        <gauge_id>default_gauge_field</gauge_id>
        <Prop>
          <operator_file>meson_operator.dat</operator_file>
          <operator>
            <elem>
              <soln_files>
**EOF**

set t = 0
while ($t < $nrow[4])
cat << **EOF** 
                <elem>./zN_prop_q1_e_t${t}</elem>
                <elem>./zN_prop_q1_o_t${t}</elem>
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
                <elem>./zN_prop_q2_e_t${t}</elem>
                <elem>./zN_prop_q2_o_t${t}</elem>
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



