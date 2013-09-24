#!/usr/bin/perl
#
# Used to generate the input file  stoch_group_meson.ini.xml
#

$gauge_type = "WEAK_FIELD";
$gauge_cfg = "dummy";

$anisoP = "false";
$mass = 0.0;
$xi_0 = 1;
$nu = 1;
@bc = (1, 1, 1, -1);

$N = 4;

@nrow = (2, 2, 2, 4);

$RsdCG = 1.0e-10;
$MaxCG = 2000;

$Lt = $nrow[3];

#
# Build chroma input
#

print <<EOF;
<?xml version="1.0"?>
<chroma>
<annotation>
;
; 
;
</annotation>
<Param> 
  <InlineMeasurements>

EOF

foreach $n (1, 2)
{
  if ($n == 1) {@rnd = (471, 1694, 3965, 563);}
  if ($n == 2) {@rnd = (714, 1573, 3042, 517);}

#  foreach $t (0 .. $Lt-1)
  foreach $t (0 .. 0)
  {
print <<EOF;
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
              <elem>$rnd[0]</elem>
              <elem>$rnd[1]</elem>
              <elem>$rnd[2]</elem>
              <elem>$rnd[3]</elem>
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
         <Mass>${mass}</Mass>
         <AnisoParam>
           <anisoP>${anisoP}</anisoP>
           <t_dir>3</t_dir>
           <xi_0>${xi_0}</xi_0>
           <nu>${nu}</nu>
         </AnisoParam>
         <FermionBC>
           <FermBC>SIMPLE_FERMBC</FermBC>
           <boundary>@{bc}</boundary>
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
EOF

} # foreach(t < Lt)

} # foreach n

# More stuff
print <<EOF;

    <elem>
      <annotation>
      ; STOCH_GROUP_MESON input file.
      </annotation>
      
      <Name>STOCH_GROUP_MESON</Name>
      <Frequency>1</Frequency>
      <Param> 
        <version>1</version>
        <creationOpContractType>SOURCE-SOURCE</creationOpContractType>
        <annihilationOpContractType>SOLUTION-SOLUTION</annihilationOpContractType>
        <mom2_max>0</mom2_max>
        <displacement_length>1</displacement_length>

        <QuarkSmearing>
          <!-- wvf_kind>NONE</wvf_kind -->
          <wvf_kind>GAUGE_INV_GAUSSIAN</wvf_kind>
          <wvf_param>2.0</wvf_param>
          <wvfIntPar>5</wvfIntPar>
          <no_smear_dir>3</no_smear_dir>
        </QuarkSmearing>

        <LinkSmearing>
          <!-- LinkSmearingType>NONE</LinkSmearingType -->
          <LinkSmearingType>STOUT_SMEAR</LinkSmearingType>
          <link_smear_fact>0.1625</link_smear_fact>
          <link_smear_num>4</link_smear_num>
          <no_smear_dir>3</no_smear_dir>
        </LinkSmearing>

        <QuarkDilutions>
	  <!-- First quark -->
          <elem>
            <version>1</version>
            <DilutionType>DILUTION_QUARK_SOURCE_CONST_FERM</DilutionType>
            <QuarkFiles>
              <TimeSliceFiles>
EOF

#  foreach $t (0 .. $Lt-1)
  foreach $t (0 .. 0)
  {
print<<EOF;
                <elem>
                  <DilutionFiles>
                    <elem>./zN_prop_q1_t${t}</elem>
                  </DilutionFiles>
                </elem>
EOF

}  # while(t < Lt)

print<<EOF;
              </TimeSliceFiles>
            </QuarkFiles>
          </elem>
	  <!-- Second quark -->
          <elem>
            <version>1</version>
            <DilutionType>DILUTION_QUARK_SOURCE_CONST_FERM</DilutionType>
            <QuarkFiles>
              <TimeSliceFiles>
EOF

#  foreach $t (0 .. $Lt-1)
  foreach $t (0 .. 0)
  {
print<<EOF;
                <elem>
                  <DilutionFiles>
                    <elem>./zN_prop_q2_t${t}</elem>
                  </DilutionFiles>
                </elem>
EOF

}  # while(t < Lt)

print <<EOF;
              </TimeSliceFiles>
            </QuarkFiles>
          </elem>
        </QuarkDilutions>
      </Param>
      <NamedObject>
        <gauge_id>default_gauge_field</gauge_id>
        <Quark_ids>ud</Quark_ids>
        <operators_file>
          <ops_file>two_displace</ops_file>
          <id>Double_Displaced</id>
        </operators_file>
      </NamedObject>
    </elem>

  </InlineMeasurements>
  <nrow>@{nrow}</nrow>
</Param>
<Cfg>
  <cfg_type>${gauge_type}</cfg_type>
  <cfg_file>${gauge_cfg}</cfg_file>
</Cfg>
</chroma>
EOF

exit(0);



