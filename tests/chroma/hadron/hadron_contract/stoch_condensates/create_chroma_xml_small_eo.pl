#!/usr/bin/perl
#
# Used to generate the input file  stoch_condensates.ini.xml
#

$gauge_type = "WEAK_FIELD";
$gauge_cfg = "dummy";

$anisoP = "false";
$Kappa = 0.11;
$xi_0 = 1;
$nu = 1;
@bc = (1, 1, 1, -1);

$N = 4;

@nrow = (2, 2, 2, 4);

$RsdCG = 1.0e-10;
$MaxCG = 2000;

$Lt = $nrow[3];
$Ltm1 = $Lt - 1;

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

@rnd = (471, 1694, 3965, 563);

foreach $t (0 .. $Ltm1)
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
         <Kappa>0.11</Kappa>
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
        <file_name>./zN_prop_q1_t${t}</file_name>
        <file_volfmt>MULTIFILE</file_volfmt>
      </File>
    </elem>
EOF

} # foreach(t < Lt)

# More stuff
print <<EOF;

    <elem>
      <Name>HADRON_CONTRACT</Name>
      <Frequency>1</Frequency>
      <NamedObject>
        <gauge_id>default_gauge_field</gauge_id>
        <output_file>stoch_cond.lime</output_file>
        <Contractions>
          <annotation>
           ; STOCH_CONDENSATES
          </annotation>
          <elem>
            <version>1</version>
            <ContractionType>stoch_diagonal_gamma_condensates</ContractionType>
            <mom2_max>3</mom2_max>
            <avg_equiv_mom>true</avg_equiv_mom>
            <mom_origin>0 0 0</mom_origin>

            <LinkSmearing>
              <!-- LinkSmearingType>NONE</LinkSmearingType -->
              <LinkSmearingType>STOUT_SMEAR</LinkSmearingType>
              <link_smear_fact>0.1625</link_smear_fact>
              <link_smear_num>4</link_smear_num>
              <no_smear_dir>3</no_smear_dir>
            </LinkSmearing>
            <soln_files>
EOF

  foreach $t (0 .. $Ltm1)
{
print<<EOF;
              <elem>./zN_prop_q1_t${t}</elem>
EOF
}  # while(t < Lt)

# Finish up
print <<EOF;
            </soln_files>
          </elem>
        </Contractions>
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



