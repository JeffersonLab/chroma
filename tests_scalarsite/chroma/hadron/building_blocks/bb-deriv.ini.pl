#!/usr/bin/perl
#
# Used to generate the input file  bb-deriv.ini.xml
#


print <<"EOF";
<?xml version="1.0"?>
<chroma>
<annotation>
Sequential source
</annotation>
<Param> 
  <InlineMeasurements>

    <elem>
      <Name>MAKE_SOURCE</Name>
      <Frequency>1</Frequency>
      <Param>
        <version>6</version>
        <Source>
          <version>2</version>
          <SourceType>SHELL_SOURCE</SourceType>
          <j_decay>3</j_decay>
          <t_srce>0 0 0 6</t_srce>

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
        </Source>

      </Param>
      <NamedObject>
        <gauge_id>default_gauge_field</gauge_id>
        <source_id>sh_source_0</source_id>
      </NamedObject>
    </elem>

    <elem>
      <Name>PROPAGATOR</Name>
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
           <anisoP>false</anisoP>
           <t_dir>3</t_dir>
           <xi_0>1.0</xi_0>
           <nu>1.0</nu>
         </AnisoParam>
         <FermionBC>
           <FermBC>SIMPLE_FERMBC</FermBC>
           <boundary>1 1 1 0</boundary>
         </FermionBC>
        </FermionAction>
        <InvertParam>
          <invType>CG_INVERTER</invType>
          <RsdCG>1.0e-12</RsdCG>
          <MaxCG>1000</MaxCG>
        </InvertParam>
      </Param>
      <NamedObject>
        <gauge_id>default_gauge_field</gauge_id>
        <source_id>sh_source_0</source_id>
        <prop_id>sh_prop_0</prop_id>
      </NamedObject>
    </elem>

    <elem>
      <Name>SINK_SMEAR</Name>
      <Frequency>1</Frequency>
      <Param>
        <version>5</version>
        <Sink>
          <version>1</version>
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
      </Param>
      <NamedObject>
        <gauge_id>default_gauge_field</gauge_id>
        <prop_id>sh_prop_0</prop_id>
        <smeared_prop_id>sh_pt_prop_0</smeared_prop_id>
      </NamedObject>
    </elem>
EOF

# Loop over sequential sources
foreach $deriv_name ("a0-pionxNABLA_T1", "a0-a0xNABLA_T1", "a0-a0_2xNABLA_T1", "a0-pion_2xNABLA_T1", 
		     "a0-rhoxNABLA_A1", "a0-rhoxNABLA_T1", "a0-rhoxNABLA_T2", "a0-rhoxNABLA_E", 
		     "a0-rho_2xNABLA_A1", "a0-rho_2xNABLA_T1", "a0-rho_2xNABLA_T2", "a0-rho_2xNABLA_E", 
		     "a0-a1xNABLA_A1", "a0-a1xNABLA_T1", "a0-a1xNABLA_T2", "a0-a1xNABLA_E", 
		     "a0-b1xNABLA_A1", "a0-b1xNABLA_T1", "a0-b1xNABLA_T2", "a0-b1xNABLA_E", 
		     "a0-pionxD_T2", "a0-a0xD_T2", "a0-a0_2xD_T2", "a0-pion_2xD_T2",
		     "a0-a1xD_A2", "a0-a1xD_T1", "a0-a1xD_T2", "a0-a1xD_E", 
		     "a0-b1xD_A2", "a0-b1xD_T1", "a0-b1xD_T2", "a0-b1xD_E", 
		     "a0-rhoxD_A2", "a0-rhoxD_T1", "a0-rhoxD_T2", "a0-rhoxD_E", 
		     "a0-rho_2xD_A2", "a0-rho_2xD_T1", "a0-rho_2xD_T2", "a0-rho_2xD_E", 
		     "a0-pionxB_T1", "a0-a0xB_T1", "a0-a0_2xB_T1", "a0-pion_2xB_T1", 
		     "a0-rhoxB_A1", "a0-rhoxB_T1", "a0-rhoxB_T2", "a0-rhoxB_E", 
		     "a0-rho_2xB_A1", "a0-rho_2xB_T1", "a0-rho_2xB_T2", "a0-rho_2xB_E", 
		     "a0-a1xB_A1", "a0-a1xB_T1", "a0-a1xB_T2", "a0-a1xB_E",
		     "a0-b1xB_A1", "a0-b1xB_T1", "a0-b1xB_T2", "a0-b1xB_E")
{
#    print STDERR "deriv_name = $deriv_name\n";

    $rep = `echo $deriv_name |sed 's/-DERIV//'|sed 's/^.*_T1/T1/'|sed 's/^.*_T2/T2/'|sed 's/^.*_A1/A1/'|sed 's/^.*_A2/A2/'|sed 's/^.*_E/E/'`;

#    print STDERR "rep = $rep\n";

    if ($rep =~ /A1/)
    {
	@deriv_dirs = (-17);
    }
    elsif ($rep =~ /A2/)
    {
	@deriv_dirs = (-17);
    }
    elsif ($rep =~ /E/)
    {
	@deriv_dirs = (0, 1);
    }
    elsif ($rep =~ /T1/)
    {
	@deriv_dirs = (0, 1, 2);
    }
    elsif ($rep =~ /T2/)
    {
	@deriv_dirs = (0, 1, 2);
    }
    else
    {
	die "illegal rep = $rep\n";
    }

#    print STDERR "deriv_dirs = @deriv_dirs\n";

#    foreach $deriv_dir (@deriv_dirs) 
    $deriv_dir = 0;
    {
#	print STDERR "deriv_dir = $deriv_dir\n";


print <<"EOF";
    <elem>
      <annotation>
       ${deriv_name} seqsource
      </annotation>

      <Name>SEQSOURCE</Name>
      <Frequency>1</Frequency>
      <Param>
        <version>2</version>
        <SeqSource>
          <version>1</version>
          <SeqSourceType>${deriv_name}</SeqSourceType>
          <deriv_dir>${deriv_dir}</deriv_dir>
          <deriv_length>1</deriv_length>
          <j_decay>3</j_decay>
          <t_sink>6</t_sink>
          <sink_mom>1 0 0</sink_mom>
        </SeqSource>
      </Param>
      <PropSink>
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
      </PropSink>
      <NamedObject>
        <gauge_id>default_gauge_field</gauge_id>
        <prop_ids>
          <elem>sh_prop_0</elem>
        </prop_ids>
        <seqsource_id>seqsource_${deriv_name}</seqsource_id>
      </NamedObject>
    </elem>

    <elem>
      <annotation>
       ${deriv_name} seqprop
      </annotation>

      <Name>PROPAGATOR</Name>
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
           <anisoP>false</anisoP>
           <t_dir>3</t_dir>
           <xi_0>1.0</xi_0>
           <nu>1.0</nu>
         </AnisoParam>
         <FermionBC>
           <FermBC>SIMPLE_FERMBC</FermBC>
           <boundary>1 1 1 0</boundary>
         </FermionBC>
        </FermionAction>
        <InvertParam>
          <invType>CG_INVERTER</invType>
          <RsdCG>1.0e-12</RsdCG>
          <MaxCG>1000</MaxCG>
        </InvertParam>
      </Param>
      <NamedObject>
        <gauge_id>default_gauge_field</gauge_id>
        <source_id>seqsource_${deriv_name}</source_id>
        <prop_id>seqprop_${deriv_name}</prop_id>
      </NamedObject>
    </elem>

    <elem>
      <annotation>
        Erase the object
      </annotation>
      <Name>ERASE_NAMED_OBJECT</Name>
      <Frequency>1</Frequency>
      <NamedObject>
        <object_id>seqsource_${deriv_name}</object_id>
      </NamedObject>
    </elem>

    <elem>
      <annotation>
       BuildingBlock input file.
      </annotation>

      <Name>BUILDING_BLOCKS</Name>
      <Frequency>1</Frequency>
      <Param>
        <version>5</version>
        <use_sink_offset>true</use_sink_offset>
        <mom2_max>3</mom2_max>
        <links_max>1</links_max>
        <canonical>false</canonical>
        <time_reverse>false</time_reverse>
        <translate>false</translate>
        <FermState>
          <Name>SIMPLE_FERM_STATE</Name>
          <FermionBC>
            <FermBC>SIMPLE_FERMBC</FermBC>
            <boundary>1 1 1 -1</boundary>
          </FermionBC>
        </FermState>
      </Param>
      <BuildingBlocks>
       <OutFileName>./examplebb_v1.out</OutFileName>
       <GaugeId>default_gauge_field</GaugeId>
       <FrwdPropId>sh_prop_0</FrwdPropId>
       <BkwdProps>
         <elem>
           <BkwdPropId>seqprop_${deriv_name}</BkwdPropId>
           <BkwdPropG5Format>G5_B_G5</BkwdPropG5Format>
           <GammaInsertion>0</GammaInsertion>
           <Flavor>U</Flavor>
           <BBFileNamePattern>${deriv_name}_qz%c%1d_qy%c%1d_qx%c%1d.bb</BBFileNamePattern>
         </elem>
       </BkwdProps>
      </BuildingBlocks>
      <!-- Uncomment the xml_file line to write to a separate file -->
      <!-- xml_file>bb.dat.xml</xml_file -->
    </elem>

EOF

  }  # foreach
    
}  # foreach


# Finish up
print <<"EOF";
  </InlineMeasurements>
   <nrow>4 4 4 8</nrow>
</Param>

<RNG>
  <Seed>	
    <elem>11</elem>
    <elem>11</elem>
    <elem>11</elem>
    <elem>0</elem>
  </Seed>
</RNG>

<Cfg>
 <cfg_type>WEAK_FIELD</cfg_type>
 <cfg_file>dummy</cfg_file>
</Cfg>
</chroma>

EOF

exit(0);
