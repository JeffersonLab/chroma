#!/usr/bin/perl
#
# Used to generate the input file to run inline_stoch_group_baryon.cc which
# will create the source and sink ops for a specified 
# list of elemental operators. 
# 
#

#$gauge_type = "SZINQIO";
#$cfgNum = shift;
#$gauge_cfg = "wlq_6p1_12_48_xi3p0.lime${cfgNum}";
#$cfgNum > 0 or die "Usage: ./<progname> <cfgnum>";

$gauge_type = "WEAK_FIELD";
$gauge_cfg = "dummy";

@nrow = (2,2,2,2);

$anisoP = "false";
$mass = 0.0;
$xi_0 = 1;
$nu = 1;
@bc = (1, 1, 1, -1);
$RsdCG = 1.0e-10;
$MaxCG = 2000;

$N = 4;

$quark_id = "lll";

# Smearing Params
$sigma = 3.0;
$n_sigma = 32;

$rho = 0.15625;
$n_rho = 16;

$disp_len = 3;

$Elem_opFile = "Nucleon_elem";

#Dilution params 
$Nt = 1;
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

foreach $n (1, 2, 3)
{
  if ($n == 1) {@rnd = (471, 1694, 3965, 563);}
  if ($n == 2) {@rnd = (714, 1573, 3042, 517);}
  if ($n == 3) {@rnd = (123, 73, 42, 17);}

#  foreach $t (0 .. $Nt-1)
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
        <file_name>./zN_prop_q${n}_t${t}_sX_spc0.lime</file_name>	
        <file_volfmt>MULTIFILE</file_volfmt>
      </File>
    </elem>
EOF

} # foreach(t < Lt)

} # foreach n


# Make the operators

print <<EOF;

    <elem>
      <annotation>
      ; STOCH_GROUP_BARYON input file.
      </annotation>
     
      <Name>STOCH_GROUP_BARYON</Name>
      <Frequency>1</Frequency>
      <Param> 
        <version>2</version>
        <moms>
					<elem>0 0 0</elem>
			  </moms>
        <displacement_length>${disp_len}</displacement_length>

       	<QuarkSmearing>
 	  <wvf_kind>GAUGE_INV_GAUSSIAN</wvf_kind>
          <wvf_param>${sigma}</wvf_param>
          <wvfIntPar>${n_sigma}</wvfIntPar>
          <no_smear_dir>3</no_smear_dir>
        </QuarkSmearing>

	<LinkSmearing>
          <LinkSmearingType>STOUT_SMEAR</LinkSmearingType>
					 <link_smear_fact>${rho}</link_smear_fact>
					         <link_smear_num>${n_rho}</link_smear_num>
          <no_smear_dir>3</no_smear_dir>
	</LinkSmearing>
     	
	<QuarkDilutions>
EOF
foreach $n (1, 2, 3)
{
print <<EOF;
          <elem>
	    <DilutionType>DILUTION_QUARK_SOURCE_CONST_FERM</DilutionType>
	      <version>1</version>
	      <QuarkFiles>
	      <TimeSliceFiles>
EOF
foreach $t (0 .. $Nt-1)
{
print <<EOF;
                <elem>
		  <DilutionFiles>
EOF
#foreach $s (0 .. 0)
#{
	print <<EOF;
	            <elem>./zN_prop_q${n}_t${t}_sX_spc0.lime</elem>	
EOF
#} #s
print <<EOF;
                  </DilutionFiles>
                </elem>
EOF
} #t

print <<EOF;
              </TimeSliceFiles>
            </QuarkFiles>
          </elem>
EOF

}#foreach n

print <<EOF;
        </QuarkDilutions>	
      </Param>
      <NamedObject>
        <gauge_id>default_gauge_field</gauge_id>
        <operators_file>
          <ops_file>${Elem_opFile}</ops_file>
          <id>Nucleons</id>
        </operators_file>
        <Quark_ids>${quark_id}</Quark_ids>
      </NamedObject>
    </elem>
EOF
print <<EOF;
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



