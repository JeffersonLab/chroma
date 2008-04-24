#!/usr/bin/perl
#
# Used to generate the input file to run inline_stoch_group_baryon.cc which
# will create the source and sink ops for a specified 
# list of elemental operators. 
# 
#

$gauge_type = "SZINQIO";
$cfgNum = shift;
$gauge_cfg = "wlq_6p1_12_48_xi3p0.lime${cfgNum}";

$cfgNum > 0 or die "Usage: ./<progname> <cfgnum>";

@nrow = (12,12,12,48);

$quark_id = "lll";

$Lt = $nrow[3];
$Ltm1 = $Lt - 1;

# Smearing Params
$sigma = 3.0;
$n_sigma = 32;

$rho = .15625;
$n_rho = 16;

$disp_len = 3;

$Elem_opFile = "Nucleon_elem";

#Dilution params 
$N_t = 1;
$Ntm1 = ${N_t} - 1;
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


# Make the operators

print <<EOF;

    <elem>
      <annotation>
      ; STOCH_GROUP_BARYON input file.
      </annotation>
     
      <Name>STOCH_GROUP_BARYON</Name>
     	<Frequency>1</Frequency>
      <Param> 
        <version>1</version>
        <mom2_max>0</mom2_max>
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
foreach $t (0 .. $Ntm1)
{
print <<EOF;
								<elem>
									<DilutionFiles>
EOF
foreach $s (0 .. 3)
{
	print <<EOF;
										<elem>./zN_prop_q${n}_t${t}_s${s}_spc0.lime</elem>	
										<elem>./zN_prop_q${n}_t${t}_s${s}_spc1.lime</elem>	
EOF
} #s
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



