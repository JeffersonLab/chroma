#!/usr/bin/perl
# This version does not write diquark in disk. Instead for each diquark 
# all possible needed qqq are constructed keeping that diquark in memory. 
# This version is about 1.8 times faster than previous 
# diquark version (writing everything on disk). 
# This code needs the associated program "parse_baryons_qqqdiquark_fly.pl" from bin

#$[ = 0;		# set array base to 0
$, = ' ';		# set output field separator
$\ = "\n";		# set output record separator

die "Usage: run_qqqdiquark_fly.pl  <stem>  <spectro_stem>  <seqno>\n" unless scalar(@ARGV) == 3;

chdir("/scratch");

$stem = $ARGV[0]; chomp $stem;
$spectro_stem = $ARGV[2]; chomp $spectro_stem;
$seqno = $ARGV[1]; chomp $seqno;

$builddir = "/u/home/nilmani/chroma_test/exe_temp/";

require '/home/edwards/bin/queue_arch_type_maui.pl';
&determine_arch();

print "Arch = $arch";
print "Run = $run";
print "Host = ", `hostname`;

$qqq_props = "/home/nilmani/chroma_test/run/baryon/qqq_props";


$gauge_type = "SZINQIO";
$gauge_root = "wl_24_64_5p5_x2p38_um0p4125_cfg";
$gauge_cfg = "/scratch/${gauge_root}_${seqno}.lime";
$specdir = "/cache/LHPC/NF2/wl/aniso/5p5_24_64_xi3p0/qqq/cascade/um4125_sm3893";

`rcp /cache/LHPC/NF2/wl/aniso/5p5_24_64_xi3p0/cfgs/${gauge_root}_${seqno}.lime $gauge_cfg`;

$ud_mass = -0.4125;
$ud_nu   = 1.0;

$s_mass = -0.38922;
$s_nu   = 1.0;

$prop_root = "ud4125_sm3893";

$anisoP = "true";
#$rho  = 0.22;
#$n_rho = 2;
$xi   = 3.0;
$xi_0 = 2.38;
$u_s = 1.0;
$u_t = 1.0;
#$ud_clovCoeffR = $ud_nu / ($u_s * $u_s * $u_s);
#$ud_clovCoeffT =  0.5*($ud_nu + 1.0/$xi)/ ($u_t * $u_s * $u_s);

#$s_clovCoeffR = $s_nu / ($u_s * $u_s * $u_s);
#$s_clovCoeffT =  0.5*($s_nu + 1.0/$xi)/ ($u_t * $u_s * $u_s);

@nrow = (24, 24, 24, 64);
$t_source = 0;

#@nrow = (4, 4, 4, 8);
#$gauge_type = "WEAK_FIELD";

$bc = "1 1 1 -1";
$mom2_max = 5;

$src_smear_type = "GAUGE_INV_GAUSSIAN";
$wvf_param_i = 3;
$wvf_num_i = 32;

$snk_smear_type = "GAUGE_INV_GAUSSIAN";
$wvf_param_f = 3;
$wvf_num_f = 32;

$link_smear_type = "STOUT_SMEAR";
$link_smear_num = 16;
$link_smear_fact = 0.15625;


$RsdCG = 1.0e-10;
$MaxCG = 4000;

$volfmt = "MULTIFILE";

$input = "DATA";
$output = "${stem}.chroma.xml${seqno}";


#
# Build chroma input
#
open(INIT, "> ${input}");
print INIT<<"EOF";
<?xml version="1.0"?>

<chroma>
<Param> 
  <InlineMeasurements>
    <elem>
      <Name>LINK_SMEAR</Name>
      <Frequency>1</Frequency>
      <Param>
        <version>3</version>
        <LinkSmearingType>${link_smear_type}</LinkSmearingType>
        <link_smear_fact>${link_smear_fact}</link_smear_fact>
        <link_smear_num>${link_smear_num}</link_smear_num>
        <smear_dirs>1 1 1 0</smear_dirs>
      </Param>
      <NamedObject>
        <gauge_id>default_gauge_field</gauge_id>
        <linksmear_id>stout</linksmear_id>
      </NamedObject>
    </elem>
EOF
close(INIT);


#
#  Now construct the templates needed by the QQQ generator
#
open(SRC, "> src_template");
print SRC<<"EOF";
    <elem>
      <Name>MAKE_SOURCE</Name>
      <Frequency>1</Frequency>
      <Param>
        <version>6</version>
        <Source>
          <version>2</version>
          <SourceType>SHELL_SOURCE</SourceType>
          <j_decay>3</j_decay>
          <t_srce>0 0 0 ${t_source}</t_srce>

          <SmearingParam>
            <wvf_kind>${src_smear_type}</wvf_kind>
            <wvf_param>${wvf_param_i}</wvf_param>
            <wvfIntPar>${wvf_num_i}</wvfIntPar>
            <no_smear_dir>3</no_smear_dir>
          </SmearingParam>

          <Displacement>
           <version>1</version>
           <DisplacementType>SIMPLE_DISPLACEMENT</DisplacementType>
           <disp_length>_DISP_LENGTH</disp_length>
           <disp_dir>_DISP_DIR</disp_dir>
          </Displacement>
        </Source>
      </Param>
      <NamedObject>
        <gauge_id>stout</gauge_id>
        <source_id>_SOURCE_NAME</source_id>
      </NamedObject>
    </elem>
EOF
close(SRC);

open(PROP, "> ud_prop_template");
print PROP<<"EOF";
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
         <Mass>${ud_mass}</Mass>
         <AnisoParam>
           <anisoP>${anisoP}</anisoP>
           <t_dir>3</t_dir>
           <xi_0>${xi_0}</xi_0>
           <nu>${ud_nu}</nu>
         </AnisoParam>
         <FermState>
           <Name>SIMPLE_FERM_STATE</Name>
           <FermionBC>
             <FermBC>SIMPLE_FERMBC</FermBC>
             <boundary>${bc}</boundary>
           </FermionBC>
         </FermState>
        </FermionAction>
        <InvertParam>
          <invType>CG_INVERTER</invType>
          <RsdCG>${RsdCG}</RsdCG>
          <MaxCG>${MaxCG}</MaxCG>
        </InvertParam>
       </Param>
       <NamedObject>
        <gauge_id>default_gauge_field</gauge_id>
        <source_id>_SOURCE_NAME</source_id>
        <prop_id>_PROP_NAME</prop_id>
       </NamedObject>
    </elem>
EOF
close(PROP);

open(PROP, "> s_prop_template");
print PROP<<"EOF";
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
         <Mass>${s_mass}</Mass>
         <AnisoParam>
           <anisoP>${anisoP}</anisoP>
           <t_dir>3</t_dir>
           <xi_0>${xi_0}</xi_0>
           <nu>${s_nu}</nu>
         </AnisoParam>
         <FermState>
           <Name>SIMPLE_FERM_STATE</Name>
           <FermionBC>
             <FermBC>SIMPLE_FERMBC</FermBC>
             <boundary>${bc}</boundary>
           </FermionBC>
         </FermState>
        </FermionAction>
        <InvertParam>
          <invType>CG_INVERTER</invType>
          <RsdCG>${RsdCG}</RsdCG>
          <MaxCG>${MaxCG}</MaxCG>
        </InvertParam>
       </Param>
       <NamedObject>
        <gauge_id>default_gauge_field</gauge_id>
        <source_id>_SOURCE_NAME</source_id>
        <prop_id>_PROP_NAME</prop_id>
       </NamedObject>
    </elem>
EOF
close(PROP);

open(SNK, "> snk_template");
print SNK<<"EOF";
    <elem>
      <Name>SINK_SMEAR</Name>
      <Frequency>1</Frequency>
      <Param>
        <version>5</version>
        <Sink>
          <version>2</version>
          <SinkType>SHELL_SINK</SinkType>
          <j_decay>3</j_decay>

          <SmearingParam>
            <wvf_kind>${snk_smear_type}</wvf_kind>
            <wvf_param>${wvf_param_f}</wvf_param>
            <wvfIntPar>${wvf_num_f}</wvfIntPar>
            <no_smear_dir>3</no_smear_dir>
          </SmearingParam>

          <Displacement>
           <version>1</version>
           <DisplacementType>SIMPLE_DISPLACEMENT</DisplacementType>
           <disp_length>_DISP_LENGTH</disp_length>
           <disp_dir>_DISP_DIR</disp_dir>
          </Displacement>
        </Sink>
      </Param>
      <NamedObject>
        <gauge_id>stout</gauge_id>
        <prop_id>_PROP_NAME_IN</prop_id>
        <smeared_prop_id>_PROP_NAME_OUT</smeared_prop_id>
      </NamedObject>
    </elem>
EOF
close(SNK);

open(QQ, "> qq_template");
print QQ<<"EOF";
    <elem>
    <Name>DIQUARK</Name>
      <Param>
        <version>1</version>
        <Dirac_basis>true</Dirac_basis>
      </Param>
      <NamedObject>
        <gauge_id>default_gauge_field</gauge_id>
        <prop_ids>
          <elem>_PROPAGATOR1</elem>
          <elem>_PROPAGATOR2</elem>
        </prop_ids>
        <diquark_id>sh_sh__DIQUARKNAME</diquark_id>
      </NamedObject>
    </elem>

EOF
close(QQ);


open(QQQ, "> qqq_template");
print QQQ<<"EOF";

    <elem>
      <Name>QQQ_DIQUARK</Name>
      <Frequency>1</Frequency>
      <Param>
        <version>1</version>
        <Dirac_basis>true</Dirac_basis>
        <sparseP>false</sparseP>
        <nrow>@nrow</nrow>
      </Param>
      <NamedObject>
        <gauge_id>default_gauge_field</gauge_id>
        <prop_id>_PROPAGATOR3</prop_id>
        <diquark_id>sh_sh__DIQUARKNAME</diquark_id>
        <qqq_file>sh_sh__QQQNAME</qqq_file>
      </NamedObject>
    </elem>

EOF
close(QQQ);

# Displacement table needed for qqq
open(DISP, "> disp_table");
print DISP<<"EOF";
3
3
3
EOF
close(DISP);

#
#  Now create the qqq task list
#
print "Run the parse_baryons.pl script";
$data_tmp = "DATA.tmp";

`/u/home/nilmani/chroma_test/run/baryon/tmp4/parse_baryons_qqqdiquark_fly.pl ${data_tmp} ${qqq_props} disp_table src_template ud_prop_template s_prop_template snk_template ${prop_root} $seqno qq_template qqq_template`;

#  Append the data file
print "Append the ${data_tmp} script";
system("cat $data_tmp >> $input");


#
#  Append the closing tasks
# 
open(CFG, ">> $input");
print CFG<<"EOF";

    <elem>
      <Name>MAKE_SOURCE</Name>
      <Frequency>1</Frequency>
      <Param>
        <version>6</version>
        <Source>
          <version>2</version>
          <SourceType>SHELL_SOURCE</SourceType>
          <j_decay>3</j_decay>
          <t_srce>0 0 0 ${t_source}</t_srce>

          <SmearingParam>
            <wvf_kind>${src_smear_type}</wvf_kind>
            <wvf_param>${wvf_param_i}</wvf_param>
            <wvfIntPar>${wvf_num_i}</wvfIntPar>
            <no_smear_dir>3</no_smear_dir>
          </SmearingParam>

          <Displacement>
           <version>1</version>
           <DisplacementType>NONE</DisplacementType>
          </Displacement>
        </Source>
      </Param>
      <NamedObject>
        <gauge_id>stout</gauge_id>
        <source_id>sh_source</source_id>
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
         <Mass>${ud_mass}</Mass>
         <AnisoParam>
           <anisoP>${anisoP}</anisoP>
           <t_dir>3</t_dir>
           <xi_0>${xi_0}</xi_0>
           <nu>${ud_nu}</nu>
         </AnisoParam>
         <FermState>
           <Name>SIMPLE_FERM_STATE</Name>
           <FermionBC>
             <FermBC>SIMPLE_FERMBC</FermBC>
             <boundary>${bc}</boundary>
           </FermionBC>
         </FermState>
        </FermionAction>
        <InvertParam>
          <invType>CG_INVERTER</invType>
          <RsdCG>${RsdCG}</RsdCG>
          <MaxCG>${MaxCG}</MaxCG>
        </InvertParam>
       </Param>
       <NamedObject>
        <gauge_id>default_gauge_field</gauge_id>
        <source_id>sh_source</source_id>
        <prop_id>sh_ud_prop</prop_id>
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
         <Mass>${s_mass}</Mass>
         <AnisoParam>
           <anisoP>${anisoP}</anisoP>
           <t_dir>3</t_dir>
           <xi_0>${xi_0}</xi_0>
           <nu>${ud_nu}</nu>
         </AnisoParam>
         <FermState>
           <Name>SIMPLE_FERM_STATE</Name>
           <FermionBC>
             <FermBC>SIMPLE_FERMBC</FermBC>
             <boundary>${bc}</boundary>
           </FermionBC>
         </FermState>
        </FermionAction>
        <InvertParam>
          <invType>CG_INVERTER</invType>
          <RsdCG>${RsdCG}</RsdCG>
          <MaxCG>${MaxCG}</MaxCG>
        </InvertParam>
       </Param>
       <NamedObject>
        <gauge_id>default_gauge_field</gauge_id>
        <source_id>sh_source</source_id>
        <prop_id>sh_s_prop</prop_id>
       </NamedObject>
    </elem>

    <elem>
      <Name>ERASE_NAMED_OBJECT</Name>
      <NamedObject>
        <object_id>sh_source</object_id>
      </NamedObject>
    </elem>

    <elem>
      <Name>SINK_SMEAR</Name>
      <Frequency>1</Frequency>
      <Param>
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
      </Param>
      <NamedObject>
        <gauge_id>stout</gauge_id>
        <prop_id>sh_ud_prop</prop_id>
        <smeared_prop_id>sh_pt_ud_sink</smeared_prop_id>
      </NamedObject>
    </elem>

    <elem>
      <Name>SINK_SMEAR</Name>
      <Frequency>1</Frequency>
      <Param>
        <version>5</version>
        <Sink>
          <version>2</version>
          <SinkType>SHELL_SINK</SinkType>
          <j_decay>3</j_decay>

          <SmearingParam>
            <wvf_kind>${snk_smear_type}</wvf_kind>
            <wvf_param>${wvf_param_f}</wvf_param>
            <wvfIntPar>${wvf_num_f}</wvfIntPar>
            <no_smear_dir>3</no_smear_dir>
          </SmearingParam>

          <Displacement>
           <version>1</version>
           <DisplacementType>NONE</DisplacementType>
          </Displacement>
        </Sink>
      </Param>
      <NamedObject>
        <gauge_id>stout</gauge_id>
        <prop_id>sh_ud_prop</prop_id>
        <smeared_prop_id>sh_sh_ud_sink</smeared_prop_id>
      </NamedObject>
    </elem>

    <elem>
      <Name>SINK_SMEAR</Name>
      <Frequency>1</Frequency>
      <Param>
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
      </Param>
      <NamedObject>
        <gauge_id>stout</gauge_id>
        <prop_id>sh_s_prop</prop_id>
        <smeared_prop_id>sh_pt_s_sink</smeared_prop_id>
      </NamedObject>
    </elem>

    <elem>
      <Name>SINK_SMEAR</Name>
      <Frequency>1</Frequency>
      <Param>
        <version>5</version>
        <Sink>
          <version>2</version>
          <SinkType>SHELL_SINK</SinkType>
          <j_decay>3</j_decay>

          <SmearingParam>
            <wvf_kind>${snk_smear_type}</wvf_kind>
            <wvf_param>${wvf_param_f}</wvf_param>
            <wvfIntPar>${wvf_num_f}</wvfIntPar>
            <no_smear_dir>3</no_smear_dir>
          </SmearingParam>

          <Displacement>
           <version>1</version>
           <DisplacementType>NONE</DisplacementType>
          </Displacement>
        </Sink>
      </Param>
      <NamedObject>
        <gauge_id>stout</gauge_id>
        <prop_id>sh_s_prop</prop_id>
        <smeared_prop_id>sh_sh_s_sink</smeared_prop_id>
      </NamedObject>
    </elem>

    <elem>
      <Name>ERASE_NAMED_OBJECT</Name>
      <NamedObject>
        <object_id>sh_ud_prop</object_id>
      </NamedObject>
    </elem>

    <elem>
      <Name>ERASE_NAMED_OBJECT</Name>
      <NamedObject>
        <object_id>sh_s_prop</object_id>
      </NamedObject>
    </elem>

    <elem>
      <annotation>
        Current list of objects before hadspec
      </annotation>
      <Name>LIST_NAMED_OBJECT</Name>
    </elem>

    <elem>
      <Name>HADRON_SPECTRUM</Name>
      <Frequency>1</Frequency>
      <Param>
        <version>1</version>
        <MesonP>true</MesonP>
        <CurrentP>true</CurrentP>
        <BaryonP>true</BaryonP>
        <time_rev>true</time_rev>
        <mom2_max>${mom2_max}</mom2_max>
        <avg_equiv_mom>true</avg_equiv_mom>
      </Param>
      <NamedObject>
        <gauge_id>default_gauge_field</gauge_id>
        <sink_pairs>
          <elem>
            <first_id>sh_pt_ud_sink</first_id>
            <second_id>sh_pt_ud_sink</second_id>
          </elem>
          <elem>
            <first_id>sh_sh_ud_sink</first_id>
            <second_id>sh_sh_ud_sink</second_id>
          </elem>
          <elem>
            <first_id>sh_pt_ud_sink</first_id>
            <second_id>sh_pt_s_sink</second_id>
          </elem>
          <elem>
            <first_id>sh_sh_ud_sink</first_id>
            <second_id>sh_sh_s_sink</second_id>
          </elem>
        </sink_pairs>
      </NamedObject>
        <xml_file>${spectro_stem}.hadspec.xml${seqno}</xml_file>
    </elem>

  </InlineMeasurements>
  <nrow>@nrow</nrow>
</Param>
<Cfg>
  <cfg_type>${gauge_type}</cfg_type>
  <cfg_file>${gauge_cfg}</cfg_file>
</Cfg>
</chroma>
EOF
close(CFG);


#
# Run the chroma program
#
print "Before chroma: ", `date`;
system("$run $builddir/$arch/chroma -i $input -o $output < /dev/null");
print "After chroma: ", `date`;

die "Output file missing = $output\n" unless -f "$output";

$err = 0xffff & system("xmllint $output > /dev/null");
if ($err > 0x00)
{
    print "Some error with output file";
    exit(1);
}

#`mkdir -p ${specdir}`;

system("gzip -f9 $output");
system("mkdir ${seqno}");
system("tar -cf ${seqno}.tar sh_sh_qqq_*.lime${seqno}");
system("mv ${seqno}.tar ${seqno}");
system("gzip -f9 ${spectro_stem}.hadspec.xml${seqno}");
`rcp ${output}.gz ${specdir}`;
`rcp ${spectro_stem}.hadspec.xml${seqno}.gz ${specdir}`;
`rcp -r ${seqno} ${specdir}`;
#system("/bin/rm -f /scratch/*sh_sh*");
#chdir("/scratch");
#system("rm *.lime${seqno}*" );
#system("find . -name "*sh_sh*" -exec rm '{}' \;");

