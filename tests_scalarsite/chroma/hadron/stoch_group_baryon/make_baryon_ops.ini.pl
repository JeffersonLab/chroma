#!/usr/bin/perl
#
# Used to generate the input file stoch_group_baryon_ops.ini.xml which
# will create the operator source and sink functions for a specified 
# list of operators. 
# 
#

@nrow = (2,2,2,2);

#$elemOpFilesRoot[0] = "/home/jbulava/all-to-all/tests";
$elemOpFilesRoot[0] = ".";

$elemFiles[0] = "Nuc_323_001";
$elemFiles[1] = "Nuc_323_002";
$elemFiles[2] = "Nuc_323_003";
$elemFiles[3] = "Nuc_323_00m1";
$elemFiles[4] = "Nuc_323_00m2";
$elemFiles[5] = "Nuc_323_00m3";
$elemFiles[6] = "Nuc_332_000";
$elemFiles[7] = "Nuc_413_000";
$elemFiles[8] = "Nuc_413_001";
$elemFiles[9] = "Nuc_413_002";
$elemFiles[10] = "Nuc_413_003";
$elemFiles[11] = "Nuc_413_00m1";
$elemFiles[12] = "Nuc_413_00m2";
$elemFiles[13] = "Nuc_413_00m3";
$elemFiles[14] = "Nuc_431_000";

#$coeffFiles[0] = "/home/jbulava/all-to-all/tests/Nucleon_G1g_1";
$coeffFiles[0] = "./Nucleon_G1g_1";

#$outFilesRoot[0] = "/home/jbulava/all-to-all/tests/";
$outFilesRoot[0] = "./";


#
# Build chroma input
#

print <<EOF;
<?xml version="1.0"?>
<MakeOps>
<annotation>
;Input file for make_ops.cc 
; 
;
</annotation>
<Param> 
  <version>1</version>
  <Layout>@{nrow}</Layout> 
  <Decay_dir>3</Decay_dir>
</Param>
<InputFiles>
  <CoeffFiles>
EOF

#Files of Ops to make 
$NcoeffFilesm1 = $#coeffFiles;
foreach $l (0 .. $NcoeffFilesm1)
{
  print<<EOF;
    <elem>${coeffFiles[$l]}</elem>
EOF
}

print<<EOF;
  </CoeffFiles>
  <ElementalOpFiles>
EOF

$NelemM1 = $#elemFiles;
foreach $l (0 .. $NelemM1)
{
  print<<EOF;
    <elem>
      <Configs>
        <elem>
          <DilutionTimeSlices>
            <elem>
              <CreationOperatorFile>${elemOpFilesRoot[0]}/${elemFiles[$l]}_t0_src.lime</CreationOperatorFile>
              <AnnihilationOperatorFile>${elemOpFilesRoot[0]}/${elemFiles[$l]}_t0_snk.lime</AnnihilationOperatorFile>
            </elem>
          </DilutionTimeSlices>
        </elem>
      </Configs>
    </elem>
EOF
} #l

print<<EOF;
  </ElementalOpFiles>
</InputFiles>
<OutputInfo>
  <CfgOutputPaths>
    <elem>
      <SourceOpOutputPath>${outFilesRoot[0]}</SourceOpOutputPath>
      <SinkOpOutputPath>${outFilesRoot[0]}</SinkOpOutputPath>
    </elem>
  </CfgOutputPaths>
</OutputInfo>
</MakeOps>
EOF
exit(0);
