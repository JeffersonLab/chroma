#!/usr/bin/perl
#
# Used to generate the input file make_meson_ops.ini.xml which
# will create the operator source and sink functions for a specified 
# list of operators. 
# 
#

@nrow = (2, 2, 2, 4);
$Lt = $nrow[3];

#$elemOpFilesRoot[0] = "/home/jbulava/all-to-all/tests";
$elemOpFilesRoot[0] = ".";

push(@elemFiles, "m11_1_-1_t0");
push(@elemFiles, "m11_1_+1_t0");
push(@elemFiles, "m12_1_-3_t0");
push(@elemFiles, "m21_1_-3_t0");

#$coeffFiles[0] = "/home/jbulava/all-to-all/tests/Nucleon_G1g_1";
$coeffFiles[0] = "./test_meson";

$outFilesRoot[0] = "./";


#
# Build chroma input
#

print <<EOF;
<?xml version="1.0"?>
<MakeMesonOps>
<annotation>
;Input file make_meson_ops.cc 
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
$NcoeffFiles = $#coeffFiles;
foreach $l (0 .. $NcoeffFiles)
{
  print<<EOF;
    <elem>${coeffFiles[$l]}</elem>
EOF
}

print<<EOF;
  </CoeffFiles>
  <ElementalOpFiles>
EOF

$Nelem = $#elemFiles;
foreach $l (0 .. $Nelem)
{
  print<<EOF;
    <elem>
      <Configs>
        <elem>
          <DilutionTimeSlices>
            <elem>
              <CreationOperatorFile>${elemOpFilesRoot[0]}/${elemFiles[$l]}_src.lime</CreationOperatorFile>
              <AnnihilationOperatorFile>${elemOpFilesRoot[0]}/${elemFiles[$l]}_snk.lime</AnnihilationOperatorFile>
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
</MakeMesonOps>
EOF
exit(0);
