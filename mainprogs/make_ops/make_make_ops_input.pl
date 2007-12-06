#!/usr/bin/perl
#
# Used to generate the input file stoch_group_baryon_ops.ini.xml which
# will create the operator source and sink functions for a specified 
# list of operators. 
# 
#

@nrow = (8,8,8,16);


$Lt = $nrow[3];
$Ltm1 = $Lt - 1;

$NcoeffFiles = 1;
$NcoeffFilesm1 = $NcoeffFiles - 1;

$Nbins = 1;
$Nbinsm1 = $Nbins - 1;

$coeffFiles[0] = "/home/jbulava/all-to-all/tests/Nucleon_G1g_1";

$elemOpFilesRoot[0] = "/home/jbulava/all-to-all/tests/";

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
<Files>
	<Coeff_files>
EOF

#Files of Ops to make 
foreach $l (0 .. $NcoeffFilesm1)
{
	print<<EOF;
		<elem>${coeffFiles[$l]}</elem>
EOF
}

print<<EOF;
	</Coeff_files>
	<ElemOp_files>
EOF

#Roots of elemental Ops for each cfg
foreach $l (0 .. $Nbinsm1)
{
	print<<EOF;
		<elem>${elemOpFilesRoot[$l]}</elem>
EOF
}

print<<EOF;
	</ElemOp_files>
	<OpOutput_files>
EOF

#Roots of Ops to be ouput for each cfg
foreach $l (0 .. $Nbinsm1)
{
	print<<EOF;
		<elem>${outFilesRoot[$l]}</elem>
EOF
}

print<<EOF;
	</OpOutput_files>
</Files>
</MakeOps>
EOF
exit(0);
