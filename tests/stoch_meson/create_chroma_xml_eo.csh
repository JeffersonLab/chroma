#!/bin/tcsh

#set gauge_type = SZIN
set gauge_type = UNIT
set gauge_cfg = ../test_purgaug.cfg1

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

foreach c (0 1 2)
foreach s (0 1 2 3)

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

          <color_mask>${c}</color_mask>
          <spin_mask>${s}</spin_mask>
        </Source>
      </Param>
      <NamedObject>
        <source_id>zN_source</source_id>
      </NamedObject>
    </elem>

    <elem>
      <Name>PROPAGATOR_FERM</Name>
      <Frequency>1</Frequency>
      <Param>
        <version>8</version>
        <nonRelProp>false</nonRelProp>
        <obsvP>false</obsvP>
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
        <file_name>./zN_prop_q${n}_e_t${t}_c${c}_s${s}</file_name>
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

          <color_mask>${c}</color_mask>
          <spin_mask>${s}</spin_mask>
        </Source>
      </Param>
      <NamedObject>
        <source_id>zN_source</source_id>
      </NamedObject>
    </elem>

    <elem>
      <Name>PROPAGATOR_FERM</Name>
      <Frequency>1</Frequency>
      <Param>
        <version>8</version>
        <nonRelProp>false</nonRelProp>
        <obsvP>false</obsvP>
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
        <file_name>./zN_prop_q${n}_o_t${t}_c${c}_s${s}</file_name>
        <file_volfmt>MULTIFILE</file_volfmt>
      </File>
    </elem>

**EOF**

@ t ++

end  # while(t < Lt)

end  # foreach s
end  # foreach c
end  # foreach n

# Close up data file
cat << **EOF** 
  </InlineMeasurements>
  <nrow>${nrow}</nrow>
</Param>
<Cfg>
  <cfg_type>${gauge_type}</cfg_type>
  <cfg_file>${gauge_cfg}</cfg_file>
</Cfg>
</chroma>
**EOF**




