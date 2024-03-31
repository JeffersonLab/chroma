[//]: # (Please format the XML snippets with xmllint --format)

# Chroma's superb tasks

The following tasks related with distillation are available in Chroma when compiling with superbblas' support. An easy way to compile chroma with superbblas support is with [chromaform] (https://github.com/eromero-vlc/chromaform):

	git clone https://github.com/eromero-vlc/chromaform
	cd chromaform
	chromaform --mg --superb chroma # install chroma with superbblas support for CPU
	chromaform --mg --superb-next chroma # install chroma with superbblas and mgproton support 
	chromaform --mg --superb --cuda chroma # install chroma with superbblas support for CUDA
	chromaform --mg --superb --hip chroma # install chroma with superbblas support for AMD GPUs
	
## Creation of the distillation basis

Example:
```
<?xml version="1.0"?>
<chroma>
  <Param>
    <InlineMeasurements>
      <elem>
        <Name>CREATE_COLORVECS_SUPERB</Name>
        <Frequency>1</Frequency>
        <Param>
          <num_vecs>16</num_vecs>
          <decay_dir>3</decay_dir>
          <t_start>0</t_start>
          <Nt_forward>4</Nt_forward>
          <phase>0 0 0</phase>
          <write_fingerprint>false</write_fingerprint>
          <LinkSmearing>
            <LinkSmearingType>STOUT_SMEAR</LinkSmearingType>
            <link_smear_fact>0.1</link_smear_fact>
            <link_smear_num>10</link_smear_num>
            <no_smear_dir>3</no_smear_dir>
          </LinkSmearing>
        </Param>
        <NamedObject>
          <gauge_id>default_gauge_field</gauge_id>
          <colorvec_out>eigs.sdb</colorvec_out>
        </NamedObject>
      </elem>
    </InlineMeasurements>
    <nrow>2 2 2 4</nrow>
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
    <cfg_type>DISORDERED</cfg_type>
    <cfg_file>dummy</cfg_file>
  </Cfg>
</chroma>
```

Main options:

* `Param/num_vecs`: distillation basis rank.
* `Param/LinkSmearing`: smearing applied to the links before computing the distillation basis.
* `NamedObject/colorvec_out`: file to store the distillation basis.
* `Param/write_fingerprint`: (optional, default `false`) whether to write down the field value for a few sites for each eigenvector instead of for all sites.
* `Param/t_start`: (optional, default `0`) first time slice to compute.
* `Param/Nt_forward`: (optional, default all time slices) number of time slices to compute after `t_start`.
* `Param/phase`: (optional, default `0 0 0`) spatial phasing applied to the computed eigenvectors.

### Notes

* When `write_fingerprint` is `false`, the whole vectors are stored in _filehash format_, and they can be used as input in any Chroma task. For large lattice size and number of eigenvectors, you may consider to set `write_fingerprint` to `true`. In that case, only a few sites for each eigenvector are going to be stored, which is enough to recover the eigenvectors on the fly (note that eigenvectors are unique up to a phase). Only superb tasks support distillation vectors in this format.

* Also all superb tasks support spatial phasing on the fly, so the option `phase` isn't usually needed.

* The smearing information is stored in the created file's metadata. When a filed created with `write_fingerprint` being `true` is passed as input to another Chroma tasks, the smearing information passed as input to the task should be the same as the smearing used for generating the eigenvectors. In case is not, the Chroma task will emit a runtime error.

# Creation of mesons

Example:
```
<?xml version="1.0"?>
<chroma>
  <Param>
    <InlineMeasurements>
      <elem>
        <Name>MESON_MATELEM_COLORVEC_SUPERB</Name>
        <Frequency>1</Frequency>
        <Param>
          <version>4</version>
          <use_derivP>true</use_derivP>
          <t_source>0</t_source>
          <Nt_forward>4</Nt_forward>
          <mom2_min>1</mom2_min>
          <mom2_max>2</mom2_max>
          <mom_list>
            <elem>0 0 0</elem>
            <elem>1 0 0</elem>
            <elem>-1 0 0</elem>
            <elem>0 1 0</elem>
            <elem>0 -1 0</elem>
            <elem>0 0 1</elem>
            <elem>0 0 -1</elem>
            <elem>2 0 0</elem>
            <elem>-2 0 0</elem>
            <elem>0 2 0</elem>
            <elem>0 -2 0</elem>
            <elem>0 0 2</elem>
            <elem>0 0 -2</elem>
            <elem>3 0 0</elem>
            <elem>-3 0 0</elem>
            <elem>0 3 0</elem>
            <elem>0 -3 0</elem>
            <elem>0 0 3</elem>
            <elem>0 0 -3</elem>
          </mom_list>
          <num_vecs>2</num_vecs>
          <phase>0 0 0</phase>
          <decay_dir>3</decay_dir>
          <displacement_list>
            <elem></elem>
            <elem>1</elem>
            <elem>2</elem>
            <elem>3</elem>
            <elem>1 1</elem>
            <elem>2 2</elem>
            <elem>3 3</elem>
            <elem>1 2</elem>
            <elem>1 3</elem>
            <elem>2 1</elem>
            <elem>2 3</elem>
            <elem>3 1</elem>
            <elem>3 2</elem>
          </displacement_list>
          <max_tslices_in_contraction>1</max_tslices_in_contraction>
          <LinkSmearing>
            <LinkSmearingType>STOUT_SMEAR</LinkSmearingType>
            <link_smear_fact>0.1</link_smear_fact>
            <link_smear_num>10</link_smear_num>
            <no_smear_dir>3</no_smear_dir>
          </LinkSmearing>
        </Param>
        <NamedObject>
          <gauge_id>default_gauge_field</gauge_id>
          <colorvec_files>
            <elem>eigs.sdb</elem>
          </colorvec_files>
          <meson_op_file>meson.sdb</meson_op_file>
        </NamedObject>
      </elem>
    </InlineMeasurements>
    <nrow>2 2 2 4</nrow>
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
    <cfg_type>DISORDERED</cfg_type>
    <cfg_file>dummy</cfg_file>
  </Cfg>
</chroma>
```

Main options:

* `Param/mom2_min` and `Param/mom2_max`: (ignored if `Param/mom_list` is given, default `0`) minimum and maximum modulus of the momenta to compute.
* `Param/mom_list`: (optional, alternative to `mom2_min` and `mom2_max`) list of the momenta to compute.
* `Param/use_derivP`: (optional, default `true`) whether the displacements are derivatives.
* `Param/displacement_list`: list of displacement/derivatives to compute.
* `Param/num_vecs`: rank of the distillation basis to use.
* `Param/decay_dir`: (should be 3) 
* `Param/t_source`: (optional, default `0`) first time slice to compute.
* `Param/Nt_forward`: (optional, default number of time slices) number of time slices to compute starting from `t_source`.
* `Param/max_tslices_in_contraction`: (optional, default all of them) maximum number of time slices to contract at once.
* `Param/max_moms_in_contraction`: (optional, default `1`) maximum number of momenta to contract at once.
* `Param/phase`: (optional, default `0 0 0`) spatial phasing on source and sink eigenvectors.
* `Param/quarkPhase` (optional, default `0 0 0`) spatial phasing on the source eigenvectors.
* `Param/aQuarkPhase` (optional, default `0 0 0`) spatial phasing on the sink eigenvectors.
* `Param/LinkSmearing`: smearing applied to the gauge links before doing the computation.
* `NamedObject/colorvec_files`: list of files with distillation basis.
* `NamedObject/meson_op_file`: file to store the results.

### Notes:

* The product of `Param/max_tslices_in_contraction` and `Param/max_moms_in_contraction` controls how much memory is used for contractions. Their tuning is more critical when using GPUs. The optimal value is usually the largest values of both parameters that the CPU or device can handle.

## Creation of baryon elementals

Example: 
```
<?xml version="1.0"?>
<chroma>
  <Param>
    <InlineMeasurements>
      <elem>
        <Name>BARYON_MATELEM_COLORVEC_SUPERB</Name>
        <Frequency>1</Frequency>
        <Param>
          <max_vecs>16</max_vecs>
          <t_source>0</t_source>
          <Nt_forward>4</Nt_forward>
          <mom_list>
            <elem>0 0 0</elem>
            <elem>1 0 0</elem>
            <elem>-1 0 0</elem>
            <elem>0 1 0</elem>
            <elem>0 -1 0</elem>
            <elem>0 0 1</elem>
            <elem>0 0 -1</elem>
            <elem>2 0 0</elem>
            <elem>-2 0 0</elem>
            <elem>0 2 0</elem>
            <elem>0 -2 0</elem>
            <elem>0 0 2</elem>
            <elem>0 0 -2</elem>
            <elem>3 0 0</elem>
            <elem>-3 0 0</elem>
            <elem>0 3 0</elem>
            <elem>0 -3 0</elem>
            <elem>0 0 3</elem>
            <elem>0 0 -3</elem>
          </mom_list>
          <decay_dir>3</decay_dir>
          <use_derivP>true</use_derivP>
          <displacement_list>
            <elem>
              <left>0</left>
              <middle>0</middle>
              <right>0</right>
            </elem>
            <elem>
              <left>0</left>
              <middle>0</middle>
              <right>1</right>
            </elem>
            <elem>
              <left>0</left>
              <middle>0</middle>
              <right>2</right>
            </elem>
            <elem>
              <left>0</left>
              <middle>0</middle>
              <right>3</right>
            </elem>
            <elem>
              <left>0</left>
              <middle>0</middle>
              <right>1 1</right>
            </elem>
            <elem>
              <left>0</left>
              <middle>0</middle>
              <right>2 2</right>
            </elem>
            <elem>
              <left>0</left>
              <middle>0</middle>
              <right>3 3</right>
            </elem>
            <elem>
              <left>0</left>
              <middle>0</middle>
              <right>1 2</right>
            </elem>
            <elem>
              <left>0</left>
              <middle>0</middle>
              <right>1 3</right>
            </elem>
            <elem>
              <left>0</left>
              <middle>0</middle>
              <right>2 1</right>
            </elem>
            <elem>
              <left>0</left>
              <middle>0</middle>
              <right>2 3</right>
            </elem>
            <elem>
              <left>0</left>
              <middle>0</middle>
              <right>3 1</right>
            </elem>
            <elem>
              <left>0</left>
              <middle>0</middle>
              <right>3 2</right>
            </elem>
            <elem>
              <left>0</left>
              <middle>1</middle>
              <right>1</right>
            </elem>
            <elem>
              <left>0</left>
              <middle>1</middle>
              <right>2</right>
            </elem>
            <elem>
              <left>0</left>
              <middle>1</middle>
              <right>3</right>
            </elem>
            <elem>
              <left>0</left>
              <middle>2</middle>
              <right>2</right>
            </elem>
            <elem>
              <left>0</left>
              <middle>2</middle>
              <right>3</right>
            </elem>
            <elem>
              <left>0</left>
              <middle>3</middle>
              <right>3</right>
            </elem>
          </displacement_list>
          <max_tslices_in_contraction>1</max_tslices_in_contraction>
          <max_moms_in_contraction>1</max_moms_in_contraction>
          <max_vecs>0</max_vecs>
          <use_superb_format>true</use_superb_format>
          <LinkSmearing>
            <LinkSmearingType>STOUT_SMEAR</LinkSmearingType>
            <link_smear_fact>0.1</link_smear_fact>
            <link_smear_num>10</link_smear_num>
            <no_smear_dir>3</no_smear_dir>
          </LinkSmearing>
        </Param>
        <NamedObject>
          <gauge_id>default_gauge_field</gauge_id>
          <colorvec_files>
            <elem>eigs.sdb</elem>
          </colorvec_files>
          <baryon_op_file>baryons.sdb</baryon_op_file>
        </NamedObject>
      </elem>
    </InlineMeasurements>
    <nrow>2 2 2 4</nrow>
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
    <cfg_type>DISORDERED</cfg_type>
    <cfg_file>dummy</cfg_file>
  </Cfg>
</chroma>
```

Main options:

* `Param/mom2_min` and `Param/mom2_max`: (ignored if `Param/mom_list` is given, default `0`) minimum and maximum modulus of the momenta to compute.
* `Param/mom_list`: (optional, alternative to `mom2_min` and `mom2_max`) list of the momenta to compute.
* `Param/use_derivP`: (optional, default `true`) whether the displacements are derivatives.
* `Param/displacement_list`: list of displacement/derivatives to compute.
* `Param/num_vecs`: rank of the distillation basis to use.
* `Param/decay_dir`: (should be 3) 
* `Param/t_source`: (ignored if `t_slices` is given, default `0`) first time slice to compute.
* `Param/Nt_forward`: (ignored if `t_slices` is given, default number of time slices) number of time slices to compute starting from `t_source`.
* `Param/t_slices`: (optional, default all time slices) list of time slices to compute.
* `Param/phase`: (optional, default `0 0 0`) spatial phasing on the eigenvectors.
* `Param/max_tslices_in_contraction`: (optional, default all of them) maximum number of time slices to contract at once.
* `Param/max_moms_in_contraction`: (optional, default all momenta) maximum number of momenta to contract at once.
* `Param/max_vecs`: (optional, default all eigenvectors) maximum number of eigenvectors to contract at once.
* `Param/use_superb_format`: (optional, default `false`) whether to use the superb file format for storing the elementals, instead of the filehash format.
* `Param/LinkSmearing`: smearing applied to the gauge links before doing the computation.
* `NamedObject/colorvec_files`: list of files with distillation basis.
* `NamedObject/baryon_op_file`: file to store the results.

### Notes:

* The product of `Param/max_tslices_in_contraction` and `Param/max_moms_in_contraction` and `Param/max_vecs` controls how much memory is used for contractions. Their tuning is more critical when using GPUs. The optimal value is usually the largest values of the parameters that the CPU or device can handle.

* In general storing the elementals in the superb file format is faster than in the traditional filehash format. The superb file format allows for multiple processes writing directly into the same file, require less IO operations, and is more efficient in storing the keys than filehash.

* Recent versions of redstar support baryons stored in the new superb format.

## Creation of perambulators

Example:
```
<?xml version="1.0"?>
<chroma>
  <Param>
    <InlineMeasurements>
      <elem>
        <Name>PROP_AND_MATELEM_DISTILLATION_SUPERB</Name>
        <Frequency>1</Frequency>
        <Param>
          <Contractions>
            <mass_label>U-0.2390</mass_label>
            <num_vecs>16</num_vecs>
            <t_sources>1</t_sources>
            <Nt_forward>1</Nt_forward>
            <Nt_backward>1</Nt_backward>
            <decay_dir>3</decay_dir>
            <use_device_for_contractions>true</use_device_for_contractions>
            <max_rhs>1</max_rhs>
          </Contractions>
          <Propagator>
            <version>10</version>
            <quarkSpinType>FULL</quarkSpinType>
            <obsvP>false</obsvP>
            <numRetries>1</numRetries>
            <FermionAction>
              <FermAct>CLOVER</FermAct>
              <Mass>0.1</Mass>
              <clovCoeff>1.20536588031793</clovCoeff>
              <FermState>
                <Name>STOUT_FERM_STATE</Name>
                <rho>0.125</rho>
                <n_smear>1</n_smear>
                <orthog_dir>-1</orthog_dir>
                <FermionBC>
                  <FermBC>SIMPLE_FERMBC</FermBC>
                  <boundary>1 1 1 -1</boundary>
                </FermionBC>
              </FermState>
            </FermionAction>
            <InvertParam>
              <invType>BICGSTAB_INVERTER</invType>
              <RsdBiCGStab>1.0e-7</RsdBiCGStab>
              <MaxBiCGStab>20000</MaxBiCGStab>
              <Verbose>true</Verbose>
            </InvertParam>
          </Propagator>
        </Param>
        <NamedObject>
          <gauge_id>default_gauge_field</gauge_id>
          <colorvec_files>
            <elem>eigs.sdb</elem>
          </colorvec_files>
          <prop_op_file>props.sdb</prop_op_file>
        </NamedObject>
      </elem>
    </InlineMeasurements>
    <nrow>2 2 2 4</nrow>
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
    <cfg_type>DISORDERED</cfg_type>
    <cfg_file>dummy</cfg_file>
  </Cfg>
</chroma>
```

Main options:

* `Param/Contractions/num_vecs`: rank of the distillation basis to use.
* `Param/Contractions/t_sources`: list of time sources to compute.
* `Param/Contractions/decay_dir`: (should be 3) 
* `Param/Contractions/Nt_forward`: number of sink time slices to compute starting from each value in `t_sources`.
* `Param/Contractions/Nt_backward`: number of sink time slices to compute previous to each value in `t_sources`.
* `Param/Contractions/mass_label`: mass label, used by `redstar`.
* `Param/Contractions/max_rhs`: (optional, default is 8) number of right-hand-sides to invert at once.
* `Param/Contractions/phase`: (optional, default `0 0 0`) spatial phasing on the eigenvectors.
* `Param/Contractions/use_device_for_contractions`: (optional, default `true`) whether to use the GPUs for the contractions if available.
* `Param/Propagator`: propagator configuration.
* `NamedObject/colorvec_files`: list of files with distillation basis.
* `NamedObject/prop_op_file`: file to store the results.

### Notes:

* The option `Param/Contractions/max_rhs` being larger than 1 may accelerate the invertor when using `MGPROTO` or `mgproton` invertors; `quda` is in progress.

* As a last resource, set `Param/Contractions/use_device_for_contractions` to `false` to save GPU memory, at the cost of slower contractions.

## Generate the generalize perambulator

Example:
```
<?xml version="1.0"?>
<chroma>
  <Param>
    <InlineMeasurements>
      <elem>
        <Name>UNSMEARED_HADRON_NODE_DISTILLATION_SUPERB</Name>
        <Frequency>1</Frequency>
        <Param>
          <DispGammaMomList>
            <elem>
              <gamma>0</gamma>
              <mom>0 0 0</mom>
              <displacement/>
            </elem>
            <elem>
              <gamma>0</gamma>
              <mom>0 0 0</mom>
              <displacement>1</displacement>
            </elem>
          </DispGammaMomList>
          <SinkSourcePairs>
            <elem>
              <t_source>0</t_source>
              <t_sink>0</t_sink>
              <Nt_forward>1</Nt_forward>
              <Nt_backward>0</Nt_backward>
            </elem>
            <elem>
              <t_source>0</t_source>
              <t_sink>1</t_sink>
              <Nt_forward>1</Nt_forward>
              <Nt_backward>0</Nt_backward>
            </elem>
          </SinkSourcePairs>
          <Contractions>
            <mass_label>U-0.2390</mass_label>
            <num_vecs>16</num_vecs>
            <Nt_forward>1</Nt_forward>
            <decay_dir>3</decay_dir>
            <max_rhs>8</max_rhs>
            <use_device_for_contractions>true</use_device_for_contractions>
            <use_genprop5_format>true</use_genprop5_format>
            <use_derivP>false</use_derivP>
          </Contractions>
          <Propagator>
            <version>10</version>
            <quarkSpinType>FULL</quarkSpinType>
            <obsvP>false</obsvP>
            <numRetries>1</numRetries>
            <FermionAction>
              <FermAct>CLOVER</FermAct>
              <Mass>0.1</Mass>
              <clovCoeff>1.20536588031793</clovCoeff>
              <FermState>
                <Name>STOUT_FERM_STATE</Name>
                <rho>0.125</rho>
                <n_smear>1</n_smear>
                <orthog_dir>-1</orthog_dir>
                <FermionBC>
                  <FermBC>SIMPLE_FERMBC</FermBC>
                  <boundary>1 1 1 -1</boundary>
                </FermionBC>
              </FermState>
            </FermionAction>
            <InvertParam>
              <invType>BICGSTAB_INVERTER</invType>
              <RsdBiCGStab>1.0e-7</RsdBiCGStab>
              <MaxBiCGStab>20000</MaxBiCGStab>
              <Verbose>true</Verbose>
            </InvertParam>
          </Propagator>
        </Param>
        <NamedObject>
          <gauge_id>default_gauge_field</gauge_id>
          <colorvec_files>
            <elem>eigs.sdb</elem>
          </colorvec_files>
          <dist_op_file>gprops.sdb</dist_op_file>
        </NamedObject>
      </elem>
    </InlineMeasurements>
    <nrow>2 2 2 4</nrow>
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
    <cfg_type>DISORDERED</cfg_type>
    <cfg_file>dummy</cfg_file>
  </Cfg>
</chroma>
```

Computing generalized perambulators is challenging and the XML interface reflects years of accumulated hacks.

Main options:

* `Param/Contractions/decay_dir`: (should be 3) 
* `Param/Contractions/use_derivP`: whether the displacements are derivatives (usually `false`).
* `Param/Contractions/num_vecs`: (ignored if `` is given) rank of the distillation basis to use.
* `Param/Contractions/mass_label`: mass label, used by `redstar`.
* `Param/Contractions/max_rhs`: (optional, default is 8) number of right-hand-sides to invert at once.
* `Param/Contractions/max_tslices_in_contraction`: (optional, default all of them) maximum number of time slices to contract at once.
* `Param/Contractions/max_moms_in_contraction`: (optional, default `1`) maximum number of momenta to contract at once.
* `Param/Contractions/use_device_for_contractions`: (optional, default `true`) whether to use the GPUs for the contractions if available.
* `Param/Contractions/use_genprop4_format`: (optional, default `false`) whether to use genprop4 filehash format to store the results.
* `Param/Contractions/use_genprop5_format`: (optional, default `false`) whether to use superb file format to store the results.
* `Param/Contractions/phase`: (optional, default `0 0 0`) spatial phasing on source and sink eigenvectors.
* `Param/Contractions/quarkPhase` (optional, default `0 0 0`) spatial phasing on the source eigenvectors.
* `Param/Contractions/aQuarkPhase` (optional, default `0 0 0`) spatial phasing on the sink eigenvectors.
* `Param/Propagator`: propagator configuration.
* `NamedObject/colorvec_files`: list of files with distillation basis.
* `NamedObject/dist_op_file`: file to store the results.

Momenta, gammas, displacements:

* `Param/DispGammaMomList`: (optional, default compute all gammas) list of gamma, momentum and displacement to compute.
* `Param/Displacements`: (optional, complement to `Param/DispGammaMomList`) list of displacements to compute.
* `Param/Moms`: (optional, complement to `Param/DispGammaMomList`) list of momenta to compute.

Source, middle, and sink time slices:

* `Param/SinkSourcePairs`: (optional) list of source-sink time slices and the range of middle time slices to compute.
* `Param/SinkSourcePairs/elem/t_source`: time slice source to compute.
* `Param/SinkSourcePairs/elem/t_sink`: time slice sink to compute.
* `Param/SinkSourcePairs/elem/Nt_backward`: middle time slices to compute since source backward.
* `Param/SinkSourcePairs/elem/Nt_forward`: middle time slices to compute since source forward.

Alternative way to express the source, middle, and sink time slices to compute:

* `Param/Contractions/t_start`: (optional, default `0`) origin time slice for `Param/SinkSources`.
* `Param/Contractions/Nt_forward`: (optional, default all time slices) number of time slices to compute starting at each source.
* `Param/SinkSources`: (optional, complement to `Param/SinkSourcePairs`) list of sink-source pairs as an adjacency map.
* `Param/SinkSources/elem/key`: time source to compute.
* `Param/SinkSources/elem/value`: list of time sinks to compute for the time source. The number of middle slices is given by `Param/Contractions/Nt_forward`.

### Notes:

* The option `Param/Contractions/max_rhs` being larger than 1 may accelerate the invertor when using `MGPROTO` or `mgproton` invertors; `quda` is in progress.

* The product of `Param/max_tslices_in_contraction` and `Param/max_moms_in_contraction` controls how much memory is used for contractions. Their tuning is more critical when using GPUs. The optimal value is usually the largest values of both parameters that the CPU or device can handle.

* As a last resource, set `Param/Contractions/use_device_for_contractions` to `false` to save GPU memory, at the cost of slower contractions.

## Generate disconnected diagrams

Example:
```
<?xml version="1.0"?>
<chroma>
  <Param>
    <InlineMeasurements>
      <elem>
        <Name>DISCO_PROBING_DEFLATION_SUPERB</Name>
        <Param>
          <max_path_length>1</max_path_length>
          <mom2_max>2</mom2_max>
          <mass_label>U-0.2070</mass_label>
          <probing_distance>0</probing_distance>
          <probing_power>2</probing_power>
          <noise_vectors>10</noise_vectors>
          <max_rhs>1</max_rhs>
          <Propagator>
            <version>10</version>
            <quarkSpinType>FULL</quarkSpinType>
            <obsvP>false</obsvP>
            <numRetries>1</numRetries>
            <FermionAction>
              <FermAct>CLOVER</FermAct>
              <Mass>0.1</Mass>
              <clovCoeff>1.20536588031793</clovCoeff>
              <FermState>
                <Name>STOUT_FERM_STATE</Name>
                <rho>0.125</rho>
                <n_smear>1</n_smear>
                <orthog_dir>-1</orthog_dir>
                <FermionBC>
                  <FermBC>SIMPLE_FERMBC</FermBC>
                  <boundary>1 1 1 -1</boundary>
                </FermionBC>
              </FermState>
            </FermionAction>
            <InvertParam>
              <invType>BICGSTAB_INVERTER</invType>
              <RsdBiCGStab>1.0e-7</RsdBiCGStab>
              <MaxBiCGStab>20000</MaxBiCGStab>
              <Verbose>true</Verbose>
            </InvertParam>
          </Propagator>
          <Projector>
            <projectorType>MG_PROTO_QPHIX_CLOVER_PROJECTOR</projectorType>
            <CloverParams>
              <Mass>0.1</Mass>
              <clovCoeff>1.20536588031793</clovCoeff>
              <AnisoParam>
                <anisoP>false</anisoP>
                <t_dir>3</t_dir>
                <xi_0>1</xi_0>
                <nu>1</nu>
              </AnisoParam>
            </CloverParams>
            <AntiPeriodicT>true</AntiPeriodicT>
            <MGLevels>3</MGLevels>
            <Blocking>
              <elem>4 4 4 4</elem>
              <elem>2 2 2 2</elem>
            </Blocking>
            <NullVecs>24 32</NullVecs>
            <NullSolverMaxIters>800 800</NullSolverMaxIters>
            <NullSolverRsdTarget>-0.002 -0.0002</NullSolverRsdTarget>
            <NullSolverVerboseP>0 0</NullSolverVerboseP>
            <EigenSolverBlockSize>1</EigenSolverBlockSize>
            <EigenSolverMaxRestartSize>32</EigenSolverMaxRestartSize>
            <EigenSolverMaxRank>1600</EigenSolverMaxRank>
            <EigenSolverRsdTarget>1.0e-3</EigenSolverRsdTarget>
            <EigenSolverMaxIters>0</EigenSolverMaxIters>
            <EigenSolverVerboseP>true</EigenSolverVerboseP>
            <BottomSolverNKrylov>40</BottomSolverNKrylov>
            <BottomSolverRsdTarget>1.0e-4</BottomSolverRsdTarget>
            <BottomSolverMaxIters>10000</BottomSolverMaxIters>
            <BottomSolverVerboseP>false</BottomSolverVerboseP>
            <SubspaceId>foo_eo_proj</SubspaceId>
          </Projector>
          <use_ferm_state_link>true</use_ferm_state_link>
        </Param>
        <NamedObject>
          <gauge_id>default_gauge_field</gauge_id>
          <sdb_file>disco.sdb</sdb_file>
        </NamedObject>
      </elem>
    </InlineMeasurements>
    <nrow>2 2 2 4</nrow>
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
    <cfg_type>DISORDERED</cfg_type>
    <cfg_file>dummy</cfg_file>
  </Cfg>
</chroma>
```

Main options:

* `Param/max_path_length`: maximum displacement in the z direction to compute.
* `Param/mom2_min` and `Param/mom2_max`: (ignored if `Param/mom_list` is given, default `0`) minimum and maximum modulus of the momenta to compute.
* `Param/mom_list`: (optional, alternative to `mom2_min` and `mom2_max`) list of the momenta to compute.
* `Param/mass_label`: mass label, used by `redstar`.
* `Param/probing_distance`: shifting used in probing.
* `Param/probing_power`: use a probing scheme that zeroes all neighbors up to `probing_power` distance after being shifted in the z direction `probing_distance` units.
* `Param/noise_vectors`: (optional, default `1`) number of noise vectors in estimating the traces.
* `Param/max_rhs`: (optional, default is `1`) number of right-hand-sides to invert at once.
* `Param/use_ferm_state_link`: (should be `true`).
* `Param/Propagator`: propagator configuration.
* `Param/Projector`: projector configuration.

## Singular value/eigenvalue partial decomposition

Compute the right singular vectors and values of the lower part of the spectrum of an operator by computing the eigenpairs of the Hermitian operator $D^{-1}\gamma_5$ for the largest eigenvalues in magnitude.

Example:
```
<?xml version="1.0"?>
<chroma>
  <Param>
    <InlineMeasurements>
      <elem>
        <Name>EIGENVALUES_SUPERB</Name>
        <Param>
          <Contractions>
            <mass_label>U-0.2070</mass_label>
            <num_vecs>5</num_vecs>
            <tolerance>0.1</tolerance>
            <max_rhs>1</max_rhs>
          </Contractions>
          <eigensolver>
            <max_block_size>1</max_block_size>
            <verbosity>detailed</verbosity>
          </eigensolver>
          <Propagator>
            <version>10</version>
            <quarkSpinType>FULL</quarkSpinType>
            <obsvP>false</obsvP>
            <numRetries>1</numRetries>
            <FermionAction>
              <FermAct>CLOVER</FermAct>
              <Mass>0.1</Mass>
              <clovCoeff>1.20536588031793</clovCoeff>
              <FermState>
                <Name>STOUT_FERM_STATE</Name>
                <rho>0.125</rho>
                <n_smear>1</n_smear>
                <orthog_dir>-1</orthog_dir>
                <FermionBC>
                  <FermBC>SIMPLE_FERMBC</FermBC>
                  <boundary>1 1 1 -1</boundary>
                </FermionBC>
              </FermState>
            </FermionAction>
            <InvertParam>
              <invType>BICGSTAB_INVERTER</invType>
              <RsdBiCGStab>1.0e-7</RsdBiCGStab>
              <MaxBiCGStab>20000</MaxBiCGStab>
              <Verbose>true</Verbose>
            </InvertParam>
          </Propagator>
        </Param>
        <NamedObject>
          <gauge_id>default_gauge_field</gauge_id>
          <eigs_file>eigenpairs.sdb</eigs_file>
        </NamedObject>
      </elem>
    </InlineMeasurements>
    <nrow>2 2 2 4</nrow>
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
    <cfg_type>DISORDERED</cfg_type>
    <cfg_file>dummy</cfg_file>
  </Cfg>
</chroma>
```

Main options:

* `Param/Contractions/num_vecs`: number of eigenpairs to compute.
* `Param/Contractions/tolerance`: $\|\gamma_5 D^{-1} x - \lambda x \|_2 \leq \text{tol}\ \|D^{-1}\|_2$.
* `Param/Contractions/mass_label`: label indicating the used mass.
* `Param/Contractions/max_rhs`: (optional, default is 8) number of right-hand-sides to invert at once.
* `Param/eigensolver/max_block_size`: (optional, default is `1`) maximum number of vectors expanding the search subspace in each iteration.
* `Param/eigensolver/max_basis_size`: (optional, default is `PRIMME`'s default) maximum rank of the search subspace.
* `Param/eigensolver/verbosity`: (optional, default `nooutput`) one of `nooutput`, `summary`, `detailed`.
* `Param/Propagator`: propagator configuration.
* `NamedObject/eigs_file`: file to store the results.

# `mgproton` solver collection

## Flexible GMRES

Example:
```
<InvertParam>
  <invType>MGPROTON</invType>
  <type>fgmres</type>
  <tol>1e-7</tol>
  <max_basis_size>3</max_basis_size>
  <max_its>20000</max_its>
  <verbosity>Detailed</verbosity>
</InvertParam>
```

Main options:

* `tol`: residual norm tolerance, stopping when $\|Dx-b\|_2 \leq \text{tol}\ \|b\|_2$.
* `max_basis_size`: (optional, default `5`) maximum size of the search subspace.
* `max_its`: (optional, default infinity) maximum number of iterations.
* `error_if_not_converged`: (optional, default `true`) whether to complain if tolerance is not achieved.
* `prec`: (optional, default none) left preconditioning, does not affect residual norm.
* `ortho_each_its`: (optional, default `8` for double and `4` for single precision) orthogonalize the basis every this number of iterations.
* `max_residual_updates`: (optional, default `4` for double and `2` for single precision) recompute the residual vector every this number of iterations.
* `max_simultaneous_rhs`: (optional, default infinity) solver this many right-hand-sides at once.
* `verbosity`: (optional, default `nooutput`) level of verbosity, one of `nonoutput`, `summary`, `detailed`.
* `prefix`: (optional, default none) prefix output related with this solver with this string.

## BiCGstab

Example:
```
<InvertParam>
  <invType>MGPROTON</invType>
  <type>bicgstab</type>
  <tol>1e-7</tol>
  <max_its>20000</max_its>
  <verbosity>Detailed</verbosity>
</InvertParam>
```

Main options:

* `tol`: residual norm tolerance, stopping when $\|Dx-b\|_2 \leq \text{tol}\ \|b\|_2$.
* `max_its`: (optional, default infinity) maximum number of iterations.
* `error_if_not_converged`: (optional, default `true`) whether to complain if tolerance is not achieved.
* `prec`: (optional, default none) right preconditioning, may affect residual norm.
* `max_simultaneous_rhs`: (optional, default infinity) solver this many right-hand-sides at once.
* `verbosity`: (optional, default `nooutput`) level of verbosity, one of `nonoutput`, `summary`, `detailed`.
* `prefix`: (optional, default none) prefix output related with this solver with this string.

## Minimum Residual (MR)

Example:
```
<InvertParam>
  <invType>MGPROTON</invType>
  <type>mr</type>
  <tol>1e-7</tol>
  <max_its>20000</max_its>
  <verbosity>Detailed</verbosity>
</InvertParam>
```

Main options:

* `tol`: residual norm tolerance, stopping when $\|Dx-b\|_2 \leq \text{tol}\ \|b\|_2$.
* `max_its`: (optional, default infinity) maximum number of iterations.
* `error_if_not_converged`: (optional, default `true`) whether to complain if tolerance is not achieved.
* `prec`: (optional, default none) left preconditioning, does not affect residual norm.
* `max_simultaneous_rhs`: (optional, default infinity) solver this many right-hand-sides at once.
* `verbosity`: (optional, default `nooutput`) level of verbosity, one of `nonoutput`, `summary`, `detailed`.
* `prefix`: (optional, default none) prefix output related with this solver with this string.

## Generalized Conjugate Residual (GCR)

Example:
```
<InvertParam>
  <invType>MGPROTON</invType>
  <type>gcr</type>
  <max_basis_size>3</max_basis_size>
  <tol>1e-7</tol>
  <max_its>20000</max_its>
  <verbosity>Detailed</verbosity>
</InvertParam>
```

Main options:

* `tol`: residual norm tolerance, stopping when $\|Dx-b\|_2 \leq \text{tol}\ \|b\|_2$.
* `max_basis_size`: (optional, default `3`) maximum size of the search subspace.
* `max_its`: (optional, default infinity) maximum number of iterations.
* `error_if_not_converged`: (optional, default `true`) whether to complain if tolerance is not achieved.
* `prec`: (optional, default none) left preconditioning, does not affect residual norm.
* `max_simultaneous_rhs`: (optional, default infinity) solver this many right-hand-sides at once.
* `verbosity`: (optional, default `nooutput`) level of verbosity, one of `nonoutput`, `summary`, `detailed`.
* `prefix`: (optional, default none) prefix output related with this solver with this string.

## Even-odd preconditioning

Approximate $D^{-1}$ by splitting the sites of $D$ into two colors (red-black, even-odd) and solving the Schur complement iteratively.

Example:
```
<InvertParam>
  <invType>MGPROTON</invType>
  <type>eo</type>
  <solver>
    <type>gcr</type>
    <tol>1e-7</tol>
    <max_basis_size>3</max_basis_size>
    <max_its>20000</max_its>
    <verbosity>Detailed</verbosity>
  </solver>
</InvertParam>
```

Main options:

* `use_Aee_prec`: (optional, default `false`) whether to preconditioning the Schur complement with $D_{ee}^{-1}$, resulting in solving $(I-D_{eo}D_{oo}^{-1}D_{oe}D_{ee}^{-1})^{-1}$, insted of $(D_{ee}-D_{eo}D_{oo}^{-1}D_{oe})^{-1}$.
* `prec_ee`: (optional, default none): left preconditioning acting on the even sites for solving the Schur complement.
* `solver`: solver for the Schur complement.
* `prefix`: (optional, default none) prefix output related with this solver with this string.

## Multigrid preconditioner

Example:
```
<InvertParam>
  <invType>MGPROTON</invType>
  <type>eo</type>
  <solver>
    <type>mr</type>
    <tol>1e-7</tol>
    <max_its>20000</max_its>
    <prefix>l0</prefix>
    <verbosity>Detailed</verbosity>
  </solver>
  <prec_ee>
    <type>mg</type>
    <num_null_vecs>24</num_null_vecs>
    <blocking>4 4 4 4</blocking>
    <solver_null_vecs>
      <type>eo</type>
      <solver>
        <type>mr</type>
        <tol>1e-3</tol>
        <max_its>50</max_its>
        <error_if_not_converged>false</error_if_not_converged>
        <prefix>nv0</prefix>
        <verbosity>Detailed</verbosity>
        <prec>
          <type>dd</type>
          <solver>
            <type>mr</type>
            <max_its>2</max_its>
            <tol>1e-1</tol>
            <error_if_not_converged>false</error_if_not_converged>
            <verbosity>Detailed</verbosity>
            <prefix>nv0_dd</prefix>
          </solver>
        </prec>
      </solver>
    </solver_null_vecs>
    <solver_coarse>
      <type>eo</type>
      <solver>
        <type>mr</type>
        <tol>1e-1</tol>
        <verbosity>Detailed</verbosity>
        <prefix>c0</prefix>
        <prec>
          <type>dd</type>
          <solver>
            <type>mr</type>
            <max_its>2</max_its>
            <tol>1e-1</tol>
            <error_if_not_converged>false</error_if_not_converged>
            <verbosity>Detailed</verbosity>
            <prefix>c0_dd</prefix>
          </solver>
        </prec>
      </solver>
    </solver_coarse>
    <solver_smoother>
      <type>eo</type>
      <solver>
        <type>mr</type>
        <tol>1e-1</tol>
        <max_its>5</max_its>
        <error_if_not_converged>false</error_if_not_converged>
        <prefix>s0</prefix>
        <prec>
          <type>dd</type>
          <solver>
            <type>mr</type>
            <max_its>4</max_its>
            <tol>1e-1</tol>
            <error_if_not_converged>false</error_if_not_converged>
            <verbosity>Detailed</verbosity>
            <prefix>s0_dd</prefix>
          </solver>
        </prec>
      </solver>
    </solver_smoother>
  </prec_ee>
</InvertParam>
```

Main options:

* `num_null_vecs`: number of null vectors.
* `blocking`: sites factor reduction in each lattice direction for producing the prolongator $V$.
* `solver_null_vecs`: solver to compute the null vectors, approximate solutions of $Dx=0$.
* `solver_coarse`: solver for the coarse operator, $V^* D V$, where $V$ is the prolongator.
* `solver_smoother`: solver for the correction (post-smoothing).

# `mgproton` projection collection

## Deflation projector

If $U$ and $V$ and $\Sigma$ are the smallest singular triplets of $D$, that is $DV=\Sigma U$, then this builds the
following oblique projector, $V(U^*D*V)^{-1}U^*D$.

```
<Projector>
  <projectorType>MGPROTON</projectorType>
  <type>defl</type>
  <rank>200</rank>
  <tol>1e-1</tol>
  <solver>
    <type>bicgstab</type>
    <tol>1e-2</tol>
    <max_its>20000</max_its>
    <verbosity>Detailed</verbosity>
    <prefix>eig</prefix>
  </solver>
  <eigensolver>
    <verbosity>Detailed</verbosity>
  </eigensolver>
</Projector>
```

Main options:
* `rank`: number of singular triplets to compute.
* `tol`: (optional, default 0.1): relative error of the singular triplets relative to $\|D^{-1}\|_2$, $\|\gamma_5 D^{-1} v - \sigma v \|_2 \leq \text{tol}\ \|D^{-1}\|_2$.
* `solver`: solver to estimate $D^{-1}$ used by the eigensolver.
* `eigensolver/max_block_size`: (optional, default is `1`) maximum number of vectors expanding the search subspace in each iteration.
* `eigensolver/max_basis_size`: (optional, default is `PRIMME`'s default) maximum rank of the search subspace.
* `eigensolver/verbosity`: (optional, default `nooutput`) one of `nooutput`, `summary`, `detailed`.

## Multigrid-based deflation projector

If $U$ and $V$ and $\Sigma$ are the smallest singular triplets of $P^*DP$, that is $DV=\Sigma U$, and $P$ is a multigrid prolongator, then this builds the
following oblique projector, $PV(U^*P^*D*PV)^{-1}U^*P^*D$.

```
<?xml version="1.0"?>
<Projector>
  <projectorType>MGPROTON</projectorType>
  <type>mg</type>
  <prolongator>
    <num_null_vecs>24</num_null_vecs>
    <blocking>4 4 4 4</blocking>
    <solver_null_vecs>
      <type>eo</type>
      <solver>
        <type>mr</type>
        <tol>1e-3</tol>
        <max_its>50</max_its>
        <error_if_not_converged>false</error_if_not_converged>
        <prefix>nv0</prefix>
        <verbosity>Detailed</verbosity>
        <prec>
          <type>dd</type>
          <solver>
            <type>mr</type>
            <max_its>2</max_its>
            <tol>1e-1</tol>
            <error_if_not_converged>false</error_if_not_converged>
            <verbosity>Detailed</verbosity>
            <prefix>nv0_dd</prefix>
          </solver>
        </prec>
      </solver>
    </solver_null_vecs>
  </prolongator>
  <rank>200</rank>
  <tol>1e-1</tol>
  <solver>
    <type>bicgstab</type>
    <tol>1e-2</tol>
    <max_its>20000</max_its>
    <verbosity>Detailed</verbosity>
    <prefix>eig</prefix>
  </solver>
  <eigensolver>
    <verbosity>Detailed</verbosity>
  </eigensolver>
</Projector>
```

Main options:
* `rank`: number of singular triplets to compute.
* `tol`: (optional, default 0.1): relative error of the singular triplets relative to $\|D^{-1}\|_2$, $\|\gamma_5 D^{-1} v - \sigma v \|_2 \leq \text{tol}\ \|D^{-1}\|_2$.
* `solver`: solver to estimate $D^{-1}$ used by the eigensolver.
* `eigensolver/max_block_size`: (optional, default is `1`) maximum number of vectors expanding the search subspace in each iteration.
* `eigensolver/max_basis_size`: (optional, default is `PRIMME`'s default) maximum rank of the search subspace.
* `eigensolver/verbosity`: (optional, default `nooutput`) one of `nooutput`, `summary`, `detailed`.
