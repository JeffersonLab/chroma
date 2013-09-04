#
#  $Id: regres.pl,v 3.60 2009-07-31 14:09:31 bjoo Exp $
#
#  This is the top-level script used by chroma/scripts/run_chroma_xmldiff.pl
#

sub regresDirs
{
    # 
    # I'd really like to include the regres.pl scripts recursively down 
    # inside the chroma/tests directory, but I'm having difficulty
    # convincing perl to do it. The problem seems to be that the model
    # I want, namely like the c-preprocessor including files that recursively
    # includes other files (in subdirs) is not how the perl "do" works.
    # So, spell out all the many regression dirs and source them individually.
    #
#   return (
#	"$test_dir/chroma/hadron/hadron_contract/regres.pl",
#	"$test_dir/chroma/hadron/colorvec_matelem/regres.pl"
# 	"$test_dir/t_leapfrog/regres.pl",
#     	"$test_dir/chroma/hadron/propagator/regres.pl",
#   );

     return ( 
 	    "$test_dir/chroma/io/szin_write_obj/regres.pl",
 	    "$test_dir/chroma/io/szin_read_obj/regres.pl",
 	    "$test_dir/chroma/io/nersc_write_obj/regres.pl",
 	    "$test_dir/chroma/io/nersc_read_obj/regres.pl",
	    "$test_dir/chroma/io/qio_write_obj/regres.pl",
	    "$test_dir/chroma/io/qio_read_obj/regres.pl",
	    "$test_dir/chroma/io/usqcd_ddpairs_prop/regres.pl",
	    "$test_dir/chroma/gfix/coulgauge/regres.pl",
	    "$test_dir/chroma/glue/gaugestate/regres.pl",
	    "$test_dir/chroma/glue/fuzwilp/regres.pl",
	    "$test_dir/chroma/glue/wilslp/regres.pl",
	    "$test_dir/chroma/glue/qtop_naive/regres.pl",
	    "$test_dir/chroma/eig/regres.pl",
	    "$test_dir/chroma/hadron/fermstate/regres.pl",
	    "$test_dir/chroma/hadron/make_source/regres.pl",
	    "$test_dir/chroma/hadron/propagator/regres.pl",
	    "$test_dir/chroma/hadron/sink_smear/regres.pl",
	    "$test_dir/chroma/hadron/seqsource/regres.pl",
	    "$test_dir/chroma/hadron/building_blocks/regres.pl",
	    "$test_dir/chroma/hadron/noisy_building_blocks/regres.pl",
	    "$test_dir/chroma/hadron/bar3ptfn/regres.pl",
	    "$test_dir/chroma/hadron/gauge_transf/regres.pl",
	    "$test_dir/chroma/hadron/hadspec/regres.pl",
	    "$test_dir/chroma/hadron/mesonspec/regres.pl",
	    "$test_dir/chroma/hadron/hadron_contract/regres.pl",
	    "$test_dir/chroma/hadron/colorvec_matelem/regres.pl",
	    "$test_dir/chroma/hadron/qqq/regres.pl",
	    "$test_dir/chroma/hadron/qqqNucNuc/regres.pl",
	    "$test_dir/chroma/hadron/static_light/regres.pl",
	    "$test_dir/chroma/hadron/stoch_group_meson/regres.pl",
	    "$test_dir/chroma/hadron/stoch_group_baryon/regres.pl",
	    "$test_dir/chroma/hadron_s/make_source/regres.pl",
	    "$test_dir/chroma/hadron_s/propagator/regres.pl",
	    "$test_dir/chroma/hadron_s/sink_smear/regres.pl",
	    "$test_dir/chroma/schrfun/regres.pl",
	    "$test_dir/chroma/smear/link_smear/regres.pl",
	    "$test_dir/chroma/pbp/regres.pl",
	    "$test_dir/t_leapfrog/regres.pl",
	    "$test_dir/spectrum_s/regres.pl",
	    "$test_dir/hmc/regres.pl",
	    "$test_dir/purgaug/regres.pl",
     );
}

# Return a true value for a require statement
1;
