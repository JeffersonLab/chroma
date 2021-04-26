// -*- C++ -*-
/*! \file
 * \brief Inline measurement that construct unsmeared hadron nodes using distillation
 */

#ifndef __inline_unsmeared_hadron_node_distillation_superb_h__
#define __inline_unsmeared_hadron_node_distillation_superb_h__

#include "chromabase.h"

#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"
#include "io/xml_group_reader.h"

#ifdef BUILD_SB

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineUnsmearedHadronNodeDistillationSuperbEnv 
  {
    bool registerAll();


    //! Parameter structure
    /*! \ingroup inlinehadron */
    struct Params
    {
      Params();
      Params(XMLReader& xml_in, const std::string& path);

      unsigned long     frequency;

      //! Parameters
      struct Param_t
      {
	GroupXML_t                  link_smearing;          /*!< link smearing xml */

	struct DispGammaMom_t
	{
	  int                       gamma;                  /*!< The gamma matrix for this displacement */
	  std::vector<int>          displacement;           /*!< The displacement path for this gamma */
	  multi1d<int>              mom;                    /*!< Array of momenta to generate */
	};

	struct KeySolnProp_t
	{
	  bool                      cacheP;
	  int                       num_vecs;               /*!< Number of color vectors to use */
	  int                       t_source;               /*!< Time slice source for props */
	};

	struct SinkSource_t
	{
	  int                       t_sink;                 /*!< Time slice for sinks */
	  int                       t_source;               /*!< Time slice source for props */
	  int                       Nt_backward;            /*!< Backward relative to source */
 	  int                       Nt_forward;             /*!< Forward relative to source */
	};
	
	struct Contract_t
	{
	  int                       alt_t_start;            /*!< Alternative Starting time-slice for genprops */
 	  int                       alt_Nt_forward;         /*!< Alternative Forward relative to t_start */
 	  int                       alt_num_vecs;           /*!< Number of color vectors to use */
	  bool                      use_derivP;             /*!< Use derivatives */
	  int                       decay_dir;              /*!< Decay direction */
	  int                       displacement_length;    /*!< Displacement length for insertions */
	  std::string               mass_label;             /*!< Some kind of mass label */
	  int                       num_tries;              /*!< In case of bad things happening in the solution vectors, do retries */
	  int                       max_rhs;                /*! maximum number of linear systems solved simultaneously */
	  int                       max_tslices_in_contraction;  /*! maximum number of contracted tslices simultaneously */
	  int                       max_moms_in_contraction;/*! maximum number of contracted momenta simultaneously */
	  bool                      use_genprop4_format;    /*!< Use the efficient genprop4 format instead of the traditional one */
	  bool                      use_multiple_writers;   /*!< Whether several processes are going to write down the elementals on separate files */
	  multi1d<float>            phase;                   /*!< Phase to apply to colorvecs */
	};

	std::vector<KeySolnProp_t>  prop_sources;           /*!< Sources */
	std::vector<SinkSource_t>   sink_source_pairs;      /*!< Combos */
	std::vector<DispGammaMom_t> disp_gamma_mom_list;    /*!< Array of displacements, gammas, and moms to generate */	
	std::map<int, std::list<int>> alt_sink_sources;     /*!< Alternative source-sink list {t_source -> list[t_sinks]} */
	std::vector< std::vector<int>> alt_displacements;   /*!< Alternative displacement paths */
	std::vector< multi1d<int>>  alt_moms;               /*!< Alternative array of momenta to generate */
	ChromaProp_t                prop;                   /*!< Propagator input */
	Contract_t                  contract;               /*!< Backward propagator and contraction pieces */
      };

      struct NamedObject_t
      {
 	std::string                 gauge_id;               /*!< Gauge field */
	std::vector<std::string>    colorvec_files;         /*!< Eigenvectors in mod format */
	std::string                 dist_op_file;           /*!< File name for propagator matrix elements */
      };

      Param_t                       param;                  /*!< Parameters */    
      NamedObject_t                 named_obj;              /*!< Named objects */
      std::string                   xml_file;               /*!< Alternate XML file pattern */
    };


    //! Inline measurement that construct hadron nodes using distillution
    /*! \ingroup inlinehadron */
    class InlineMeas : public AbsInlineMeasurement 
    {
    public:
      ~InlineMeas() {}
      InlineMeas(const Params& p) : params(p) {}
      InlineMeas(const InlineMeas& p) : params(p.params) {}

      unsigned long getFrequency(void) const {return params.frequency;}

      //! Do the measurement
      void operator()(const unsigned long update_no,
		      XMLWriter& xml_out); 

    protected:
      //! Do the measurement
      void func(const unsigned long update_no,
		XMLWriter& xml_out); 

    private:
      Params params;
    };

  } // namespace InlineUnsmearedHadronNodeDistillationSuperbEnv
}

#endif // BUILD_SB
#endif
