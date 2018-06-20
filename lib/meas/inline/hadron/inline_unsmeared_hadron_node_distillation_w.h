// -*- C++ -*-
/*! \file
 * \brief Inline measurement that construct unsmeared hadron nodes using distillation
 */

#ifndef __inline_unsmeared_hadron_node_distillation_h__
#define __inline_unsmeared_hadron_node_distillation_h__

#ifndef QDP_IS_QDPJIT_NO_NVPTX

#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"
#include "io/xml_group_reader.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineUnsmsearedHadronNodeDistillationEnv 
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

	struct Contract_t
	{
	  bool                      use_derivP;             /*!< Use derivatives */
	  int                       srce_num_vecs;          /*!< Number of color vectors to use */
	  int                       sink_num_vecs;          /*!< Number of color vectors to use */
	  int                       t_source;               /*!< Time slice source for props */
	  int                       t_sink;                 /*!< Time slice sink for props */
	  int                       decay_dir;              /*!< Decay direction */
	  int                       displacement_length;    /*!< Displacement length for insertions */
	  std::string               mass_label;             /*!< Some kind of mass label */
	  int                       num_tries;              /*!< In case of bad things happening in the solution vectors, do retries */
	};

	std::vector<DispGammaMom_t> disp_gamma_mom_list;    /*!< Array of displacements, gammas, and moms to generate */	
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

  } // namespace InlineUnsmsearedHadronNodeDistillationEnv
}
#endif
#endif
