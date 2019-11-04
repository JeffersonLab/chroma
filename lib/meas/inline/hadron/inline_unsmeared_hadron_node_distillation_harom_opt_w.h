// -*- C++ -*-
/*! \file
 * \brief Inline measurement that construct unsmeared hadron nodes using distillation
 */

#ifndef __inline_unsmeared_hadron_node_distillation_harom_opt_h__
#define __inline_unsmeared_hadron_node_distillation_harom_opt_h__

#include "chromabase.h"

#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"
#include "io/xml_group_reader.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineUnsmearedHadronNodeDistillationHaromOptEnv 
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
	struct Contract_t
	{
	  int                       num_vecs;               /*!< Number of color vectors to use */
	  int                       t_start;                /*!< Starting time-slice for genprops */
 	  int                       Nt_forward;             /*!< Forward relative to t_start */
	  int                       decay_dir;              /*!< Decay direction */
	  int                       displacement_length;    /*!< Displacement length for insertions */
	  std::string               mass_label;             /*!< Some kind of mass label */
	  int                       num_tries;              /*!< In case of bad things happening in the solution vectors, do retries */
	};

	std::vector<int>            prop_sources;           /*!< Sources */
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

  } // namespace InlineUnsmearedHadronNodeDistillationEnv
}

#endif
