// -*- C++ -*-
/*! \file
 * \brief Compute the propagator from distillation
 *
 * Propagator calculation in distillation
 */

#ifndef __inline_prop_distillation_allspin_w_h__
#define __inline_prop_distillation_allspin_w_h__

#ifndef QDP_IS_QDPJIT_NO_NVPTX

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"
#include "io/xml_group_reader.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlinePropDistillationAllSpinEnv 
  {
    bool registerAll();

    //! Parameter structure
    /*! \ingroup inlinehadron */ 
    struct Params 
    {
      Params();
      Params(XMLReader& xml_in, const std::string& path);

      unsigned long     frequency;

      struct Param_t
      {
	struct Contract_t
	{
	  int           num_vecs;       /*!< Number of color vectors to use */
	  int           decay_dir;      /*!< Decay direction */
	  multi1d<int>  t_sources;      /*!< Array of time slice sources for props */
	  int           Nt_forward;     /*!< Time-slices in the forward direction */
	  int           Nt_backward;    /*!< Time-slices in the backward direction */
	  std::string   mass_label;     /*!< Some kind of mass label */

	  int           num_tries;      /*!< In case of bad things happening in the solution vectors, do retries */
	  bool          zero_colorvecs;
	  bool          fuse_timeloop;  
	};

	ChromaProp_t    prop;
	Contract_t      contract;
      };

      struct NamedObject_t
      {
	std::string                 gauge_id;        /*!< Gauge field */
	std::vector<std::string>    colorvec_files;  /*!< Eigenvectors in mod format */
	std::string                 soln_file;       /*!< File name for propagator matrix elements */
      };

      Param_t           param;
      NamedObject_t     named_obj;
      std::string       xml_file;       /*!< Alternate XML file pattern */
    };


    //! Inline task for the propagator from distillation
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

  } // namespace PropColorVec


}

#endif // QDPJIT

#endif // Chroma
