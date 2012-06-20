// -*- C++ -*-
/*! \file
 * \brief Compute the propagator from distillution
 *
 * Propagator calculation in distillution
 */

#ifndef __inline_prop_distillution_w_h__
#define __inline_prop_distillution_w_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlinePropDistillutionEnv 
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
	  int             num_vecs;       /*!< Number of vectors */
	  int             num_vec_dils;   /*!< Number of eigenvector dilutions to use */
	  multi1d<int>    t_sources;      /*!< Array of time slice sources for props */
	  multi1d<int>    quark_lines;    /*!< Quark line numbers */
	  std::string     mass;           /*!< Some kind of mass label */
	  int             Nt_forward;     /*!< Time-slices in the forward direction */
	  int             Nt_backward;    /*!< Time-slices in the backward direction */
	};

	ChromaProp_t    prop;
	Contract_t      contract;
      };

      struct NamedObject_t
      {
	std::string     gauge_id;       /*!< Gauge field */
	std::string     distillution_id;/*!< Distillution factory */
	std::string     colorvec_file;  /*!< Eigenvector file */
	std::string     prop_file;      /*!< Map for output propagator solutions */
      };


      Param_t           param;
      NamedObject_t     named_obj;
      std::string       xml_file;       /*!< Alternate XML file pattern */
    };


    //! Inline task for the propagator from distillution
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

#endif
