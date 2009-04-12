// -*- C++ -*-
// $Id: inline_annih_prop_matelem_colorvec_w.h,v 3.2 2009-04-12 22:06:13 edwards Exp $
/*! \file
 * \brief Compute the annihilation diagram propagator elements    M^-1 * multi1d<LatticeColorVector>
 *
 * Annihilation diagrams version of propagator calculation on a colorvector
 */

#ifndef __inline_annih_prop_matelem_colorvec_w_h__
#define __inline_annih_prop_matelem_colorvec_w_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineAnnihPropMatElemColorVecEnv 
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
	  int           dt;             /*!< Separation of sources */
	  multi1d<int>  t_sources_start; /*!< Starting location for diluted time sources */
	  int           N;              /*!< N in Z(N) source */
	  Seed          ran_seed;       /*!< RNG seed for this source */
	  std::string   mass_label;     /*!< Some kind of mass label */
	};

	Contract_t      contract;       /*!< Contraction parameters*/
	ChromaProp_t    prop;           /*!< Parameters for prop solutions */
      };

      struct NamedObject_t
      {
	std::string     gauge_id;       /*!< Gauge field */
	std::string     colorvec_id;    /*!< LatticeColorVector EigenInfo */
	std::string     prop_op_file;   /*!< File name for propagator matrix elements */
      };

      Param_t           param;          /*!< Parameters */
      NamedObject_t     named_obj;      /*!< Input and output objects */
      std::string       xml_file;       /*!< Alternate XML file pattern */
    };


    //! Inline task for compute LatticeColorVector matrix elements of a propagator
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
