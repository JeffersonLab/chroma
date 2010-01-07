// -*- C++ -*-
/*! \file
 * \brief Compute the perambulators
 * 
 * Compute propagator elements    M^-1 * multi1d<LatticeColorVector>
 * and the matrix element of  LatticeColorVector*M^-1*LatticeColorVector
 */

#ifndef __inline_prop_and_matelem_colorvec_w_h__
#define __inline_prop_and_matelem_colorvec_w_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"
#include "io/xml_group_reader.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlinePropAndMatElemColorVecEnv 
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
	  std::string   mass_label;     /*!< Some kind of mass label */
	};

	ChromaProp_t    prop;
	Contract_t      contract;
      };

      struct NamedObject_t
      {
	std::string     gauge_id;       /*!< Gauge field */
	std::string     colorvec_id;    /*!< LatticeColorVector EigenInfo */
	std::string     prop_id;        /*!< Output propagator solution id */
	GroupXML_t      prop_obj;       /*!< Map for output propagator solutions */
	std::string     prop_op_file;   /*!< File name for propagator matrix elements */
      };

      Param_t           param;
      NamedObject_t     named_obj;
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
