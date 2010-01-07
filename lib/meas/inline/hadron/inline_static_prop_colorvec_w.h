// -*- C++ -*-
/*! \file
 * \brief Compute a static prop  (1/2)*(1+gamma_4)U*U*...U * multi1d<LatticeColorVector>
 *
 * Static propagator calculation on a colorvector
 */

#ifndef __inline_static_prop_colorvec_w_h__
#define __inline_static_prop_colorvec_w_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"
#include "io/xml_group_reader.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineStaticPropColorVecEnv 
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
	  int num_vecs;             /*!< Number of color vectors to use */
	  int decay_dir;            /*!< Decay direction */
	  multi1d<int> t_sources;   /*!< Array of time slice sources for static props */
	};

	Contract_t      contract;
      } param;
      
      GroupXML_t        map_obj_params; /*!< Parameters for MapObj factory */
      struct NamedObject_t
      {
	std::string     gauge_id;       /*!< Gauge field */
	std::string     colorvec_id;    /*!< LatticeColorVector EigenInfo */
	std::string     prop_id;        /*!< Id for output propagator solutions */
	GroupXML_t      prop_obj;       /*!< Map for output propagator solutions */
      } named_obj;

      std::string xml_file;  // Alternate XML file pattern
    };


    //! Inline task for compute LatticeColorVector matrix elements of a static propagator
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
