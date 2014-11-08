// -*- C++ -*-
// 
/*! \file
 * \brief Compute the matrix element of  LatticeColorVector*M^-1*Gamma*displace*M^-1**LatticeColorVector
 *
 * Generalized propagator calculation on a colorstd::vector for meson distribution amplitudes and parton dstribution functions
 */

#ifndef __inline_genprop_matelem_da_colorvec_h__
#define __inline_genprop_matelem_da_colorvec_h__

#include "meas/inline/abs_inline_measurement.h"
#include "io/xml_group_reader.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineGenPropMatElemDAColorVecEnv 
  {
    bool registerAll();

    //! Parameter structure
    /*! \ingroup inlinehadron */
    struct Params
    {
      Params();
      Params(XMLReader& xml_in, const std::string& path);
      void writeXML(XMLWriter& xml_out, const std::string& path) const;

      unsigned long      frequency;

      struct Param_t
      {
	int          t_source   ; /*!< Source time slice for props            */
	int          t_sink     ; /*!< Sink time slice for props              */
	int          mom_max    ; /*!< (mom) <= mom_max                       */
	int          boost_dir  ; /*!< Direction of the boost                 */
	int          num_vecs   ; /*!< Number of color vectors to use         */
	int          decay_dir  ; /*!< Decay direction                        */
	std::string  mass_label ; /*!< Some kind of mass label                */
	int          gamma      ; /*!< The gamma matrix for this displacement */
	bool         restrict_plateau;
      };

      struct NamedObject_t
      {
	std::string gauge_id;         /*!< Gauge field                       */
	std::string sm_gauge_id;      /*!< smeared Gauge field               */
	std::string source_prop_id;   /*!< Id for input propagator solutions */
	std::string sink_prop_id;     /*!< Id for input propagator solutions */
	std::string genprop_op_file;  /*!< File for generalized propagators operators */
      };

      Param_t        param;      /*!< Parameters */    
      NamedObject_t  named_obj;  /*!< Named objects */
      std::string    xml_file;   /*!< Alternate XML file pattern */
    };


    //! Compute the matrix element of  LatticeColorVector*M^-1*Gamma*M^-1**LatticeColorVector
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

  } // namespace InlineGenPropMatElemDAColorVecEnv 
}

#endif
