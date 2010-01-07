// -*- C++ -*-
/*! \file
 * \brief Compute the matrix element of  LatticeColorVector*M^-1*Gamma*M^-1**LatticeColorVector
 *
 * Generalized propagator calculation on a colorvector
 */

#ifndef __inline_genprop_matelem_pt_colorvec_h__
#define __inline_genprop_matelem_pt_colorvec_h__

#include "meas/inline/abs_inline_measurement.h"
#include "io/xml_group_reader.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineGenPropMatElemPtColorVecEnv 
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
	bool                    restrict_plateau;       /*!< Restrict the time slices save to the plateau region */
	bool                    avg_equiv_mom;          /*!< Average the genprop over rotations of momenta at fixed mom^2 */
	int                     t_sink;                 /*!< Sink time slice for props */
	int                     mom2_max;               /*!< (mom)^2 <= mom2_max */
	multi1d<int>            mom_offset;             /*!< Momentum origin - the mom2_max will be around here */
	int                     displacement_length;    /*!< Displacement length for creat. and annih. ops */
	int                     num_vecs;               /*!< Number of color vectors to use */
	std::string             mass_label;             /*!< Some kind of mass label */

	struct DispGamma_t
	{
	  int                   gamma;                  /*!< The gamma matrix for this displacement */
	  multi1d<int>          displacement;           /*!< The displacement path for this gamma */
	};

	multi1d<DispGamma_t>    disp_gamma_list;        /*!< Array of displacements and gammas to generate */

	GroupXML_t              link_smearing;          /*!< link smearing xml */
      };

      struct NamedObject_t
      {
	std::string         gauge_id;               /*!< Gauge field */
	std::string         source_prop_id;         /*!< Id for input propagator solutions */
	std::string         sink_prop_id;           /*!< Id for input propagator solutions */
	std::string         genprop_op_file;        /*!< File for generalized propagators operators */
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

  } // namespace InlineGenPropMatElemPtColorVecEnv 
}

#endif
