// -*- C++ -*-
// $Id: inline_stoch_group_baryon_w.h,v 1.1 2007-06-18 19:40:03 edwards Exp $
/*! \file
 * \brief Inline measurement of stochastic group baryon operator
 */

#ifndef __inline_stoch_group_baryon_h__
#define __inline_stoch_group_baryon_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineStochGroupBaryonEnv 
  {
    extern const std::string name;
    bool registerAll();


    //! Parameter structure
    /*! \ingroup inlinehadron */
    struct Params
    {
      Params();
      Params(XMLReader& xml_in, const std::string& path);
      void writeXML(XMLWriter& xml_out, const std::string& path);

      unsigned long      frequency;

      struct Param_t
      {
	int              mom2_max;               /*!< (mom)^2 <= mom2_max */
	int              displacement_length;    /*!< Displacement length for creat. and annih. ops */
	GroupXML_t       source_quark_smearing;  /*!< xml string holding smearing params */
	GroupXML_t       sink_quark_smearing;    /*!< xml string holding smearing params */
	GroupXML_t       link_smearing;          /*!< link smearing xml */
      } param;

      struct NamedObject_t
      {
	multi1d<std::string>  operator_coeff_files;   /*!< Files holding group coefficients */

	//! Operators
	struct Operator_t
	{
	  multi1d<std::string> soln_files;
	};

	multi1d<Operator_t>   op;                     /*!< All the quarks and their solutions that are needed */

	std::string           gauge_id;
	std::string           operator_file;
      } named_obj;

      std::string xml_file;  // Alternate XML file pattern
    };


    //! Inline measurement of stochastic group baryon operators
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

  } // namespace InlineStochGroupBaryonEnv 
}

#endif
