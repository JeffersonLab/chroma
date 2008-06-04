// -*- C++ -*-
// $Id: inline_stoch_group_meson_w.h,v 1.3 2008-06-04 03:19:53 edwards Exp $
/*! \file
 * \brief Inline measurement of stochastic group meson operators
 */

#ifndef __inline_stoch_group_meson_h__
#define __inline_stoch_group_meson_h__

#include "meas/inline/abs_inline_measurement.h"
#include "io/xml_group_reader.h"
#include "io/enum_io/enum_mesonop_io.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineStochGroupMesonEnv 
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
	MesonOpType         creat_op_contract_type; /*<! Contraction type for creation op */
	MesonOpType         annih_op_contract_type; /*<! Contraction type for annihilation op */
	int                 mom2_max;               /*!< (mom)^2 <= mom2_max */
	int                 displacement_length;    /*!< Displacement length for creat. and annih. ops */
	GroupXML_t          quark_smearing;         /*!< xml string holding smearing params */
	GroupXML_t          link_smearing;          /*!< link smearing xml */

	multi1d<GroupXML_t> quark_dils;             /*!< Dilutions for each quark */
      };

      struct NamedObject_t
      {
	struct TwoQuarkOpsFile_t
	{
	  std::string        ops_file;       /*!< Coefficient file name */
	  std::string        id;             /*!< ID/tag used in analysis codes*/
	};

	TwoQuarkOpsFile_t    operators_file; /*!< Files holding 2-quark ops to make*/
	std::string          quark_ids;      /*!< 2 character string indicating which quarks are degenerate */
	std::string          gauge_id;       /*!< Gauge field */
      };

      Param_t        param;      /*!< Parameters */    
      NamedObject_t  named_obj;  /*!< Named objects */
      std::string    xml_file;   /*!< Alternate XML file pattern */
    };


    //! Inline measurement of stochastic group meson operators
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

  } // namespace InlineStochGroupMesonEnv 
}

#endif
