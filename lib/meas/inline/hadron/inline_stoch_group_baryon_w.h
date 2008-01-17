// -*- C++ -*-
// $Id: inline_stoch_group_baryon_w.h,v 1.5 2008-01-17 15:09:26 jbulava Exp $
/*! \file
 * \brief Inline measurement of stochastic group baryon operator
 */

#ifndef __inline_stoch_group_baryon_h__
#define __inline_stoch_group_baryon_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"

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
				int                 mom2_max;              /*!< (mom)^2 <= mom2_max */
				int                 displacement_length;   /*!< Displacement length for creat. and annih. ops */
				GroupXML_t          source_quark_smearing; /*!< xml string holding smearing params */
				GroupXML_t          sink_quark_smearing;   /*!< xml string holding smearing params */
				GroupXML_t          link_smearing;         /*!< link smearing xml */

				multi1d<GroupXML_t> quark_dils;             /*!< Dilutions for each quark */
      } param;

      struct NamedObject_t
      {
				struct ThreeQuarkOpsFile_t
				{
	  			std::string        ops_file;       /*!< Coefficient file name */
	  			std::string        id;             /*!< ID/tag used in analysis codes*/
				};

				ThreeQuarkOpsFile_t  operators_file; /*!< Files holding 3-quark ops to make*/

				std::string          quark_ids;      /*!< 3 character string indicating which quarks are degenerate */
	
				std::string          gauge_id;       /*!< Gauge field */
      
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
