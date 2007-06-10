// -*- C++ -*-
// $Id: inline_hadron_contract.h,v 3.1 2007-06-10 14:40:23 edwards Exp $
/*! \file
 * \brief Inline hadron contractions - for correlators
 *
 * Hadron spectrum calculations. The general version that write output
 * into lime files.
 */

#ifndef __inline_hadron_contract_h__
#define __inline_hadron_contract_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/xml_group_reader.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineHadronContractEnv 
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

      struct Param_t
      {
	int           mom2_max;           /*!< (mom - mom_origin)^2 <= mom2_max */
	multi1d<int>  mom_origin;         /*!< Origin for the momentum */
	bool          avg_equiv_mom;      /*!< average over equivalent momenta */
      };

      struct NamedObject_t
      {
	std::string  gauge_id;            /*!< Input gauge field */
	std::string  output_file;         /*!< Output file for data in QIO/LIME format */

	struct Correlator_t
	{
	  GroupXML_t xml;                 /*!< XML used for each correlator */
	};

	multi1d<Correlator_t> correlators;  /*!< The correlators */
      };

      unsigned long   frequency;

      Param_t         param;       /*!< Params common to all correlators */
      NamedObject_t   named_obj;   /*!< Named objects */
      std::string     xml_file;    /*!< Alternate XML file pattern */
    };


    //! Inline measurement of hadron contraction functions
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

  }

}

#endif
