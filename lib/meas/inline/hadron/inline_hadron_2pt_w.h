// -*- C++ -*-
// $Id: inline_hadron_2pt_w.h,v 3.1 2007-04-18 02:32:26 edwards Exp $
/*! \file
 * \brief Inline hadron spectrum calculations
 *
 * Hadron spectrum calculations. The general version that write output
 * into lime files.
 */

#ifndef __inline_hadron_2pt_h__
#define __inline_hadron_2pt_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineHadron2PtEnv 
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
	bool          time_rev;           /*!< Use time reversal in baryon spectroscopy */
	int           mom2_max;           /*!< (mom - mom_origin)^2 <= mom2_max */
	multi1d<int>  mom_origin;         /*!< Origin for the momentum */
	bool          avg_equiv_mom;      /*!< average over equivalent momenta */
      };

      struct NamedObject_t
      {
	std::string  gauge_id;            /*!< Input gauge field */

	struct Correlators_t
	{
	  GroupXML_t             corr_xml;   /*!< xml string holding source spin insertion params */
	  multi1d<std::string>   prop_ids;   /*!< All the prop ids needed for this correlator */
	};

	multi1d<Correlators_t> correlators;  /*!< The correlatords */
      };

      unsigned long   frequency;

      Param_t         param;       /*!< Params common to all correlators */
      NamedObject_t   named_obj;   /*!< Named objects */
      std::string     data_file;   /*!< Output file for data in QIO/LIME format */
      std::string     xml_file;    /*!< Alternate XML file pattern */
    };


    //! Inline measurement of hadron 2pt functions
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
