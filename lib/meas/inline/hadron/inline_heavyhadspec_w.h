// -*- C++ -*-
// $Id: inline_heavyhadspec_w.h,v 1.1 2008-07-21 18:15:36 kostas Exp $
/*! \file
 * \brief Inline hadron spectrum calculations
 *
 * Hadron spectrum calculations
 */

#ifndef __inline_heavyhadspec_h__
#define __inline_heavyhadspec_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineHeavyHadSpecEnv 
  {
    extern const std::string name;
    bool registerAll();
  }

  //! Parameter structure
  /*! \ingroup inlinehadron */
  struct InlineHeavyHadSpecParams 
  {
    InlineHeavyHadSpecParams();
    InlineHeavyHadSpecParams(XMLReader& xml_in, const std::string& path);
    void write(XMLWriter& xml_out, const std::string& path);

    unsigned long frequency;

    struct Param_t
    {
      bool time_rev;           // Use time reversal in baryon spectroscopy

      int mom2_max;            // (mom)^2 <= mom2_max. mom2_max=7 in szin.
      bool avg_equiv_mom;      // average over equivalent momenta
    } param;

    struct NamedObject_t
    {
      std::string  gauge_id;           /*!< Input gauge field */

      struct Props_t
      {
	std::string  first_id;
	std::string  second_id;
      };

      multi1d<Props_t> sink_pairs;
    } named_obj;

    std::string xml_file;  // Alternate XML file pattern
  };


  //! Inline measurement of hadron spectrum
  /*! \ingroup inlinehadron */
  class InlineHeavyHadSpec : public AbsInlineMeasurement 
  {
  public:
    ~InlineHeavyHadSpec() {}
    InlineHeavyHadSpec(const InlineHeavyHadSpecParams& p) : params(p) {}
    InlineHeavyHadSpec(const InlineHeavyHadSpec& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out); 

  protected:
    //! Do the measurement
    void func(const unsigned long update_no,
	      XMLWriter& xml_out); 

  private:
    InlineHeavyHadSpecParams params;
  };

};

#endif
