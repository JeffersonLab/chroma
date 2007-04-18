// -*- C++ -*-
// $Id: inline_hadspec_w.h,v 3.5 2007-04-18 02:32:26 edwards Exp $
/*! \file
 * \brief Inline hadron spectrum calculations
 *
 * Hadron spectrum calculations
 */

#ifndef __inline_hadspec_h__
#define __inline_hadspec_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineHadSpecEnv 
  {
    extern const std::string name;
    bool registerAll();
  }

  //! Parameter structure
  /*! \ingroup inlinehadron */
  struct InlineHadSpecParams 
  {
    InlineHadSpecParams();
    InlineHadSpecParams(XMLReader& xml_in, const std::string& path);
    void write(XMLWriter& xml_out, const std::string& path);

    unsigned long frequency;

    struct Param_t
    {
      bool MesonP;             // Meson spectroscopy
      bool CurrentP;           // Meson currents
      bool BaryonP;            // Baryons spectroscopy

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
  class InlineHadSpec : public AbsInlineMeasurement 
  {
  public:
    ~InlineHadSpec() {}
    InlineHadSpec(const InlineHadSpecParams& p) : params(p) {}
    InlineHadSpec(const InlineHadSpec& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out); 

  protected:
    //! Do the measurement
    void func(const unsigned long update_no,
	      XMLWriter& xml_out); 

  private:
    InlineHadSpecParams params;
  };

};

#endif
