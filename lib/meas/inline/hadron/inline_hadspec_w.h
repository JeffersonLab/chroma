// -*- C++ -*-
// $Id: inline_hadspec_w.h,v 1.2 2005-09-25 20:41:09 edwards Exp $
/*! \file
 * \brief Inline hadron spectrum calculations
 *
 * Hadron spectrum calculations
 */

#ifndef __inline_hadspec_h__
#define __inline_hadspec_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/smearing_io.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineHadSpecEnv 
  {
    extern const std::string name;
    extern const bool registered;
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

      bool HybMesP;            // Hybrid meson spectroscopy
      int  numb_sm;            // number of smearing levels for E- and B-fields
      Real fact_sm;            // Smearing factor for "smeared" E- and B-fields 

      bool time_rev;           // Use time reversal in baryon spectroscopy

      int mom2_max;            // (mom)^2 <= mom2_max. mom2_max=7 in szin.
      bool avg_equiv_mom;      // average over equivalent momenta
      multi1d<int> nrow;
    } param;

    struct NamedObject_t
    {
      multi1d<std::string> prop_ids;
    } named_obj;

    std::string xml_file;  // Alternate XML file pattern
  };


  //! Inline measurement of Wilson loops
  /*! \ingroup inlinehadron */
  class InlineHadSpec : public AbsInlineMeasurement 
  {
  public:
    ~InlineHadSpec() {}
    InlineHadSpec(const InlineHadSpecParams& p) : params(p) {}
    InlineHadSpec(const InlineHadSpec& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const multi1d<LatticeColorMatrix>& u,
		    XMLBufferWriter& gauge_xml,
		    const unsigned long update_no,
		    XMLWriter& xml_out); 

  protected:
    //! Do the measurement
    void func(const multi1d<LatticeColorMatrix>& u,
	      XMLBufferWriter& gauge_xml,
	      const unsigned long update_no,
	      XMLWriter& xml_out); 

  private:
    InlineHadSpecParams params;
  };

};

#endif
