// -*- C++ -*-
// $Id: inline_spectrum_s.h,v 1.2 2005-08-24 16:39:23 mcneile Exp $
/*! \file
 * \brief Inline staggered spectrum calculations
 *
 * Staggered spectrum calculations
 */

#ifndef __inline_spectrum_s_h__
#define __inline_spectrum_s_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "invtype.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineSpectrumEnv 
  {
    extern const std::string name;
    extern const bool registered;
  }

  //! Parameter structure
  /*! \ingroup inlinehadron */
  struct InlineSpectrumParams_s 
  {
    InlineSpectrumParams_s();
    InlineSpectrumParams_s(XMLReader& xml_in, const std::string& path);
    void write(XMLWriter& xml_out, const std::string& path);

    unsigned long frequency;

    struct Param_t
    {
      bool Meson_local;             // Meson spectroscopy
      bool Baryon_local;            // Baryons spectroscopy
      bool disconnected_local;            // disconnected loops

      multi1d<int> boundary;
      multi1d<int> nrow;
      multi1d<int> t_srce;

    } param;

    struct Quark_Prop_t
    {
      Real         Mass;      // Staggered mass
      Real         u0;        // Tadpole Factor

      InvertParam_t  invParam; // inversion parameters

    } prop_param ;

    std::string xml_file;  // Alternate XML file pattern
  };


  //! Inline measurement of Wilson loops
  /*! \ingroup inlinehadron */
  class InlineSpectrum_s : public AbsInlineMeasurement 
  {
  public:
    ~InlineSpectrum_s() {}
    InlineSpectrum_s(const InlineSpectrumParams_s& p) : params(p) {}
    InlineSpectrum_s(const InlineSpectrum_s& p) : params(p.params) {}

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
    InlineSpectrumParams_s params;
  };

};

#endif
