// -*- C++ -*-
// $Id: inline_spectrum_s.h,v 1.1 2005-08-23 17:29:55 mcneile Exp $
/*! \file
 * \brief Inline staggered spectrum calculations
 *
 * Staggered spectrum calculations
 */

#ifndef __inline_spectrum_s_h__
#define __inline_spectrum_s_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/smearing_io.h"

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
      bool Pt_snk;             // point sink
      bool Sl_snk;             // shell sink
      bool Wl_snk;             // wall sink

      bool MesonP;             // Meson spectroscopy
      bool CurrentP;           // Meson currents
      bool BaryonP;            // Baryons spectroscopy

      bool HybMesP;            // Hybrid meson spectroscopy
      int  numb_sm;            // number of smearing levels for E- and B-fields
      Real fact_sm;            // Smearing factor for "smeared" E- and B-fields 

      bool time_rev;           // Use time reversal in baryon spectroscopy

      int mom2_max;            // (mom)^2 <= mom2_max. mom2_max=7 in szin.
      bool avg_equiv_mom;      // average over equivalent momenta
      WvfKind       wvf_kind;  // Wave function kind: gauge invariant
      multi1d<Real> wvf_param; // Array of width's or other parameters
      //   for "shell" source/sink wave function
      multi1d<int> wvfIntPar;  // Array of iter numbers to approx. Gaussian or
      //   terminate CG inversion for Wuppertal smearing
      Real          link_smear_fact; // smearing factor
      int           link_smear_num;  // number of smearing hits

      multi1d<int> nrow;
    } param;

    struct Prop_t
    {
      multi1d<std::string> prop_files;  // The files are expected to be in SciDAC format!
    } prop;

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
