// -*- C++ -*-
// $Id: inline_spectrum_w.h,v 1.1 2005-04-06 04:34:54 edwards Exp $
/*! \file
 * \brief Inline spectrum calculations
 *
 * Spectrum calculations
 */

#ifndef __inline_spectrum_h__
#define __inline_spectrum_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/smearing_io.h"

namespace Chroma 
{ 
  namespace InlineSpectrumEnv 
  {
    extern const std::string name;
    extern const bool registered;
  }

  //! Parameter structure
  struct InlineSpectrumParams 
  {
    InlineSpectrumParams();
    InlineSpectrumParams(XMLReader& xml_in, const std::string& path);
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

      multi1d<int> nrow;
    } param;

    struct Prop_t
    {
      multi1d<std::string> prop_files;  // The files are expected to be in SciDAC format!
    } prop;
  };


  //! Inline measurement of Wilson loops
  class InlineSpectrum : public AbsInlineMeasurement 
  {
  public:
    ~InlineSpectrum() {}
    InlineSpectrum(const InlineSpectrumParams& p) : params(p) {}
    InlineSpectrum(const InlineSpectrum& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const multi1d<LatticeColorMatrix>& u,
		    XMLBufferWriter& gauge_xml,
		    const unsigned long update_no,
		    XMLWriter& xml_out); 

  private:
    InlineSpectrumParams params;
  };

};

#endif
