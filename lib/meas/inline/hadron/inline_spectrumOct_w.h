// -*- C++ -*-
// $Id: inline_spectrumOct_w.h,v 3.0 2006-04-03 04:59:03 edwards Exp $
/*! \file
 * \brief Inline heavy light spectrum calculations
 *
 * Octet baryons and mesons spectrum calculations 
 */

#ifndef __inline_spectrumOct_h__
#define __inline_spectrumOct_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/smearing_io.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineSpectrumOctEnv 
  {
    extern const std::string name;
    extern const bool registered;
  }

  //! Parameter structure
  /*! \ingroup inlinehadron */
  struct InlineSpectrumOctParams 
  {
    InlineSpectrumOctParams();
    InlineSpectrumOctParams(XMLReader& xml_in, const std::string& path);
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

    } param;

    struct NamedObject_t
    {
      std::string          gauge_id;
      multi1d<std::string> prop_ids; 
    } named_obj;

    std::string xml_file;  // Alternate XML file pattern
  };


  //! Inline measurement of Wilson loops
  /*! \ingroup inlinehadron */
  class InlineSpectrumOct : public AbsInlineMeasurement 
  {
  public:
    ~InlineSpectrumOct() {}
    InlineSpectrumOct(const InlineSpectrumOctParams& p) : params(p) {}
    InlineSpectrumOct(const InlineSpectrumOct& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out); 

  protected:
    //! Do the measurement
    void func(const unsigned long update_no,
	      XMLWriter& xml_out); 

  private:
    InlineSpectrumOctParams params;
  };

};

#endif
