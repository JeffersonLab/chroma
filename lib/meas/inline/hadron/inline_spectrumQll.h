// -*- C++ -*-
// $Id: inline_spectrumQll.h,v 1.4 2007-04-18 02:32:26 edwards Exp $
/*! \file
 * \brief Inline spectrum calculations
 *
 * Heavy-light baryon Spectrum calculations (infinitely heavy quark: Qll)
 */

#error "DEPRECATED - TAKEN OUT OF MAKEFILE"


#ifndef __inline_spectrumQll_h__
#define __inline_spectrumQll_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineSpectrumQllEnv 
  {
    extern const std::string name;
    bool registerAll();
  }

  //! Parameter structure
  /*! \ingroup inlinehadron */
  struct InlineSpectrumQllParams 
  {
    InlineSpectrumQllParams();
    InlineSpectrumQllParams(XMLReader& xml_in, const std::string& path);
    void write(XMLWriter& xml_out, const std::string& path);

    unsigned long frequency;

    struct Param_t
    {
      bool Pt_snk;             // point sink
      bool Sl_snk;             // shell sink
      bool Wl_snk;             // wall sink

      WvfKind       wvf_kind;  // Wave function kind: gauge invariant
      multi1d<Real> wvf_param; // Array of width's or other parameters
      //   for "shell" source/sink wave function
      multi1d<int> wvfIntPar;  // Array of iter numbers to approx. Gaussian or
      //   terminate CG inversion for Wuppertal smearing
      Real          link_smear_fact; // smearing factor
      int           link_smear_num;  // number of smearing hits

      multi1d<int> Qsrc_coord ;
    } param;

    struct NamedObject_t
    {
      std::string          gauge_id;
      multi1d<std::string> prop_ids;
    } named_obj;

    std::string xml_file;  // Alternate XML file pattern
  };


  //! Inline measurement of static-light-light spectrum
  /*! \ingroup inlinehadron */
  class InlineSpectrumQll : public AbsInlineMeasurement 
  {
  public:
    ~InlineSpectrumQll() {}
    InlineSpectrumQll(const InlineSpectrumQllParams& p) : params(p) {}
    InlineSpectrumQll(const InlineSpectrumQll& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out); 

  protected:
    //! Do the measurement
    void func(const unsigned long update_no,
	      XMLWriter& xml_out); 

  private:
    InlineSpectrumQllParams params;
  };

};

#endif
