// -*- C++ -*-
// $Id: inline_spectrum_s.h,v 2.1 2006-02-02 16:23:14 egregory Exp $
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
#include "meas/hadron/enum_loops_s.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineStaggeredSpectrumEnv 
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
      bool Meson_local;                 // Meson spectroscopy
      bool Baryon_local;                // Baryons spectroscopy
      bool Baryon_vary;                 // Baryons spectroscopy variational
      bool LocalPion_vary;              // local pion spectroscopy variational
      bool disconnected_local ;         // disconnected loops local
      bool disconnected_fuzz  ;         // disconnected loops fuzz
      bool ps4link_singlet_conn ;  
      bool eight_pions;                 // all pseudoscalar meson tastes
      bool eight_scalars;               // all scalar meson tastes
      bool eight_rhos;                  // all vector meson tastes

 
      // parameters for disconnected loops
      int Nsamp;
      int CFGNO ;
      VolSrc_type volume_source ;
      bool gauge_invar_oper ; 
      bool sym_shift_oper ; 

      int fuzz_width ; 

      int src_seperation ;

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
