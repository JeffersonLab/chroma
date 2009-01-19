// -*- C++ -*-
// $Id: inline_static_light_cont_w.h,v 1.2 2009-01-19 20:02:35 caubin Exp $
/*! \file
 * \brief Inline static light contractions for weak three and four point functions
 *
 * Static light contractions for weak three and four point functions
 */

#ifndef __inline_static_light_cont_h__
#define __inline_static_light_cont_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineStaticLightContEnv 
  {
    extern const std::string name;
    bool registerAll();
  }

  //! Parameter structure
  /*! \ingroup inlinehadron */
  struct InlineStaticLightContParams 
  {
    InlineStaticLightContParams();
    InlineStaticLightContParams(XMLReader& xml_in, const std::string& path);
    void write(XMLWriter& xml_out, const std::string& path);

    unsigned long frequency;

    struct Param_t
    {

      bool         Pt_snk;
      bool         Sl_snk;
      bool         Wl_snk;
      string       wvf_kind;
      multi1d<Real>  wvf_param;
      multi1d<int>   wvfIntPar;

      bool MesonP;             // Do Meson Spect.
      bool FourPt;             // Calculate four-point function in addition to 3pt
    } param;

    struct NamedObject_t
    {
      std::string  gauge_id;           /*!< Input gauge field */

      struct Props_t
      {
	std::string  first_id;        
	std::string  second_id;      
	std::string  third_id;      //This is the spect. quark and is ignored if FourPt==false
      };

      multi1d<Props_t> sink_pairs;
    } named_obj;

    std::string xml_file;  // Alternate XML file pattern
  };


  //! Inline measurement of static-light quark spectroscopy
  /*! \ingroup inlinehadron */
  class InlineStaticLightCont : public AbsInlineMeasurement 
  {
  public:
    ~InlineStaticLightCont() {}
    InlineStaticLightCont(const InlineStaticLightContParams& p) : params(p) {}
    InlineStaticLightCont(const InlineStaticLightCont& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out); 

  protected:
    //! Do the measurement
    void func(const unsigned long update_no,
	      XMLWriter& xml_out); 

  private:
    InlineStaticLightContParams params;
  };

};

#endif
