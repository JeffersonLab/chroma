// -*- C++ -*-
// $Id: inline_heavy_light_cont_w.h,v 3.1 2009-01-20 16:16:58 caubin Exp $
/*! \file
 * \brief Inline heavy light contractions for weak three and four point functions
 *
 * Heavy light contractions for weak three and four point functions
 */

#ifndef __inline_heavy_light_cont_h__
#define __inline_heavy_light_cont_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineHeavyLightContEnv 
  {
    extern const std::string name;
    bool registerAll();
  }

  //! Parameter structure
  /*! \ingroup inlinehadron */
  struct InlineHeavyLightContParams 
  {
    InlineHeavyLightContParams();
    InlineHeavyLightContParams(XMLReader& xml_in, const std::string& path);
    void write(XMLWriter& xml_out, const std::string& path);

    unsigned long frequency;

    struct Param_t
    {
      
      
      bool MesonP;             // Do Meson Spect.
      bool FourPt;             // Calculate four-point function in addition to 3pt
    } param;

    struct NamedObject_t
    {
      std::string  gauge_id;           /*!< Input gauge field */

      struct Props_t
      {
	std::string  first_id;        // Light quark 1
	std::string  second_id;      // Light quark 2
	std::string  third_id;      //This is the spect. quark and is ignored if FourPt==false

	// Following are optional arguments...if not input, static quarks are used.
	std::string heavy_id1;     // Heavy quark 1
	std::string heavy_id2;     // Heavy quark 2
      };

      multi1d<Props_t> sink_pairs;
    } named_obj;

    std::string xml_file;  // Alternate XML file pattern
  };


  //! Inline measurement of heavy-light quark spectroscopy
  /*! \ingroup inlinehadron */
  class InlineHeavyLightCont : public AbsInlineMeasurement 
  {
  public:
    ~InlineHeavyLightCont() {}
    InlineHeavyLightCont(const InlineHeavyLightContParams& p) : params(p) {}
    InlineHeavyLightCont(const InlineHeavyLightCont& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out); 

  protected:
    //! Do the measurement
    void func(const unsigned long update_no,
	      XMLWriter& xml_out); 

  private:
    InlineHeavyLightContParams params;
  };

};

#endif
