// -*- C++ -*-
/*! \file
 * \brief Inline static light spectrum calculations
 *
 * Static light spectrum calculations
 */

#ifndef __inline_static_light_spec_h__
#define __inline_static_light_spec_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineStaticLightSpecEnv 
  {
    extern const std::string name;
    bool registerAll();
  }

  //! Parameter structure
  /*! \ingroup inlinehadron */
  struct InlineStaticLightSpecParams 
  {
    InlineStaticLightSpecParams();
    InlineStaticLightSpecParams(XMLReader& xml_in, const std::string& path);
    void write(XMLWriter& xml_out, const std::string& path);

    unsigned long frequency;

    struct Param_t
    {
      bool MesonP;             // Meson spectroscopy
      bool BaryonP;            // Baryons spectroscopy
      bool MesonPot;           // Meson potential 
      bool BaryonPot;          // Baryon potential 
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


  //! Inline measurement of static-light quark spectroscopy
  /*! \ingroup inlinehadron */
  class InlineStaticLightSpec : public AbsInlineMeasurement 
  {
  public:
    ~InlineStaticLightSpec() {}
    InlineStaticLightSpec(const InlineStaticLightSpecParams& p) : params(p) {}
    InlineStaticLightSpec(const InlineStaticLightSpec& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out); 

  protected:
    //! Do the measurement
    void func(const unsigned long update_no,
	      XMLWriter& xml_out); 

  private:
    InlineStaticLightSpecParams params;
  };

};

#endif
