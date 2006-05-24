// -*- C++ -*-
// $Id: inline_mesonspec_w.h,v 3.1 2006-05-24 21:09:41 edwards Exp $
/*! \file
 * \brief Inline meson spectrum calculations
 *
 * Meson spectrum calculations
 */

#ifndef __inline_mesonspec_h__
#define __inline_mesonspec_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/xml_group_reader.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineMesonSpecEnv 
  {
    extern const std::string name;
    extern const bool registered;
  }

  //! Parameter structure
  /*! \ingroup inlinehadron */
  struct InlineMesonSpecParams 
  {
    InlineMesonSpecParams();
    InlineMesonSpecParams(XMLReader& xml_in, const std::string& path);
    void write(XMLWriter& xml_out, const std::string& path);

    unsigned long frequency;

    struct Param_t
    {
      int mom2_max;            /*!< (mom)^2 <= mom2_max. mom2_max=7 in szin. */
      bool avg_equiv_mom;      /*!< average over equivalent momenta */
    } param;

    struct NamedObject_t
    {
      std::string    gauge_id; 

      struct Props_t
      {
	struct Sinks_t
	{
	  std::string   first_id;                    /*!< first quark */
	  std::string   second_id;                   /*!< second quark (this one is adjointed) */

	  GroupXML_t    source_spin_insertion;       /*!< xml string holding source spin insertion params */
	  GroupXML_t    sink_spin_insertion;         /*!< xml string holding sink spin insertion params */

	  PlusMinus     operation;
	};
	multi1d<Sinks_t> sink_ids;
      };

      multi1d<Props_t> prop_ids;
    } named_obj;

    std::string xml_file;  // Alternate XML file pattern
  };


  //! Inline measurement of Wilson loops
  /*! \ingroup inlinehadron */
  class InlineMesonSpec : public AbsInlineMeasurement 
  {
  public:
    ~InlineMesonSpec() {}
    InlineMesonSpec(const InlineMesonSpecParams& p) : params(p) {}
    InlineMesonSpec(const InlineMesonSpec& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out); 

  protected:
    //! Do the measurement
    void func(const unsigned long update_no,
	      XMLWriter& xml_out); 

  private:
    InlineMesonSpecParams params;
  };

};

#endif
