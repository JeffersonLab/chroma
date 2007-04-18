// -*- C++ -*-
// $Id: inline_stag_to_wils.h,v 3.2 2007-04-18 02:32:26 edwards Exp $
/*! \file
 * \brief Inline construction of propagator
 *
 * Propagator calculations
 */

#ifndef __inline_StagToWils_h__
#define __inline_StagToWils_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineStagToWilsEnv 
  {
    extern const std::string name;
    bool registerAll();
  }

  //! Parameter structure
  /*! \ingroup inlinehadron */ 
  struct InlineStagToWilsParams 
  {
    InlineStagToWilsParams();
    InlineStagToWilsParams(XMLReader& xml_in, const std::string& path);
    void writeXML(XMLWriter& xml_out, const std::string& path);

    unsigned long     frequency;

    struct NamedObject_t
    {
      std::string     stag_id;
      std::string     wils_id;
    } named_obj;

    std::string xml_file;  // Alternate XML file pattern
  };

  //! Inline measurement of staggered-to-wilson conversion
  /*! \ingroup inlinehadron */
  class InlineStagToWils : public AbsInlineMeasurement 
  {
  public:
    ~InlineStagToWils() {}
    InlineStagToWils(const InlineStagToWilsParams& p) : params(p) {}
    InlineStagToWils(const InlineStagToWils& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out); 

  protected:
    //! Do the measurement
    void func(const unsigned long update_no,
	      XMLWriter& xml_out); 

  private:
    InlineStagToWilsParams params;
  };

}

#endif
