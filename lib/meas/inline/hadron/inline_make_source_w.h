// -*- C++ -*-
// $Id: inline_make_source_w.h,v 3.2 2007-04-18 02:32:26 edwards Exp $
/*! \file
 * \brief Inline construction of make_source
 *
 * Construct source for propagator calculations
 */

#ifndef __inline_make_source_h__
#define __inline_make_source_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineMakeSourceEnv 
  {
    extern const std::string name;
    bool registerAll();
  }

  //! Parameter structure
  /*! \ingroup inlinehadron */
  struct InlineMakeSourceParams 
  {
    InlineMakeSourceParams();
    InlineMakeSourceParams(XMLReader& xml_in, const std::string& path);
    void write(XMLWriter& xml_out, const std::string& path);

    unsigned long frequency;

    PropSourceConst_t  param;

    struct NamedObject_t
    {
      std::string     gauge_id;
      std::string     source_id;
    } named_obj;

    std::string xml_file;  // Alternate XML file pattern
  };

  //! Inline task creating sources for quark inversion
  /*! \ingroup inlinehadron */
  class InlineMakeSource : public AbsInlineMeasurement 
  {
  public:
    ~InlineMakeSource() {}
    InlineMakeSource(const InlineMakeSourceParams& p) : params(p) {}
    InlineMakeSource(const InlineMakeSource& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out); 

  protected:
    //! Do the measurement
    void func(const unsigned long update_no,
	      XMLWriter& xml_out); 

  private:
    InlineMakeSourceParams params;
  };

};

#endif
