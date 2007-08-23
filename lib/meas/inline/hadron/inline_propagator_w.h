// -*- C++ -*-
// $Id: inline_propagator_w.h,v 3.5 2007-08-23 19:02:45 edwards Exp $
/*! \file
 * \brief Inline construction of propagator
 *
 * Propagator calculations
 */

#ifndef __inline_propagator_h__
#define __inline_propagator_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlinePropagatorEnv 
  {
    extern const std::string name;
    bool registerAll();
  }

  //! Parameter structure
  /*! \ingroup inlinehadron */ 
  struct InlinePropagatorParams 
  {
    InlinePropagatorParams();
    InlinePropagatorParams(XMLReader& xml_in, const std::string& path);
    void writeXML(XMLWriter& xml_out, const std::string& path);

    unsigned long     frequency;

    ChromaProp_t      param;

    struct NamedObject_t
    {
      std::string     gauge_id;
      std::string     source_id;
      std::string     prop_id;
    } named_obj;

    std::string xml_file;  // Alternate XML file pattern
  };

  //! Inline propagator calculation
  /*! \ingroup inlinehadron */
  class InlinePropagator : public AbsInlineMeasurement 
  {
  public:
    ~InlinePropagator() {}
    InlinePropagator(const InlinePropagatorParams& p) : params(p) {}
    InlinePropagator(const InlinePropagator& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out); 

  protected:
    //! Do the measurement
    void func(const unsigned long update_no,
	      XMLWriter& xml_out); 

  private:
    InlinePropagatorParams params;
  };

}

#endif
