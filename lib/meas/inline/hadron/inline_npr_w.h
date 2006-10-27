// -*- C++ -*-
// $Id: inline_npr_w.h,v 1.2 2006-10-27 04:12:20 kostas Exp $
/*! \file
 * \brief Inline construction of Landau gauge propagator
 *
 * Landau gauge Propagators for NPR  calculations
 */

#ifndef __inline_npr_h__
#define __inline_npr_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineNprEnv 
  {
    extern const std::string name;
    bool registerAll();
  }

  //! Parameter structure
  /*! \ingroup inlinehadron */ 
  struct InlineNprParams 
  {
    InlineNprParams();
    InlineNprParams(XMLReader& xml_in, const std::string& path);
    void writeXML(XMLWriter& xml_out, const std::string& path);

    unsigned long     frequency;

    ChromaProp_t      param;
    std::string       stateInfo;

    struct NamedObject_t
    {
      std::string     gauge_id;
      //std::string     source_id;
      std::string     prop_id;
    } named_obj;

    std::string xml_file;  // Alternate XML file pattern
  };

  //! Inline measurement of Wilson loops
  /*! \ingroup inlinehadron */
  class InlineNpr : public AbsInlineMeasurement 
  {
  public:
    ~InlineNpr() {}
    InlineNpr(const InlineNprParams& p) : params(p) {}
    InlineNpr(const InlineNpr& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out); 

  protected:
    //! Do the measurement
    void func(const unsigned long update_no,
	      XMLWriter& xml_out); 
    
    void make_source(LatticePropagator& src,
		     const Handle<const ConnectState>& state,
		     const mult1d<ind>& t_source,
		     int mu) ;

  private:
    InlineNprParams params;
  };

}

#endif
