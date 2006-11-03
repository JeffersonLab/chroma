// -*- C++ -*-
// $Id: inline_npr_w.h,v 1.5 2006-11-03 20:21:07 hwlin Exp $
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
#include "handle.h"
#include "state.h"

using namespace QDP;

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

    multi1d<int> NprSources ; // mu>3 or mu<0 means point source
    int max_mom2 ; // max p^2

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
		     const multi1d<int>& t_source,
		     int mu) ;

  private:
    InlineNprParams params;
  };

}

#endif
