// -*- C++ -*-
// $Id: inline_usqcd_read_ddpairs_prop.h,v 3.4 2008-05-02 21:03:35 bjoo Exp $
/*! \file
 * \brief Inline task to read a USQCD DD Pairs Prop 
 *
 *
 */

#ifndef __inline_usqcd_read_ddpairs_prop_h__
#define __inline_usqcd_read_ddpairs_prop_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/xml_group_reader.h"

namespace Chroma 
{ 
  /*! \ingroup inlineio */
  namespace InlineUSQCDReadDDPairsPropEnv
  {
    extern const std::string name;
    bool registerAll();
  }

  //! Parameter structure
  /*! \ingroup inlineio */
  struct InlineUSQCDReadDDPairsPropParams
  {
    InlineUSQCDReadDDPairsPropParams();
    InlineUSQCDReadDDPairsPropParams(XMLReader& xml_in, const std::string& path);
    void write(XMLWriter& xml_out, const std::string& path);

    unsigned long frequency;

    // The DDPairs Format contains both the source and the propagator
    // We will create a named object for each for now. One can always
    // Call a delete for the source later
    
    std::string   source_id;
    std::string   prop_id;
    
    std::string   input_file_name;
    std::string   prop_xml;    
    std::string   xml_file;
  };

  //! Inline writing of memory objects
  /*! \ingroup inlineio */
  class InlineUSQCDReadDDPairsProp : public AbsInlineMeasurement 
  {
  public:
    ~InlineUSQCDReadDDPairsProp() {}
    InlineUSQCDReadDDPairsProp(const InlineUSQCDReadDDPairsPropParams& p) : params(p) {}
    InlineUSQCDReadDDPairsProp(const InlineUSQCDReadDDPairsProp& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the writing
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out); 

  private:
    void func(unsigned long update_no, XMLWriter& xml_out);
    InlineUSQCDReadDDPairsPropParams params;
  };

};

#endif
