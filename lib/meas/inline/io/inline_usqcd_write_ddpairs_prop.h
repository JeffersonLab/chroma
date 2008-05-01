// -*- C++ -*-
// $Id: inline_usqcd_write_ddpairs_prop.h,v 1.2 2008-05-01 19:32:56 bjoo Exp $
/*! \file
 * \brief Inline task to read a USQCD DD Pairs Prop 
 *
 *
 */

#ifndef __inline_usqcd_write_ddpairs_prop_h__
#define __inline_usqcd_write_ddpairs_prop_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/xml_group_reader.h"

namespace Chroma 
{ 
  /*! \ingroup inlineio */
  namespace InlineUSQCDWriteDDPairsPropEnv
  {
    extern const std::string name;
    bool registerAll();
  }

  //! Parameter structure
  /*! \ingroup inlineio */
  struct InlineUSQCDWriteDDPairsPropParams
  {
    InlineUSQCDWriteDDPairsPropParams();
    InlineUSQCDWriteDDPairsPropParams(XMLReader& xml_in, const std::string& path);
    void write(XMLWriter& xml_out, const std::string& path);

    unsigned long frequency;

    // We will recreate the source from the header
    std::string   prop_id;
    std::string   gauge_id;

    std::string   output_file_name;
    QDP_volfmt_t  qio_volfmt;
    std::string   precision;

    std::string   xml_file;
  };

  //! Inline writing of memory objects
  /*! \ingroup inlineio */
  class InlineUSQCDWriteDDPairsProp : public AbsInlineMeasurement 
  {
  public:
    ~InlineUSQCDWriteDDPairsProp() {}
    InlineUSQCDWriteDDPairsProp(const InlineUSQCDWriteDDPairsPropParams& p) : params(p) {}
    InlineUSQCDWriteDDPairsProp(const InlineUSQCDWriteDDPairsProp& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the writing
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out); 

  private:
    void func(unsigned long update_no, XMLWriter& xml_out);
    InlineUSQCDWriteDDPairsPropParams params;
  };

};

#endif
