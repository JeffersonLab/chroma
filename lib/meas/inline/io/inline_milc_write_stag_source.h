// -*- C++ -*-
/*! \file
 * \brief Inline task to read a USQCD DD Pairs Prop 
 *
 *
 */

#ifndef __inline_milc_write_stag_source_h__
#define __inline_milc_write_stag_source_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/xml_group_reader.h"

namespace Chroma 
{ 
  /*! \ingroup inlineio */
  namespace InlineMILCWriteStagSourceEnv
  {
    extern const std::string name;
    bool registerAll();
  }

  //! Parameter structure
  /*! \ingroup inlineio */
  struct InlineMILCWriteStagSourceParams
  {
    InlineMILCWriteStagSourceParams();
    InlineMILCWriteStagSourceParams(XMLReader& xml_in, const std::string& path);
    void write(XMLWriter& xml_out, const std::string& path);

    unsigned long frequency;

    // We will recreate the source from the header
    std::string   source_id;
    int 		  t_slice;

    std::string   output_file_name;
    QDP_volfmt_t  qio_volfmt;
    QDP_serialparallel_t parallel_io;
    std::string   precision;
    std::string   xml_file;
  };

  //! Inline writing of memory objects
  /*! \ingroup inlineio */
  class InlineMILCWriteStagSource : public AbsInlineMeasurement
  {
  public:
    ~InlineMILCWriteStagSource() {}
    InlineMILCWriteStagSource(const InlineMILCWriteStagSourceParams& p) : params(p) {}
    InlineMILCWriteStagSource(const InlineMILCWriteStagSource& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the writing
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out); 

  private:
    void func(unsigned long update_no, XMLWriter& xml_out);
    InlineMILCWriteStagSourceParams params;
  };

};

#endif
