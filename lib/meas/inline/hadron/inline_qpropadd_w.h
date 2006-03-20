// -*- C++ -*-
// $Id: inline_qpropadd_w.h,v 2.2 2006-03-20 04:22:03 edwards Exp $
/*! \file
 * \brief Inline measurement of qpropadd
 *
 * Form-factors
 */

#ifndef __inline_qpropadd_h__
#define __inline_qpropadd_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineQpropAddEnv 
  {
    extern const std::string name;
    extern const bool registered;
  }

  //! Parameter structure
  /*! \ingroup inlinehadron */
  struct InlineQpropAddParams 
  {
    InlineQpropAddParams();
    InlineQpropAddParams(XMLReader& xml_in, const std::string& path);
    void write(XMLWriter& xml_out, const std::string& path);

    unsigned long      frequency;

    struct NamedObject_t
    {
      string           propA;    
      string           propB;   
      string           propApB; 
    } named_obj;

  };


  //! Inline measurement of Wilson loops
  /*! \ingroup inlinehadron */
  class InlineQpropAdd : public AbsInlineMeasurement 
  {
  public:
    ~InlineQpropAdd() {}
    InlineQpropAdd(const InlineQpropAddParams& p) : params(p) {}
    InlineQpropAdd(const InlineQpropAdd& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out); 

  private:
    InlineQpropAddParams params;
  };

};

#endif
