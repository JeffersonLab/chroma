// -*- C++ -*-
// $Id: inline_qpropqio_w.h,v 3.3 2007-04-18 02:32:26 edwards Exp $
/*! \file
 * \brief Inline measurement of qpropqio
 *
 * Form-factors
 */

#ifndef __inline_qpropqio_h__
#define __inline_qpropqio_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineQpropQIOEnv 
  {
    extern const std::string name;
    bool registerAll();
  }

  //! Parameter structure
  /*! \ingroup inlinehadron */
  struct InlineQpropQIOParams 
  {
    InlineQpropQIOParams();
    InlineQpropQIOParams(XMLReader& xml_in, const std::string& path);
    void write(XMLWriter& xml_out, const std::string& path);

    unsigned long      frequency;

    struct Param_t
    {
    } param;

    struct Prop_t
    {
      string           prop_in_file;    // The files is expected to be in SciDAC format!
      string           prop_out_file;   // The files is expected to be in SciDAC format!
      QDP_volfmt_t     prop_out_volfmt; // volume format (SINGLEFILE or MULTIFILE)
    } prop;

  };


  //! Inline task for quark prop IO
  /*! \ingroup inlinehadron */
  class InlineQpropQIO : public AbsInlineMeasurement 
  {
  public:
    ~InlineQpropQIO() {}
    InlineQpropQIO(const InlineQpropQIOParams& p) : params(p) {}
    InlineQpropQIO(const InlineQpropQIO& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out); 

  private:
    InlineQpropQIOParams params;
  };

};

#endif
