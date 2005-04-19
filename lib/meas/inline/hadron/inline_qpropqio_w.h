// -*- C++ -*-
// $Id: inline_qpropqio_w.h,v 1.2 2005-04-19 20:05:22 edwards Exp $
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
    extern const bool registered;
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
      multi1d<int>     nrow;
    } param;

    struct Prop_t
    {
      string           prop_in_file;    // The files is expected to be in SciDAC format!
      string           prop_out_file;   // The files is expected to be in SciDAC format!
      QDP_volfmt_t     prop_out_volfmt; // volume format (SINGLEFILE or MULTIFILE)
    } prop;

  };


  //! Inline measurement of Wilson loops
  /*! \ingroup inlinehadron */
  class InlineQpropQIO : public AbsInlineMeasurement 
  {
  public:
    ~InlineQpropQIO() {}
    InlineQpropQIO(const InlineQpropQIOParams& p) : params(p) {}
    InlineQpropQIO(const InlineQpropQIO& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const multi1d<LatticeColorMatrix>& u,
		    XMLBufferWriter& gauge_xml,
		    const unsigned long update_no,
		    XMLWriter& xml_out); 

  private:
    InlineQpropQIOParams params;
  };

};

#endif
