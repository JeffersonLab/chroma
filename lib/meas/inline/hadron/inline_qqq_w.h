// -*- C++ -*-
// $Id: inline_qqq_w.h,v 1.2 2005-04-19 20:05:22 edwards Exp $
/*! \file
 * \brief Inline construction of qqq_w
 *
 * QQQ calcs
 */

#ifndef __inline_qqq_h__
#define __inline_qqq_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineQQQEnv 
  {
    extern const std::string name;
    extern const bool registered;
  }

  //! Parameter structure
  /*! \ingroup inlinehadron */
  struct InlineQQQParams 
  {
    InlineQQQParams();
    InlineQQQParams(XMLReader& xml_in, const std::string& path);
    void write(XMLWriter& xml_out, const std::string& path);

    unsigned long      frequency;

    struct Param_t
    {
      bool             Dirac_basis;     // Use the Dirac basis for output?
      multi1d<int>     nrow;		// Lattice dimension
    } param;

    struct Prop_t
    {
      multi1d<string>  prop_file;  // The file is expected to be in SciDAC format!
      string           qqq_file;   // The file is expected to be in SciDAC format!
    } prop;
  };


  //! Inline measurement of Wilson loops
  /*! \ingroup inlinehadron */
  class InlineQQQ : public AbsInlineMeasurement 
  {
  public:
    ~InlineQQQ() {}
    InlineQQQ(const InlineQQQParams& p) : params(p) {}
    InlineQQQ(const InlineQQQ& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const multi1d<LatticeColorMatrix>& u,
		    XMLBufferWriter& gauge_xml,
		    const unsigned long update_no,
		    XMLWriter& xml_out); 

  private:
    InlineQQQParams params;
  };

};

#endif
