// -*- C++ -*-
// $Id: inline_mres_w.h,v 1.2 2005-04-10 20:41:11 edwards Exp $
/*! \file
 * \brief Inline mres calculations
 *
 * Mres calculations
 */

#ifndef __inline_mres_h__
#define __inline_mres_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"

namespace Chroma 
{ 
  namespace InlineMresEnv 
  {
    extern const std::string name;
    extern const bool registered;
  }

  //! Parameter structure
  struct InlineMresParams 
  {
    InlineMresParams();
    InlineMresParams(XMLReader& xml_in, const std::string& path);
    void write(XMLWriter& xml_out, const std::string& path);

    unsigned long     frequency;

    struct Param_t
    {
      multi1d<int>    nrow;
    } param;

    std::string       stateInfo;

    struct Prop_t
    {
      string          prop_file;
    } prop;
  };


  //! Inline measurement of Wilson loops
  class InlineMres : public AbsInlineMeasurement 
  {
  public:
    ~InlineMres() {}
    InlineMres(const InlineMresParams& p) : params(p) {}
    InlineMres(const InlineMres& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const multi1d<LatticeColorMatrix>& u,
		    XMLBufferWriter& gauge_xml,
		    const unsigned long update_no,
		    XMLWriter& xml_out); 

  private:
    InlineMresParams params;
  };

};

#endif
