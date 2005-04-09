// -*- C++ -*-
// $Id: inline_multipole_w.h,v 1.2 2005-04-09 23:15:42 edwards Exp $
/*! \file
 *  \brief Inline multipole measurements
 */

#ifndef __inline_multipole_h__
#define __inline_multipole_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"

namespace Chroma 
{ 
  namespace InlineMultipoleEnv 
  {
    extern const std::string name;
    extern const bool registered;
  }

  //! Parameter structure
  struct InlineMultipoleParams 
  {
    InlineMultipoleParams();
    InlineMultipoleParams(XMLReader& xml_in, const std::string& path);
    void write(XMLWriter& xml, const std::string& path);

    unsigned long frequency;

    struct Param_t
    {
      int max_L;
    } param;

    //! Propagators
    struct Prop_t
    {
      std::string       seqprop_file;       // sequential prop
      int               GammaInsertion;     // second gamma insertion
    };

    //! Multipole output
    struct Multipole_out_t
    {
      std::string       prop_file;          // input forward prop
      multi1d<Prop_t>   seqprops;
    } pole;


    std::string         multipole_file;     // where to write XML - maybe empty
  };

  //! Inline measurement of Wilson loops
  class InlineMultipole : public AbsInlineMeasurement 
  {
  public:
    ~InlineMultipole() {}
    InlineMultipole(const InlineMultipoleParams& p) : params(p) {}
    InlineMultipole(const InlineMultipole& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const multi1d<LatticeColorMatrix>& u,
		    XMLBufferWriter& gauge_xml,
		    const unsigned long update_no,
		    XMLWriter& xml_out); 

  private:
    InlineMultipoleParams params;
  };

};

#endif
