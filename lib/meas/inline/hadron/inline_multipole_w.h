// -*- C++ -*-
// $Id: inline_multipole_w.h,v 1.3 2005-04-19 17:11:07 edwards Exp $
/*! \file
 *  \brief Inline multipole measurements
 */

#ifndef __inline_multipole_h__
#define __inline_multipole_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineMultipoleEnv 
  {
    extern const std::string name;
    extern const bool registered;
  }

  //! Parameter structure
  /*! \ingroup inlinehadron */
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

    std::string xml_file;  // Alternate XML file pattern
  };

  //! Inline measurement of Wilson loops
  /*! \ingroup inlinehadron */
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

  protected:
    //! Do the measurement
    void func(const multi1d<LatticeColorMatrix>& u,
	      XMLBufferWriter& gauge_xml,
	      const unsigned long update_no,
	      XMLWriter& xml_out); 

  private:
    InlineMultipoleParams params;
  };

};

#endif
