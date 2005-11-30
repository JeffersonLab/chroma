// -*- C++ -*-
// $Id: inline_make_source_ferm_w.h,v 2.2 2005-11-30 04:46:39 edwards Exp $
/*! \file
 * \brief Inline construction of make_source for lattice fermions
 *
 * Construct source for propagator calculations
 */

#ifndef __inline_make_source_ferm_h__
#define __inline_make_source_ferm_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineMakeSourceFermEnv 
  {
    extern const std::string name;
    extern const bool registered;
  }

  //! Parameter structure
  /*! \ingroup inlinehadron */
  struct InlineMakeSourceFermParams 
  {
    InlineMakeSourceFermParams();
    InlineMakeSourceFermParams(XMLReader& xml_in, const std::string& path);

    unsigned long frequency;

    PropSourceConst_t  param;

    struct NamedObject_t
    {
      std::string     source_id;
    } named_obj;

    std::string xml_file;  // Alternate XML file pattern
  };


  //! Make source input
  void read(XMLReader& xml, const string& path, InlineMakeSourceFermParams& input);

  //! Make source output
  void write(XMLWriter& xml, const string& path, const InlineMakeSourceFermParams& input);


  //! Inline measurement of Wilson loops
  /*! \ingroup inlinehadron */
  class InlineMakeSourceFerm : public AbsInlineMeasurement 
  {
  public:
    ~InlineMakeSourceFerm() {}
    InlineMakeSourceFerm(const InlineMakeSourceFermParams& p) : params(p) {}
    InlineMakeSourceFerm(const InlineMakeSourceFerm& p) : params(p.params) {}

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
    InlineMakeSourceFermParams params;
  };

};

#endif
