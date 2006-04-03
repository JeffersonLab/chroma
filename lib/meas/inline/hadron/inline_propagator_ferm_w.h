// -*- C++ -*-
// $Id: inline_propagator_ferm_w.h,v 3.0 2006-04-03 04:59:02 edwards Exp $
/*! \file
 * \brief Inline construction of propagator returning only a single lattice fermion
 *
 * Propagator calculations
 */

#ifndef __inline_propagator_ferm_h__
#define __inline_propagator_ferm_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlinePropagatorFermEnv 
  {
    extern const std::string name;
    extern const bool registered;
  }

  //! Parameter structure
  /*! \ingroup inlinehadron */ 
  struct InlinePropagatorFermParams 
  {
    InlinePropagatorFermParams();
    InlinePropagatorFermParams(XMLReader& xml_in, const std::string& path);

    unsigned long     frequency;

    ChromaProp_t      param;
    std::string       stateInfo;

    struct NamedObject_t
    {
      std::string     gauge_id;
      std::string     source_id;
      std::string     prop_id;
    } named_obj;

    std::string xml_file;  // Alternate XML file pattern
  };

  //! Reader
  /*! \ingroup inlinehadron */
  void read(XMLReader& xml, const string& path, InlinePropagatorFermParams& param);

  //! Writer
  /*! \ingroup inlinehadron */
  void write(XMLWriter& xml, const string& path, const InlinePropagatorFermParams& param);


  //! Inline measurement of Wilson loops
  /*! \ingroup inlinehadron */
  class InlinePropagatorFerm : public AbsInlineMeasurement 
  {
  public:
    ~InlinePropagatorFerm() {}
    InlinePropagatorFerm(const InlinePropagatorFermParams& p) : params(p) {}
    InlinePropagatorFerm(const InlinePropagatorFerm& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out); 

  protected:
    //! Do the measurement
    void func(const unsigned long update_no,
	      XMLWriter& xml_out); 

  private:
    InlinePropagatorFermParams params;
  };

};

#endif
