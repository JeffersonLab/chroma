// -*- C++ -*-
// $Id: inline_ritz_H_w.h,v 3.0 2006-04-03 04:59:01 edwards Exp $
/*! \file
 * \brief Inline construction of eigenvalues (Ritz)
 *
 * Eigenvalue calculations
 */

#ifndef __inline_ritz_H_w_h__
#define __inline_ritz_H_w_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/eigen_io.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineRitzEnv 
  {
    extern const std::string name;
    extern const bool registered;
  }

  //! Parameter structure
  /*! \ingroup inlinehadron */ 
  struct InlineRitzParams 
  {
    InlineRitzParams();
    InlineRitzParams(XMLReader& xml_in, const std::string& path);
    void write(XMLWriter& xml_out, const std::string& path);

    unsigned long     frequency;
    int            version;
    std::string    fermact;
    RitzParams_t   ritz_params;
    std::string       stateInfo;
    
    struct NamedObject_t
    {
      std::string     gauge_id;
      std::string     eigen_id;
    } named_obj;

    std::string xml_file;  // Alternate XML file pattern
  };

  //! Inline measurement of Wilson loops
  /*! \ingroup inlinehadron */
  class InlineRitz : public AbsInlineMeasurement 
  {
  public:
    ~InlineRitz() {}
    InlineRitz(const InlineRitzParams& p) : params(p) {}
    InlineRitz(const InlineRitz& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out); 

  protected:
    //! Do the measurement
    void func(const unsigned long update_no,
	      XMLWriter& xml_out); 

  private:
    InlineRitzParams params;
  };

};

#endif
