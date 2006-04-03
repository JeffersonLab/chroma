// -*- C++ -*-
// $Id: inline_mres_w.h,v 3.0 2006-04-03 04:59:02 edwards Exp $
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
  /*! \ingroup inlinehadron */
  namespace InlineMresEnv 
  {
    extern const std::string name;
    extern const bool registered;
  }

  //! Parameter structure
  /*! \ingroup inlinehadron */
  struct InlineMresParams 
  {
    InlineMresParams();
    InlineMresParams(XMLReader& xml_in, const std::string& path);
    void write(XMLWriter& xml_out, const std::string& path);

    unsigned long     frequency;

    struct Param_t
    {
      multi1d<int>    nrow;
      std::string     fermact;
    } param;

    std::string       stateInfo;

    struct NamedObject_t
    {
      std::string     gauge_id;
      std::string     prop_id;
    } named_obj;

    std::string xml_file;  // Alternate XML file pattern
  };


  //! Inline measurement of Wilson loops
  /*! \ingroup inlinehadron */
  class InlineMres : public AbsInlineMeasurement 
  {
  public:
    ~InlineMres() {}
    InlineMres(const InlineMresParams& p) : params(p) {}
    InlineMres(const InlineMres& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out); 

  protected:
    //! Do the measurement
    void func(const unsigned long update_no,
	      XMLWriter& xml_out); 

  private:
    InlineMresParams params;
  };

};

#endif
