// -*- C++ -*-
// $Id: inline_stoch_meson_w.h,v 3.2 2006-09-20 20:28:03 edwards Exp $
/*! \file
 * \brief Inline measurement of stochastic meson operator
 *
 * Form-factors
 */

#ifndef __inline_stoch_meson_h__
#define __inline_stoch_meson_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineStochMesonEnv 
  {
    extern const std::string name;
    bool registerAll();
  }

  //! Parameter structure
  /*! \ingroup inlinehadron */
  struct InlineStochMesonParams 
  {
    InlineStochMesonParams();
    InlineStochMesonParams(XMLReader& xml_in, const std::string& path);
    void write(XMLWriter& xml_out, const std::string& path);

    unsigned long      frequency;

    struct Param_t
    {
      int              mom2_max;           /*!< (mom)^2 <= mom2_max */
    } param;

    PropSourceSmear_t  source_smearing;
    PropSinkSmear_t    sink_smearing;

    struct Prop_t
    {
      //! Operators
      struct Operator_t
      {
	multi1d<std::string> soln_files;
      };

      std::string          op_file;
      multi1d<Operator_t>  op;
    };

    struct NamedObject_t
    {
      Prop_t                prop;
      std::string           gauge_id;
    } named_obj;

    std::string xml_file;  // Alternate XML file pattern
  };


  //! Inline measurement of stochastic meson operators
  /*! \ingroup inlinehadron */
  class InlineStochMeson : public AbsInlineMeasurement 
  {
  public:
    ~InlineStochMeson() {}
    InlineStochMeson(const InlineStochMesonParams& p) : params(p) {}
    InlineStochMeson(const InlineStochMeson& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out); 

  protected:
    //! Do the measurement
    void func(const unsigned long update_no,
	      XMLWriter& xml_out); 

  private:
    InlineStochMesonParams params;
  };

}

#endif
