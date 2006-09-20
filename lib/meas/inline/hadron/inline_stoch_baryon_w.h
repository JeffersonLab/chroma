// -*- C++ -*-
// $Id: inline_stoch_baryon_w.h,v 3.4 2006-09-20 20:28:02 edwards Exp $
/*! \file
 * \brief Inline measurement of stochastic baryon operator
 *
 * Form-factors
 */

#ifndef __inline_stoch_baryon_h__
#define __inline_stoch_baryon_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineStochBaryonEnv 
  {
    extern const std::string name;
    bool registerAll();
  }

  //! Parameter structure
  /*! \ingroup inlinehadron */
  struct InlineStochBaryonParams 
  {
    InlineStochBaryonParams();
    InlineStochBaryonParams(XMLReader& xml_in, const std::string& path);
    void write(XMLWriter& xml_out, const std::string& path);

    unsigned long      frequency;

    struct Param_t
    {
      int              mom2_max;           /*!< (mom)^2 <= mom2_max */

      std::string      baryon_operator;        /*!< baryon operator xml */
      std::string      baryon_operator_type;   /*!< baryon operator name */

    } param;

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


  //! Inline measurement of stochastic baryon operators
  /*! \ingroup inlinehadron */
  class InlineStochBaryon : public AbsInlineMeasurement 
  {
  public:
    ~InlineStochBaryon() {}
    InlineStochBaryon(const InlineStochBaryonParams& p) : params(p) {}
    InlineStochBaryon(const InlineStochBaryon& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out); 

  protected:
    //! Do the measurement
    void func(const unsigned long update_no,
	      XMLWriter& xml_out); 

  private:
    InlineStochBaryonParams params;
  };

}

#endif
