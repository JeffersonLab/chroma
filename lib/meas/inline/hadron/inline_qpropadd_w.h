// -*- C++ -*-
// $Id: inline_qpropadd_w.h,v 3.3 2007-04-18 02:32:26 edwards Exp $
/*! \file
 * \brief Inline measurement of qpropadd
 *
 * Form-factors
 */

#ifndef __inline_qpropadd_w_h__
#define __inline_qpropadd_w_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineQpropAddEnv 
  {
    extern const std::string name;
    bool registerAll();

    //! Parameter structure
    /*! \ingroup inlinehadron */
    struct Params 
    {
      Params();
      Params(XMLReader& xml_in, const std::string& path);
      void writeXML(XMLWriter& xml_out, const std::string& path);

      unsigned long      frequency;

      struct NamedObject_t
      {
	Real             factorA;
	std::string      propA;    
	Real             factorB;
	std::string      propB;   
	std::string      propApB; 
      } named_obj;
    };


    //! Inline measurement of to add two props
    /*! \ingroup inlinehadron */
    class InlineMeas : public AbsInlineMeasurement 
    {
    public:
      ~InlineMeas() {}
      InlineMeas(const Params& p) : params(p) {}
      InlineMeas(const InlineMeas& p) : params(p.params) {}

      unsigned long getFrequency(void) const {return params.frequency;}

      //! Do the measurement
      void operator()(const unsigned long update_no,
		      XMLWriter& xml_out); 

    private:
      Params params;
    };

  }

}

#endif
