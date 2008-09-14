// -*- C++ -*-
// $Id: inline_rng.h,v 3.1 2008-09-14 02:26:11 edwards Exp $
/*! \file
 * \brief Inline task to read and write RNG seed
 */

#ifndef __inline_rng_h__
#define __inline_rng_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"

namespace Chroma 
{ 
  /*! \ingroup inlineio */
  namespace InlineSetRNGEnv 
  {
    bool registerAll();

    //! Parameter structure
    /*! \ingroup inlineio */
    struct Params 
    {
      Params();
      Params(XMLReader& xml_in, const std::string& path);
      void writeXML(XMLWriter& xml_out, const std::string& path);

      unsigned long frequency;

      QDP::Seed       ran_seed;
    };

    //! Set the RNG seed
    /*! \ingroup inlineio */
    class InlineMeas : public AbsInlineMeasurement 
    {
    public:
      ~InlineMeas() {}
      InlineMeas(const Params& p) : params(p) {}
      InlineMeas(const InlineMeas& p) : params(p.params) {}

      unsigned long getFrequency(void) const {return params.frequency;}

      //! Do the writing
      void operator()(const unsigned long update_no,
		      XMLWriter& xml_out); 

    private:
      Params params;
    };

  }
}

#endif
