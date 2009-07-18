// -*- C++ -*-
// $Id: inline_stoch_laph_quark_w.h,v 3.3 2009-07-18 02:34:47 jbulava Exp $
/*! \file
 * \brief Compute the laph-diluted sources and solutions. Write them out to a single db file.  
 *
 * Propagator calculation on a laph diluted source 
 */

#ifndef __inline_stoch_laph_quark_w_h__
#define __inline_stoch_laph_quark_w_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"
#include "meas/laph/laph.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineStochLaphQuarkEnv 
  {
    bool registerAll();

    //! Inline task for compute LatticeColorVector matrix elements of a propagator
    /*! \ingroup inlinehadron */
    class StochLaphQuarkInlineMeas : public AbsInlineMeasurement 
    {
    public:
      ~StochLaphQuarkInlineMeas() {}
      
			StochLaphQuarkInlineMeas(XMLReader& xml_in, const std::string& path) 
			: xml_rdr(xml_in, path) {}

      //! Do the measurement
      void operator()(const unsigned long update_no,
		                  XMLWriter& xml_out); 

			long unsigned int getFrequency() const {return 0;}
   
		private:
			XMLReader xml_rdr;
		};
	

  } // namespace PropMatElemColorVec

}

#endif
