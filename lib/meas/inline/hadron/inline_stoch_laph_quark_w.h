// -*- C++ -*-
// $Id: inline_stoch_laph_quark_w.h,v 3.4 2009-08-21 14:54:51 colin Exp $
/*! \file
 * \brief Compute the laph-diluted quark sources and sinks. Write them 
 *  out to db files.  Uses a QuarkSourceSinkHandler.
 *
 * Propagator calculation on laph diluted sources
 */

#ifndef __INLINE_STOCH_LAPH_QUARK_W_H__
#define __INLINE_STOCH_LAPH_QUARK_W_H__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "meas/laph/laph.h"


namespace Chroma { 
  namespace InlineStochLaphQuarkEnv {

 // **************************************************************


extern const std::string name;
bool registerAll();

    /*! \ingroup inlinehadron */

class StochLaphQuarkInlineMeas : public AbsInlineMeasurement 
{

   XMLReader xml_rdr;   // holds the XML input for this inline
                        // measurement, for use by the operator()
                        // member below

   struct SinkComputation {
      LaphEnv::LaphNoiseInfo Noise;
      int SourceTime;
      int FileIndex; 

    SinkComputation(const LaphEnv::LaphNoiseInfo& in_noise, int in_time, int in_file_ind)
       : Noise(in_noise), SourceTime(in_time), FileIndex(in_file_ind) {}
   };

   struct SourceComputation {
      LaphEnv::LaphNoiseInfo Noise;
      int FileIndex; 

    SourceComputation(const LaphEnv::LaphNoiseInfo& in_noise, int in_file_ind)
       : Noise(in_noise), FileIndex(in_file_ind) {}
   };

   list<SinkComputation> sinkComputations;
   list<SourceComputation> sourceComputations;

 public:

   StochLaphQuarkInlineMeas(XMLReader& xml_in, const std::string& path) 
                              : xml_rdr(xml_in, path) {}

   ~StochLaphQuarkInlineMeas() {}
      
   void setSinkComputations(int TimeExtent);
   void clearSinkComputations();

   void setSourceComputations();
   void clearSourceComputations();

      //! Do the measurement
   void operator()(const unsigned long update_no, XMLWriter& xml_out); 

   unsigned long getFrequency() const {return 0;}
   
};
	

// ***********************************************************
  }
}

#endif
