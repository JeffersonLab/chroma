// -*- C++ -*-
// $Id: inline_hyp_smear.h,v 1.1 2005-04-07 03:23:20 edwards Exp $
/*! \file
 *  \brief Inline Hyp smearing
 */

#ifndef __inline_hyp_smear_h__
#define __inline_hyp_smear_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/enum_io/enum_cfgtype_io.h"

namespace Chroma 
{ 
  namespace InlineHypSmearEnv 
  {
    extern const std::string name;
    extern const bool registered;
  }

  //! Parameter structure
  struct InlineHypSmearParams 
  {
    InlineHypSmearParams();
    InlineHypSmearParams(XMLReader& xml_in, const std::string& path);
    void write(XMLWriter& xml_out, const std::string& path);

    unsigned long frequency;

    struct Param_t
    {
      Real alpha1;			// Smearing parameters
      Real alpha2;
      Real alpha3;

      int num_smear;                    // Number of smearing iterations

      multi1d<int> nrow;		// Lattice dimension
      int j_decay;			// Direction corresponding to time

      /*
       *  Now some various rules for truncating the configuration
       */
      int trunc;			// Whether to truncate the output
      int t_start;			// Starting time slice
      int t_end;			// Ending time slice
    } param;

    struct Hyp_t
    {
      CfgType  cfg_type;       // storage order for stored gauge configuration
      string   hyp_file;
    } hyp;

  };


  //! Inline measurement of Wilson loops
  class InlineHypSmear : public AbsInlineMeasurement 
  {
  public:
    ~InlineHypSmear() {}
    InlineHypSmear(const InlineHypSmearParams& p) : params(p) {}
    InlineHypSmear(const InlineHypSmear& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const multi1d<LatticeColorMatrix>& u,
		    XMLBufferWriter& gauge_xml,
		    const unsigned long update_no,
		    XMLWriter& xml_out); 

  private:
    InlineHypSmearParams params;
  };

};

#endif
