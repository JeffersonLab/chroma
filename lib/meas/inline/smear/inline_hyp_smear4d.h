// -*- C++ -*-
// $Id: inline_hyp_smear4d.h,v 1.2 2005-08-19 05:46:48 edwards Exp $
/*! \file
 *  \brief Inline Hyp smearing
 */

#ifndef __inline_hyp_smear4d_h__
#define __inline_hyp_smear4d_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/enum_io/enum_cfgtype_io.h"
#include "io/enum_io/enum_qdpvolfmt_io.h"

namespace Chroma 
{ 
  /*! \ingroup inlinesmear */
  namespace InlineHypSmear4dEnv 
  {
    extern const std::string name;
    extern const bool registered;
  }

  //! Parameter structure
  /*! \ingroup inlinesmear */
  struct InlineHypSmear4dParams 
  {
    InlineHypSmear4dParams();
    InlineHypSmear4dParams(XMLReader& xml_in, const std::string& path);
    void write(XMLWriter& xml_out, const std::string& path);

    unsigned long frequency;

    struct Param_t
    {
      Real alpha1;			// Smearing parameters
      Real alpha2;
      Real alpha3;

      int num_smear;                    // Number of smearing iterations

      multi1d<int> nrow;		// Lattice dimension
    } param;

    struct Hyp_t
    {
      QDP_volfmt_t  hyp_volfmt;	   // single or multi file volume format
      string   hyp_file;	   // storage for hyp config
    } hyp;

  };


  //! Inline measurement of Wilson loops
  /*! \ingroup inlinesmear */
  class InlineHypSmear4d : public AbsInlineMeasurement 
  {
  public:
    ~InlineHypSmear4d() {}
    InlineHypSmear4d(const InlineHypSmear4dParams& p) : params(p) {}
    InlineHypSmear4d(const InlineHypSmear4d& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const multi1d<LatticeColorMatrix>& u,
		    XMLBufferWriter& gauge_xml,
		    const unsigned long update_no,
		    XMLWriter& xml_out); 

  private:
    InlineHypSmear4dParams params;
  };

};

#endif
