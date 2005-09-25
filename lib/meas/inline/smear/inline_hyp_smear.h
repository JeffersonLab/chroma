// -*- C++ -*-
// $Id: inline_hyp_smear.h,v 2.0 2005-09-25 21:04:39 edwards Exp $
/*! \file
 *  \brief Inline Hyp smearing
 */

#ifndef __inline_hyp_smear_h__
#define __inline_hyp_smear_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"

namespace Chroma 
{ 
  /*! \ingroup inlinesmear */
  namespace InlineHypSmearEnv 
  {
    extern const std::string name;
    extern const bool registered;
  }

  //! Parameter structure
  /*! \ingroup inlinesmear */
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
    } param;

    struct NamedObject_t
    {
      std::string   hyp_id;             // Memory object for hyp config
    } named_obj;

  };


  //! Inline measurement of Wilson loops
  /*! \ingroup inlinesmear */
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
