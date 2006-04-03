// -*- C++ -*-
// $Id: inline_stout_smear.h,v 3.0 2006-04-03 04:59:04 edwards Exp $
/*! \file
 *  \brief Inline Stout smearing
 */

#ifndef __inline_stout_smear_h__
#define __inline_stout_smear_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"

namespace Chroma 
{ 
  /*! \ingroup inlinesmear */
  namespace InlineStoutSmearEnv 
  {
    extern const std::string name;
    extern const bool registered;
  }

  //! Parameter structure
  /*! \ingroup inlinesmear */
  struct InlineStoutSmearParams 
  {
    InlineStoutSmearParams();
    InlineStoutSmearParams(XMLReader& xml_in, const std::string& path);
    void write(XMLWriter& xml_out, const std::string& path);

    unsigned long frequency;

    struct Param_t
    {
      int link_smear_num;
      Real link_smear_fact;		// Smearing parameters

      int j_decay;			// Decay direction
      multi1d<int> nrow;		// Lattice dimension
    } param;

    struct NamedObject_t
    {
      std::string 	gauge_id;	/*!< Input gauge field */
      std::string 	stout_id;       /*!< Output memory object stout config */
    } named_obj;

  };

  //! Inline measurement of Wilson loops
  /*! \ingroup inlinesmear */
  class InlineStoutSmear : public AbsInlineMeasurement 
  {
  public:
    ~InlineStoutSmear() {}
    InlineStoutSmear(const InlineStoutSmearParams& p) : params(p) {}
    InlineStoutSmear(const InlineStoutSmear& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out); 

  private:
    InlineStoutSmearParams params;
  };

};

#endif
