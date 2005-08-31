// -*- C++ -*-
// $Id: inline_qqbar_w.h,v 1.1 2005-08-31 05:50:00 edwards Exp $
/*! \file
 * \brief Inline construction of qqbar
 *
 * QQbar calcs
 */

#ifndef __inline_qqbar_h__
#define __inline_qqbar_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineQQbarEnv 
  {
    extern const std::string name;
    extern const bool registered;
  }

  //! Parameter structure
  /*! \ingroup inlinehadron */
  struct InlineQQbarParams 
  {
    InlineQQbarParams();
    InlineQQbarParams(XMLReader& xml_in, const std::string& path);
    void write(XMLWriter& xml_out, const std::string& path);

    unsigned long      frequency;

    struct Param_t
    {
      bool             Dirac_basis;     /*!< Use the Dirac basis for output? */
      multi1d<int>     nrow;		/*!< Lattice dimension */
    } param;

    struct Prop_t
    {
      multi1d<string>  prop_file;    /*!< The file is expected to be in SciDAC format! */
      string           qqbar_file;   /*!< The file is expected to be in SciDAC format! */
    } prop;
  };


  //! Inline measurement of Wilson loops
  /*! \ingroup inlinehadron */
  class InlineQQbar : public AbsInlineMeasurement 
  {
  public:
    ~InlineQQbar() {}
    InlineQQbar(const InlineQQbarParams& p) : params(p) {}
    InlineQQbar(const InlineQQbar& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const multi1d<LatticeColorMatrix>& u,
		    XMLBufferWriter& gauge_xml,
		    const unsigned long update_no,
		    XMLWriter& xml_out); 

  private:
    InlineQQbarParams params;
  };

};

#endif
