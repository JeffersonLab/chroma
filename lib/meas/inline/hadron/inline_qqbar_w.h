// -*- C++ -*-
// $Id: inline_qqbar_w.h,v 3.3 2007-04-18 02:32:26 edwards Exp $
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
    bool registerAll();
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
    } param;

    struct NamedObject_t
    {
      string           gauge_id;        /*!< Input gauge field */
      multi1d<string>  prop_ids;        /*!< Input sink smeared propagators */
      string           qqbar_file;      /*!< qqbar output file */
    } named_obj;
  };


  //! Inline measurement of quark-antiquark 2-pt correlators
  /*! \ingroup inlinehadron */
  class InlineQQbar : public AbsInlineMeasurement 
  {
  public:
    ~InlineQQbar() {}
    InlineQQbar(const InlineQQbarParams& p) : params(p) {}
    InlineQQbar(const InlineQQbar& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out); 

  private:
    InlineQQbarParams params;
  };

};

#endif
