// -*- C++ -*-
// $Id: inline_diquark_w.h,v 1.1 2007-02-22 06:58:55 edwards Exp $
/*! \file
 * \brief Inline construction of the diquark within a QQQ
 *
 * Diquarks for QQQ calcs
 */

#ifndef __inline_diquark_w_h__
#define __inline_diquark_w_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineDiquarkEnv 
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

      struct Param_t
      {
	bool             Dirac_basis;     /*!< Use the Dirac basis for output? */
      } param;

      struct NamedObject_t
      {
	string           gauge_id;        /*!< Input gauge field */
	multi1d<string>  prop_ids;        /*!< Input sink smeared propagators */
	string           diquark_id;      /*!< Output qqq file */
      } named_obj;
    };


    //! Inline computation of diquarks for QQQ
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
