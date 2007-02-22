// -*- C++ -*-
// $Id: inline_qqq_w.h,v 3.6 2007-02-22 04:52:42 edwards Exp $
/*! \file
 * \brief Inline construction of qqq_w
 *
 * QQQ calcs
 */

#ifndef __inline_qqq_w_h__
#define __inline_qqq_w_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineQQQEnv 
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
	multi1d<QQQSpinIndices_t> spin_indices;
	bool                      sparseP;
	bool                      Dirac_basis;     // Use the Dirac basis for output?
      } param;

      struct NamedObject_t
      {
	string           gauge_id;        /*!< Input gauge field */
	multi1d<string>  prop_ids;        /*!< Input sink smeared propagators */
	string           qqq_file;        /*!< Output qqq file */
      } named_obj;
    };


    //! Inline measurement of QQQ's
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
