// -*- C++ -*-
// $Id: inline_stoch_condensates_w.h,v 3.1 2007-08-10 21:27:02 edwards Exp $
/*! \file
 * \brief Inline measurement of (stochastic) condensates
 *
 */

#ifndef __inline_stoch_condensates_h__
#define __inline_stoch_condensates_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineStochCondensatesEnv 
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
	int              mom2_max;           /*!< (mom)^2 <= mom2_max */
      } param;

      struct NamedObject_t
      {
	std::string          gauge_id;
	multi1d<std::string> soln_files;
      } named_obj;

      std::string xml_file;  // Alternate XML file pattern
    };


    //! Inline measurement of stochastic condensates
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

    protected:
      //! Do the measurement
      void func(const unsigned long update_no,
		XMLWriter& xml_out); 

    private:
      Params params;
    };

  } // namespace InlineStochCondensatesEnv

} // namespace Chroma

#endif
