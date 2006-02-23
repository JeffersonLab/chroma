// -*- C++ -*-
// $Id: inline_coulgauge.h,v 1.1 2006-02-23 20:58:00 edwards Exp $
/*! \file
 *  \brief Inline coulomb (and landau) gauge fixing loops
 */

#ifndef __inline_gfix_h__
#define __inline_gfix_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"

namespace Chroma 
{ 
  /*! \ingroup inlinegfix */
  namespace InlineCoulGaugeEnv 
  {
    extern const std::string name;
    extern const bool registered;

    //! Parameter structure
    /*! \ingroup inlinegfix */
    struct Params 
    {
      Params();
      Params(XMLReader& xml_in, const std::string& path);
      void write(XMLWriter& xml_out, const std::string& path);

      unsigned long frequency;

      struct Param_t
      {
	Real GFAccu;      /*!< desired accuracy for gauge fixing */
	int  GFMax;       /*!< maximal number of gauge fixing iterations */
	bool OrDo;        /*!< use overrelaxation or not */
	Real OrPara;      /*!< overrelaxation parameter */
	int  j_decay;     /*!< direction perpendicular to slices to be gauge fixed */
      } param;

      struct NamedObject_t
      {
	std::string     gfix_id;
      } named_obj;

    };


    //! Inline measurement of Coulomb gauge fixing
    /*! \ingroup inlinegfix */
    class InlineMeas : public AbsInlineMeasurement 
    {
    public:
      ~InlineMeas() {}
      InlineMeas(const Params& p) : params(p) {}
      InlineMeas(const InlineMeas& p) : params(p.params) {}

      unsigned long getFrequency(void) const {return params.frequency;}

      //! Do the measurement
      void operator()(const multi1d<LatticeColorMatrix>& u,
		      XMLBufferWriter& gauge_xml,
		      const unsigned long update_no,
		      XMLWriter& xml_out); 

    private:
      Params params;
    };

  }

}

#endif
