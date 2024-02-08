// -*- C++ -*-
// inline_ringed_ferm_w.h
/*! \file
 * \brief ringed fermion Inline task
 *
 *  ---More documentation here--- 
 */

#ifndef __inline_gluon_cc_h__
#define __inline_gluon_cc_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
//Do we need qprop_io.h ?
#include "io/qprop_io.h"  
#include "io/xml_group_reader.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineGluonCCEnv 
  {
    extern const std::string name;
    bool registerAll();
  }

    //! Parameter structure
    /*! \ingroup inlinehadron */ 
    struct InlineGluonCCParams 
    {
      InlineGluonCCParams();
      InlineGluonCCParams(XMLReader& xml_in, const std::string& path);
      void writeXML(XMLWriter& xml_out, const std::string& path);

      unsigned long     frequency;

      struct Param_t
      {
	int zmax;
	int zdir;
	std::string db;
	Real k5;
      } param;

      struct NamedObject_t
      {
	std::string     gauge_in;       /*!< Gauge fields */
      } named_obj;

      std::string xml_file;  // Alternate XML file pattern
    };

    /*! \ingroup inlinehadron */
    class InlineGluonCC : public AbsInlineMeasurement 
    {
    public:
      ~InlineGluonCC() {}
      InlineGluonCC(const InlineGluonCCParams& p) : params(p) {}
      InlineGluonCC(const InlineGluonCC& p) : params(p.params) {}

      unsigned long getFrequency(void) const {return params.frequency;}

      //! Do the measurement
      void operator()(const unsigned long update_no,
		      XMLWriter& xml_out); 

    protected:
      //! Do the measurement
      void func(const unsigned long update_no,
		XMLWriter& xml_out); 

    private:

      void field_str(multi2d<LatticeColorMatrix>& F, multi1d<LatticeColorMatrix>& u);

      InlineGluonCCParams params;
    };



} // namespace Chroma

#endif
