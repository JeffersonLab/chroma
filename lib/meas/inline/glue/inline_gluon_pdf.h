// -*- C++ -*-
// inline_ringed_ferm_w.h
/*! \file
 * \brief ringed fermion Inline task
 *
 *  ---More documentation here--- 
 */

#ifndef __inline_gluon_pdf_h__
#define __inline_gluon_pdf_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
//Do we need qprop_io.h ?
#include "io/qprop_io.h"  
#include "io/xml_group_reader.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineGluonPDFEnv 
  {
    extern const std::string name;
    bool registerAll();
  }

    //! Parameter structure
    /*! \ingroup inlinehadron */ 
    struct InlineGluonPDFParams 
    {
      InlineGluonPDFParams();
      InlineGluonPDFParams(XMLReader& xml_in, const std::string& path);
      void writeXML(XMLWriter& xml_out, const std::string& path);

      unsigned long     frequency;

      struct Param_t
      {
	int zmax;
	int zdir;
	std::string db;
	std::string smear;
	multi1d<int> mom;
	Real k5;
      } param;

      struct NamedObject_t
      {
	std::string     gauge_in;       /*!< Gauge fields */
      } named_obj;

      std::string xml_file;  // Alternate XML file pattern
    };

    /*! \ingroup inlinehadron */
    class InlineGluonPDF : public AbsInlineMeasurement 
    {
    public:
      ~InlineGluonPDF() {}
      InlineGluonPDF(const InlineGluonPDFParams& p) : params(p) {}
      InlineGluonPDF(const InlineGluonPDF& p) : params(p.params) {}

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

      InlineGluonPDFParams params;
    };



} // namespace Chroma

#endif
