// -*- C++ -*-
/*! \file
 * \brief Inline construction of QpropMatMul
 *
 * QpropMatMul calculations
 */

#ifndef __inline_qprop_matmul_w_h__
#define __inline_qprop_matmul_w_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineQpropMatMulEnv
  {
    extern const std::string name;
    bool registerAll();
  }

  //! Parameter structure
  /*! \ingroup inlinehadron */ 
  struct InlineQpropMatMulParams
  {
    InlineQpropMatMulParams();
    InlineQpropMatMulParams(XMLReader& xml_in, const std::string& path);
    void writeXML(XMLWriter& xml_out, const std::string& path);

    unsigned long     frequency;


    GroupXML_t      fermact;

    struct NamedObject_t
    {
      std::string     gauge_id;
      std::string     header_id;
      std::string     source_id;
      std::string     result_id;
    } named_obj;

    std::string xml_file;  // Alternate XML file pattern
  };

  //! Inline QpropMatMul calculation
  /*! \ingroup inlinehadron */
  class InlineQpropMatMul : public AbsInlineMeasurement
  {
  public:
    ~InlineQpropMatMul() {}
    InlineQpropMatMul(const InlineQpropMatMulParams& p) : params(p) {}
    InlineQpropMatMul(const InlineQpropMatMul& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out); 

  protected:
    //! Do the measurement
    void func(const unsigned long update_no,
	      XMLWriter& xml_out); 

  private:
    InlineQpropMatMulParams params;
  };

}

#endif
