// -*- C++ -*-
/*! \file
 * \brief Inline measurement of glueball operators
 */

#ifndef __inline_glueball_ops_h__
#define __inline_glueball_ops_h__

#include "meas/inline/abs_inline_measurement.h"
#include "io/xml_group_reader.h"

namespace Chroma 
{ 
  /*! \ingroup inlineglue */
  namespace InlineGlueballOpsEnv 
  {
    bool registerAll();


    //! Parameter structure
    /*! \ingroup inlineglue */
    struct Params
    {
      Params();
      Params(XMLReader& xml_in, const std::string& path);
      void writeXML(XMLWriter& xml_out, const std::string& path) const;

      unsigned long      frequency;

      struct Param_t
      {
	int                     mom2_max;               /*!< (mom)^2 <= mom2_max */
	int                     displacement_length;    /*!< Displacement length for creat. and annih. ops */
	int                     decay_dir;              /*!< Decay direction */
	multi1d< multi1d<int> > displacement_list;      /*!< Array of displacements list to generate */
	GroupXML_t              link_smearing;          /*!< link smearing xml */
      };

      struct NamedObject_t
      {
	std::string         gauge_id;               /*!< Gauge field */
	std::string         glue_op_file;           /*!< File name for creation operators */
      };

      Param_t        param;      /*!< Parameters */    
      NamedObject_t  named_obj;  /*!< Named objects */
      std::string    xml_file;   /*!< Alternate XML file pattern */
    };


    //! Inline measurement of glueball operators
    /*! \ingroup inlineglue */
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

  } // namespace InlineGlueballOpsEnv 
}

#endif
