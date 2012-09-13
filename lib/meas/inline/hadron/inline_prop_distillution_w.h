// -*- C++ -*-
/*! \file
 * \brief Compute the propagator from distillution
 *
 * Propagator calculation in distillution
 */

#ifndef __inline_prop_distillution_w_h__
#define __inline_prop_distillution_w_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"
#include "io/xml_group_reader.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlinePropDistillutionEnv 
  {
    bool registerAll();

    //! Parameter structure
    /*! \ingroup inlinehadron */ 
    struct Params 
    {
      Params();
      Params(XMLReader& xml_in, const std::string& path);

      unsigned long     frequency;

      struct Param_t
      {
	struct Contract_t
	{
	  std::vector<int>    quark_lines;        /*!< Quark line numbers */
	  std::string         mass;               /*!< Some kind of mass label */
	  GroupXML_t          quark_line_xml;     /*!< Quark line xml */
	};

	Contract_t      contract;
	ChromaProp_t    prop;
      };

      struct NamedObject_t
      {
	bool            save_solnP;         /*!< Save solutions */
	std::string     gauge_id;           /*!< Gauge field */
	std::string     distillution_id;    /*!< Distillution factory */
	std::string     src_file;           /*!< File output propagator sources */
	std::string     soln_file;          /*!< File output propagator solutions */
      };


      Param_t           param;
      NamedObject_t     named_obj;
      std::string       xml_file;           /*!< Alternate XML file pattern */
    };


    //! Inline task for the propagator from distillution
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

  } // namespace PropColorVec


}

#endif
