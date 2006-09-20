// -*- C++ -*-
// $Id: inline_link_smear.h,v 3.3 2006-09-20 20:28:03 edwards Exp $
/*! \file
 *  \brief Inline link smearing
 */

#ifndef __inline_link_smear_h__
#define __inline_link_smear_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/xml_group_reader.h"

namespace Chroma 
{ 
  /*! \ingroup inlinesmear */
  namespace InlineLinkSmearEnv 
  {
    extern const std::string name;
    bool registerAll();

    
    //! Parameter structure
    struct Params 
    {
      Params();
      Params(XMLReader& xml_in, const std::string& path);
      void writeXML(XMLWriter& xml_out, const std::string& path);

      unsigned long frequency;

      GroupXML_t       link_smearing;        /*!< link smearing xml */

      struct NamedObject_t
      {
	std::string 	gauge_id;	     /*!< Input gauge field */
	std::string 	linksmear_id;	     /*!< Output memory object ape config */

      } named_obj;
      
    };

    //! Inline link smearing
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

};

#endif
