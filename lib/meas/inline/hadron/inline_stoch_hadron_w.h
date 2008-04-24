// -*- C++ -*-
// $Id: inline_stoch_hadron_w.h,v 1.5 2008-04-24 05:20:17 kostas Exp $
/*! \file
 * \brief Inline measurement of stochastic hadron operator (mesons and baryons).
 *
 * spectroscopy
 */

#ifndef __inline_stoch_hadron_h__
#define __inline_stoch_hadron_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"
//#include <map>

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineStochHadronEnv 
  {
    extern const std::string name;
    bool registerAll();
    
    // The flavors
    enum Flavor {up, down, strange, charm, bottom};

    //! Parameter structure
    /*! \ingroup inlinehadron */
    struct Params 
    {
      Params();
      Params(XMLReader& xml_in, const std::string& path);
       
      unsigned long      frequency;
    
      struct Param_t
      {
	int          mom2_max;              /*!< (mom)^2 <= mom2_max */
	
	//Operators can be mesons or baryons (quark-antiquark or three quarks).
	multi1d<GroupXML_t> ops ; /*!< array with hadron operators xml */
	//Each operator needs its own output file
	//should be defined in the GroupXML_t
      
	multi1d<GroupXML_t>  smearing; /*!< xml holding smearing params */
	multi1d<GroupXML_t>  displace; /*!< xml holding displacement params */
	GroupXML_t   link_smear;         /*!< link smearing xml one for all*/

	multi1d<GroupXML_t> quarks ;     /*! dilutions */
	
      } param;
    
    
      struct NamedObject_t
      {
	std::string         gauge_id;
      } named_obj;
      
      std::string xml_file;  // Alternate XML file pattern

      void write(XMLWriter& xml_out, const std::string& path);

    };
  

  //! Inline measurement of stochastic baryon operators
  /*! \ingroup inlinehadron */
    class InlineMeas : public AbsInlineMeasurement{
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

  }; // name space InlineHadronEnv

};

#endif
