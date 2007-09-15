// -*- C++ -*-
// $Id: inline_stoch_hadron_w.h,v 1.2 2007-09-15 02:39:25 kostas Exp $
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
  }

  // The flavors
  enum Flavor {up, down, strange, charm, bottom};

  //! Parameter structure
  /*! \ingroup inlinehadron */
  struct InlineStochHadronParams 
  {
    InlineStochHadronParams();
    InlineStochHadronParams(XMLReader& xml_in, const std::string& path);
    void write(XMLWriter& xml_out, const std::string& path);
    
    unsigned long      frequency;
    
    struct Param_t
    {
      int          mom2_max;              /*!< (mom)^2 <= mom2_max */
      
      //map<const char *, GroupXML_t> ops /*!< map with hadron operators xml */
      // Allow for many operators at a time.
      // Operators can be mesons or baryons (quark-antiquark or three quarks).
      multi1d<GroupXML_t> ops ; /*!< array with hadron operators xml */
      //Each operator needs its own output file
      //should be defined in the GroupXML_t
      
      GroupXML_t   source_quark_smearing; /*!< xml holding smearing params */
      GroupXML_t   sink_quark_smearing;   /*!< xml holding smearing params */
      GroupXML_t   link_smearing;         /*!< link smearing xml */
      
    } param;
    
    struct Flavor_t{
      //! Flavor 0: up, 1: down, 2: strange, 3: charm, 4: bottom 
      struct Dilutions_t{
	multi1d<std::string> soln_files;
      };
    
      multi1d<Dilutions_t>  flavor;
    };
    
    //this is not needed. The output file leaves in the operator XML
    //std::string          op_file; // need to figure out how to write out
    
  };
  
  struct NamedObject_t
  {
    multi1d<Flavors_t>   flavors;
    std::string          gauge_id;
  } named_obj;
  
  std::string xml_file;  // Alternate XML file pattern
};


  //! Inline measurement of stochastic baryon operators
  /*! \ingroup inlinehadron */
  class InlineStochHadron : public AbsInlineMeasurement 
  {
  public:
    ~InlineStochHadron() {}
    InlineStochHadron(const InlineStochHadronParams& p) : params(p) {}
    InlineStochHadron(const InlineStochHadron& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out); 

  protected:
    //! Do the measurement
    void func(const unsigned long update_no,
	      XMLWriter& xml_out); 

  private:
    InlineStochHadronParams params;
  };

}

#endif
