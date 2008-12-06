// -*- C++ -*-
// $Id: inline_disco_w.h,v 1.1 2008-12-06 01:00:52 kostas Exp $
/*! \file
 * \brief Inline measurement of stochastic 3pt functions.
 *
 * spectroscopy
 */

#ifndef __inline_disco_h__
#define __inline_disco_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"
//#include <map>

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineDiscoEnv 
  {
    extern const std::string name;
    bool registerAll();
    
    // The flavors
    //! Parameter structure
    /*! \ingroup inlinehadron */
    struct Params 
    {
      Params();
      Params(XMLReader& xml_in, const std::string& path);
       
      unsigned long      frequency;
    
      struct Param_t
      {
	int max_path_length ; /*! maximum displacement path */
	int p2_max ; /*! maximum p2  */
	multi1d<GroupXML_t> chi ;     /*! dilutions */
	string mass_label ; /*! a string flag maybe used in analysis*/
      } param;
    
      struct NamedObject_t
      {
	std::string         gauge_id;
	std::string         op_db_file;
      } named_obj;
      
      std::string xml_file;  // Alternate XML file pattern

      void write(XMLWriter& xml_out, const std::string& path);

    };
  

  //! Inline measurement of stochastic baryon operators
  /*! \ingroup inlinehadron */
    class InlineMeas : public AbsInlineMeasurement{
    protected:
      //! Do the measurement
      void func(const unsigned long update_no,
		XMLWriter& xml_out); 
      
    private:
 

      Params params;

    public:
      ~InlineMeas() {}
      InlineMeas(const Params& p) : params(p) {}
      InlineMeas(const InlineMeas& p) : params(p.params) {}
      
      unsigned long getFrequency(void) const {return params.frequency;}
      
      //! Do the measurement
      void operator()(const unsigned long update_no,
		      XMLWriter& xml_out); 
      
    };

  }; // name space InlineDiscoEnv

};

#endif
