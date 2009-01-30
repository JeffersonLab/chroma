// -*- C++ -*-
// $Id: inline_disco_eigcg_w.h,v 1.5 2009-01-30 03:42:39 kostas Exp $
/*! \file
 * \brief Inline measurement of stochastic 3pt functions.
 *
 * spectroscopy
 */

#ifndef __inline_disco_eigcg_h__
#define __inline_disco_eigcg_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"
//#include <map>

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineDiscoEigCGEnv 
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
	GroupXML_t action ; /*! need to specify the action whose operator D we are computing the trace. The vectors may require manipulation if they come from an EO action and or if they are vectors that approximate the D^\dagger D invese (as those that EigCG produces are */
      } param;
    
      struct NamedObject_t
      {
	std::string         gauge_id    ;
	std::string         evecs_file  ;
	std::string         op_db_file  ;
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

  }; // name space InlineDiscoEigCGEnv

};

#endif
