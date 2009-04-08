// -*- C++ -*-
// $Id: inline_disco_eoprec_w.h,v 3.1 2009-04-08 18:34:11 caubin Exp $
/*! \file
 * \brief Inline measurement of stochastic 3pt functions.
 *
 * spectroscopy
 */

#ifndef __inline_disco_eoprec_h__
#define __inline_disco_eoprec_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"
//#include <map>

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineDiscoEOPrecEnv 
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
	GroupXML_t action ; /*! need to specify the action whose operator D we are computing the trace. \
The vectors may require manipulation if they come from an EO action and or if they are vectors that appr\
oximate the D^\dagger D invese (as those that EigCG produces are */

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

    //! Function object used for constructing the time-slice set but only with odd sites                 
    //! Has an index ultimately such that                                                                
    //! ind = 2*t + rb                                                                                   
    //! where t is the timeslice we want and rb is 0 (even) or 1 (odd), telling us what                  
    //! kind of site this will be.                                                                       
    class TimeSliceRBFunc : public SetFunc
    {
    public:
      TimeSliceRBFunc(int dir): dir_decay(dir) {}

      int operator() (const multi1d<int>& coordinate) const
      {
        int sum = 0;
        for(int m=0; m < coordinate.size(); ++m)
          sum += coordinate[m];

        if ((dir_decay<0)||(dir_decay>=Nd)) {
          return sum & 1 ; //In this case return the entire rb subset, with all timeslices               
        }
        else {
          return ( 2*(coordinate[dir_decay]) + (sum & 1) ) ;
        }
      }

      int numSubsets() const
      {
        if ((dir_decay<0)||(dir_decay>=Nd)) {
          return 2 ; //Here there's the even or the odd subset                                           
        } else {
          // There are this many subsets, because for each time direction, we have                       
          // either even or odd sites...                                                                 
          return 2*(Layout::lattSize()[dir_decay]) ;
        }
      }

    private:
      TimeSliceRBFunc() {}  // hide default constructor                                                  

      int dir_decay;
    };



  }; // name space InlineDiscoEOPrecEnv

};

#endif
