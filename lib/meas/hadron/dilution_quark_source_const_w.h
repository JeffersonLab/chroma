// -*- C++ -*-
// $Id: dilution_quark_source_const_w.h,v 1.15 2008-11-25 22:31:13 kostas Exp $
/*! \file
 * \brief Dilution scheme inferred from pre-generated solutions.
 * 
 * This dilution scheme takes pre-generated LatticeFermion solutions built
 * in the chain of MAKE_SOURCE and PROPAGATOR inline measurement calls
 */

#ifndef __dilution_quark_source_const_h__
#define __dilution_quark_source_const_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "meas/hadron/dilution_scheme.h"
#include "io/qprop_io.h"
#include "io/param_io.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace DilutionQuarkSourceConstEnv 
  {
    extern const std::string name;
    bool registerAll();

    //! Parameter structure
    /*! \ingroup inlinehadron */
    struct Params
    {
      Params();
      Params(XMLReader& xml_in, const std::string& path);
      void writeXML(XMLWriter& xml_out, const std::string& path) const;

      struct QuarkFiles_t
      {
	struct TimeSliceFiles_t 
	{
	  multi1d<std::string> dilution_files;  /*!< dilution files per timeslice*/		
	};
					
	multi1d<TimeSliceFiles_t> timeslice_files; /*!< full time dilution is assumed*/ 
      };
	
      QuarkFiles_t  quark_files;       /*!< All the solutions for a single quark */

      bool UseSourceHeaderSmearing ;

    }; // struct Params

    //! Structure holding a source and its solutions
    struct QuarkSourceSolutions_t
    {
      struct TimeSlices_t
      {
	struct Dilutions_t
	{
	  std::string        soln_file;  //Solution filename for this dilution 
	  PropSourceConst_t  source_header;
	  ChromaProp_t       prop_header;
	};

	int       t0;      //time slice on which this diluted source 
	//has support. Full time dilution is assumed
			
	multi1d<Dilutions_t>  dilutions; 
      };
				
      multi1d<TimeSlices_t> time_slices;

      int   decay_dir;
      Seed  seed;
    };
      
		

    //! Dilution scheme constructed by propagator solutions over diluted MAKE_SOURCE calls 
    class ConstDilutionScheme : public DilutionScheme<LatticeFermion>
    {
    public:

      //! Virtual destructor to help with cleanup;
      ~ConstDilutionScheme() {}
			
      //! Default constructor
      ConstDilutionScheme( const Params& p ) 
	{ 
	  params = p; 
	  init();
	}

      //! The decay direction
      int getDecayDir() const {return quark.decay_dir;}

      //! The seed identifies this quark
      const Seed& getSeed() const {return quark.seed;}

      //! The actual t0 corresponding to this time dilution element 
      int getT0( int t0 ) const {return quark.time_slices[t0].t0;}
      
      //! The number of dilutions per timeslice fo timeslice t0
      int getDilSize( int t0 ) const {return quark.time_slices[t0].dilutions.size();}

      //! The number of dilution timeslices included  
      int getNumTimeSlices() const {return quark.time_slices.size();}

      //! The kappa parameter in the wilson action 
      Real getKappa() const 
	{
	  Real kappa;
	  //Assume the kappa is the same for all dilutions
	  std::istringstream  xml_k(quark.time_slices[0].dilutions[0].prop_header.fermact.xml);
				
	  XMLReader  proptop(xml_k);
	  if ( toBool(proptop.count("/FermionAction/Kappa") != 0) )
	  {
	    read(proptop, "/FermionAction/Kappa", kappa);
	  }
	  else 
	  {
	    Real mass; 
	    read(proptop, "/FermionAction/Mass", mass);
	    kappa = massToKappa(mass);
	  }

	  return kappa;
	}
     
      //! The info from the cfg on which the inversions were performed
      std::string getCfgInfo() const
	{
	  return cfgInfo;
	}

      //! returns the prop header for a given dilution
      std::string getPropHeader(int t0, int dil) const 
	{
	  return quark.time_slices[t0].dilutions[dil].prop_header.fermact.xml;
	}

      //! returns the source header for a given dilution
      std::string getSourceHeader(int t0, int dil) const 
	{
	  return quark.time_slices[t0].dilutions[dil].source_header.source.xml;
	}

      //! Return the diluted source vector
      LatticeFermion dilutedSource(int t0, int dil) const;
    
      //! Return the solution vector corresponding to the diluted source
      LatticeFermion dilutedSolution(int t0, int dil) const;

    protected:
      //! Initialize the object
      void init();
 
      //! Hide partial constructor
      ConstDilutionScheme() {}
			
    private:
      Params params;
      QuarkSourceSolutions_t quark;
      std::string cfgInfo;
    };
    
  } // namespace DilutionQuarkSourceConstEnv


  //! Reader
  /*! @ingroup hadron */
  void read(XMLReader& xml, const string& path, DilutionQuarkSourceConstEnv::Params& param);

  //! Writer
  /*! @ingroup hadron */
  void write(XMLWriter& xml, const string& path, const DilutionQuarkSourceConstEnv::Params& param);

} // namespace Chroma

#endif
