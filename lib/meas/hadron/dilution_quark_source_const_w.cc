// $Id: dilution_quark_source_const_w.cc,v 1.18 2009-03-09 18:23:26 edwards Exp $
/*! \file
 * \brief Dilution scheme specified by MAKE_SOURCE and PROPAGATOR calls  
 *
 */

#include "handle.h"
#include "meas/hadron/dilution_quark_source_const_w.h"
#include "meas/hadron/dilution_scheme_factory.h"
#include "meas/inline/io/named_objmap.h"
#include "meas/sources/dilutezN_source_const.h"
#include "meas/sources/zN_src.h"
#include "util/ft/sftmom.h"


namespace Chroma 
{ 

  // Read parameters
  void read(XMLReader& xml, const string& path, DilutionQuarkSourceConstEnv::Params& param)
  {
    DilutionQuarkSourceConstEnv::Params tmp(xml, path);
    param = tmp;
  }


  // Writer
  void write(XMLWriter& xml, const string& path, const DilutionQuarkSourceConstEnv::Params& param)
  {
    param.writeXML(xml, path);
  }


  /*!
   * \ingroup hadron
   *
   * 
   */
  namespace DilutionQuarkSourceConstEnv
  { 
    //Read Quark dilution files per timeslice
    void read(XMLReader& xml, const string& path, DilutionQuarkSourceConstEnv::Params::QuarkFiles_t::TimeSliceFiles_t& input)
    {
      XMLReader inputtop(xml, path);
      
      read(inputtop, "DilutionFiles", input.dilution_files);
    }

    //Read Quark timeslice files 
    void read(XMLReader& xml, const string& path, DilutionQuarkSourceConstEnv::Params::QuarkFiles_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "TimeSliceFiles", input.timeslice_files);
    }


    //Write Quark dilution files 
    void write(XMLWriter& xml, const string& path, const DilutionQuarkSourceConstEnv::Params::QuarkFiles_t::TimeSliceFiles_t& input)
    {
      push(xml, path);
      write(xml, "dilution_files", input.dilution_files);
      pop(xml);
    }

    //Write Quark time slice files 
    void write(XMLWriter& xml, const string& path, const DilutionQuarkSourceConstEnv::Params::QuarkFiles_t& input)
    {
      push(xml, path);
      write(xml, "timeslice_files", input.timeslice_files);
      pop(xml);
    }


    //! Initialize
    Params::Params()
    {
      UseSourceHeaderSmearing = false ;
    }


    //! Read parameters
    Params::Params(XMLReader& xml, const string& path)
    {
      XMLReader paramtop(xml, path);

      int version;
      read(paramtop, "version", version);

      switch (version) 
      {
      case 1:
	/**************************************************************************/
	break;

      default :
	/**************************************************************************/

	QDPIO::cerr << "Input parameter version " << version << " unsupported." << endl;
	QDP_abort(1);
      }

      UseSourceHeaderSmearing = true ;
      read(paramtop, "QuarkFiles", quark_files);
      if(paramtop.count("UseSourceHeaderSmearing")!=0)
	read(paramtop, "UseSourceHeaderSmearing", UseSourceHeaderSmearing);
    }



    // Writer
    void Params::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);

      int version = 1;
      write(xml, "version", version);
      write(xml, "QuarkFiles", quark_files);
      if(!UseSourceHeaderSmearing)
	write(xml, "UseSourceHeaderSmearing", UseSourceHeaderSmearing);
      pop(xml);
    }

    // Anonymous namespace for registration
    namespace
    {
      DilutionScheme<LatticeFermion>* createScheme(XMLReader& xml_in, 
						   const std::string& path) 
      {
	return new ConstDilutionScheme(Params(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "DILUTION_QUARK_SOURCE_CONST_FERM";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
			
      if (! registered)
      {
	success &= TheFermDilutionSchemeFactory::Instance().registerObject(name, createScheme);
	registered = true;
      }
      return success;
    }



    bool operator!= (const QuarkSourceSolutions_t::TimeSlices_t::Dilutions_t & dilA,
		     const QuarkSourceSolutions_t::TimeSlices_t::Dilutions_t & dilB)
    {
      bool val = false;

      //Check if the seeds are the same 
      Seed seedA, seedB;
		
      std::istringstream  xml_a(dilA.source_header.source.xml);
      XMLReader  rdr_a(xml_a);

      std::istringstream  xml_b(dilB.source_header.source.xml);
      XMLReader  rdr_b(xml_b);

			
      read(rdr_a, "/Source/ran_seed" , seedA);
      read(rdr_b, "/Source/ran_seed" , seedB);

      val |= toBool(seedA != seedB);
      if (val) 
      {
	QDPIO::cerr<< "random seeds are not the same." <<endl;
      }
			
      //Check Spatial mask and spatial mask size
			
      multi1d<int> mask_sizeA, mask_sizeB, cmaskA, cmaskB, smaskA, smaskB;
      multi1d< multi1d<int> > maskA, maskB;

      read(rdr_a, "/Source/spatial_mask_size" , mask_sizeA);
      read(rdr_b, "/Source/spatial_mask_size" , mask_sizeB);
      val |= toBool(mask_sizeA != mask_sizeB);
      if (val)
      {
	QDPIO::cerr<< "spatial mask sizes are not the same." <<endl;
      }

      read(rdr_a, "/Source/spatial_mask" , maskA);
      read(rdr_b, "/Source/spatial_mask" , maskB);
      val |= toBool(maskA != maskB);
      if (val)
      {
	QDPIO::cerr<< "spatial masks are not the same." <<endl;
      }

      read(rdr_a, "/Source/color_mask" , cmaskA);
      read(rdr_b, "/Source/color_mask" , cmaskB);
      val |= toBool(maskA != maskB);
      if (val)
      {
	QDPIO::cerr<< "color masks are not the same." <<endl;
      }

      read(rdr_a, "/Source/spin_mask" , smaskA);
      read(rdr_b, "/Source/spin_mask" , smaskB);
      val |= toBool(maskA != maskB);
      if (val)
      {
	QDPIO::cerr<< "spin masks are not the same." <<endl;
      }

      return val;

    }
    //-------------------------------------------------------------------------------
    // Function call
    void ConstDilutionScheme::init()
    {
      START_CODE();

      StopWatch snoop;
      snoop.reset();
      snoop.start();

      //
      // Read the soultion headers
      //
      StopWatch swatch;
      swatch.reset();
      swatch.start();

      //zN source N
      int N; 

      try
      {

	bool initq = false;
	Real kappa; 

	QDPIO::cout << "Attempt to read solutions " << endl;
				
	quark.time_slices.resize( params.quark_files.timeslice_files.size() );
				
	QDPIO::cout<< "time_slices.size = " << quark.time_slices.size() << endl;

	for (int t0 = 0 ; t0 < quark.time_slices.size() ; ++t0 )
	{
	  quark.time_slices[t0].dilutions.resize(
	    params.quark_files.timeslice_files[t0].dilution_files.size() );
				
	  QDPIO::cout << "dilutions.size = " << 
	    quark.time_slices[t0].dilutions.size() << endl;
					
	  int time;

	  for(int dil = 0; dil < quark.time_slices[t0].dilutions.size(); ++dil)
	  {				
	    quark.time_slices[t0].dilutions[dil].soln_file =
	      params.quark_files.timeslice_files[t0].dilution_files[dil];

	    XMLReader file_xml, record_xml, record_xml_source;

	    QDPIO::cout << "reading file = " << 
	      quark.time_slices[t0].dilutions[dil].soln_file << endl;

	    QDPFileReader from(file_xml, 
			       quark.time_slices[t0].dilutions[dil].soln_file, QDPIO_SERIAL);

	    //Use the record xml only, throw away the lattice fermion
	    LatticeFermion junk;
	    read(from, record_xml, junk);
	    close(from);

	    read(record_xml, "/Propagator/PropSource", 
		 quark.time_slices[t0].dilutions[dil].source_header);
	    read(record_xml, "/Propagator/ForwardProp", 
		 quark.time_slices[t0].dilutions[dil].prop_header);

	    //read the current N 
	    int currN; 
	    read(record_xml, "/Propagator/PropSource/Source/N", currN);

	    if (!initq)
	    {
	      read(record_xml, "/Propagator/PropSource/Source/ran_seed",
		   quark.seed);
				
	      read(record_xml, "/Propagator/PropSource/Source/N", N);
	      quark.decay_dir = 
		quark.time_slices[t0].dilutions[dil].source_header.j_decay;

	      //Test that kappa is the same for all dilutions of this
	      //quark
	      std::istringstream  xml_k(
		quark.time_slices[t0].dilutions[dil].prop_header.fermact.xml);
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
	      		
	      //Test that config is the same for every dilution 
	      XMLReader xml_tmp(record_xml, "/Propagator/Config_info");
	      std::ostringstream os;
	      xml_tmp.print(os);

	      cfgInfo = os.str();

	      initq = true;
	    }

	    Real kappa2;
	    std::istringstream  xml_k2(
	      quark.time_slices[t0].dilutions[dil].prop_header.fermact.xml);
	    XMLReader  proptop2(xml_k2);
						
	
	    if ( toBool(proptop2.count("/FermionAction/Kappa") != 0) )
	    {
	      read(proptop2, "/FermionAction/Kappa", kappa2);
	    }
	    else 
	    {
	      Real mass; 
	      read(proptop2, "/FermionAction/Mass", mass);
	      kappa2 = massToKappa(mass);
	    }
						
	    if ( toBool(kappa != kappa2) )
	    {
	      QDPIO::cerr << "Kappa is not the same for all dilutions: t0=" <<
		t0 << " dil= "<< dil << " kappa2 = "<< kappa2 << endl;
								
	      QDP_abort(1);
	    }

	    if (currN != N)
	    {
	      QDPIO::cerr << "N is not the same for all dilutions: t0 = " <<
		t0 << " dil = " << dil << endl;
	    }

	    //Test that t0 is the same for all dilutions on this timeslice
	    //grab the first t0 to prime the process
					
	    if (dil == 0)
	    {
	      time = quark.time_slices[t0].dilutions[dil].source_header.t_source;
	    }

	    if (time != quark.time_slices[t0].dilutions[dil].source_header.t_source)
	    {
	      QDPIO::cerr << "t0's DO NOT MATCH FOR ALL DILUTIONS ON TIME SLICE "
			  << t0 << endl;

	      QDP_abort(1);
	    }

					
	    //Test that each t0 has the same dilutions per timeslice 
	    if ( quark.time_slices[t0].dilutions[dil] != 
		 quark.time_slices[0].dilutions[dil] )
	    {
	      QDPIO::cerr << "Dilutions do not match on time slice " << 
		t0 << " dil = "<< dil<< endl;

	      QDP_abort(1);
	    }

	    //Test that this dilution element was created on correct cfg
	    std::string currCfgInfo;
	    {
	      XMLReader xml_tmp(record_xml, "/Propagator/Config_info");
	      std::ostringstream os;
	      xml_tmp.print(os);

	      currCfgInfo = os.str();
	    }

	    if (cfgInfo != currCfgInfo)
	    {
	      QDPIO::cerr << "Cfgs do not match on time slice " << 
		t0 << " dil = "<< dil<< endl;

	      QDP_abort(1);
	    }
					
	  }//dil
					
	  quark.time_slices[t0].t0 = time;

	}//t0

			
#if 0
#warning "Turned off the sanity check that the dilutions summed to a unity operator on a time-slice"
	//Ensure that the given dilutions form a full dilution scheme per 
	//timeslice. Only need to check a single timeslice as 
	//we have guaranteed the same dilutions per timeslice
			
	LatticeFermion quark_noise;      // noisy source on entire lattice
	QDP::RNG::setrn(quark.seed);
	zN_src(quark_noise, N);

			
	for (int dil = 0 ; dil < quark.time_slices[0].dilutions.size() ; ++dil)
	{
	  LatticeFermion source = dilutedSource(0, dil);
	  quark_noise -= source; 
	}
	
	SftMom phases_nomom(0, true, 
			    quark.time_slices[0].dilutions[0].source_header.j_decay);

	Double dcnt = norm2(quark_noise, 
			    phases_nomom.getSet()[quark.time_slices[0].t0] );

	if (toDouble(dcnt) != 0.0)
	{
	  QDPIO::cerr << "Not a complete dilution scheme per timeslice" <<
	    endl;

	  QDP_abort(1);
	}
#endif
				
      } //try
      catch (const string& e) 
      {
	QDPIO::cerr << "Error extracting headers: " << e << endl;
	QDP_abort(1);
      }
      swatch.stop();

      QDPIO::cout << "Source and solution headers successfully read: time = "
		  << swatch.getTimeInSeconds() 
		  << " secs" << endl;


      END_CODE();
    } // init


		
    //Create and return the diluted source
    //For gauge invariance testing reasons, this routine could be changed
    //to get the source from the named object map
    LatticeFermion ConstDilutionScheme::dilutedSource(int t0, int dil ) const 
    {
      const QuarkSourceSolutions_t::TimeSlices_t::Dilutions_t &qq = 
	quark.time_slices[t0].dilutions[dil];

/*
//For now read the source 
LatticeFermion sour;

std::string filename = 
params.quark_files.timeslice_files[t0].dilution_files[dil];

filename.erase(0,9);

std::string source_filename = "zN_source" + filename;

XMLReader file_xml, record_xml;

QDPIO::cout << "reading file = " << source_filename << endl;
QDPFileReader from(file_xml, source_filename, QDPIO_SERIAL);

read(from, record_xml, sour);
*/

			
      //Dummy gauge field to send to the source construction routine;
      multi1d<LatticeColorMatrix> dummy; 
			

      // Build source construction
      QDPIO::cout << "Source_id = " << qq.source_header.source.id << endl;
      QDPIO::cout << "Source = XX" << qq.source_header.source.xml << "XX" << endl;

      std::istringstream  xml_s(qq.source_header.source.xml);
      XMLReader  sourcetop(xml_s);

      if (qq.source_header.source.id != DiluteZNQuarkSourceConstEnv::getName())
      {
	QDPIO::cerr << "Expected source_type = " << DiluteZNQuarkSourceConstEnv::getName() << endl;
	QDP_abort(1);
      }

      // Manually create the params so I can peek into them and use the source constructor
      DiluteZNQuarkSourceConstEnv::Params  srcParams(sourcetop, 
						     qq.source_header.source.path);

      if(!params.UseSourceHeaderSmearing){
	QDPIO::cout<<name
		   <<": Will NOT use Smearing/displacement options specified in the header\n"<<endl ;
	srcParams.smear = false ; 
      }
      else{
	QDPIO::cout<<name<<": Smearing/displacement not implemented\n"<<endl ;
	QDPIO::cout<<name
		   <<": Will NOT use Smearing/displacement options specified in the header\n"<<endl ;
	srcParams.smear = false ; 
      }

      DiluteZNQuarkSourceConstEnv::SourceConst<LatticeFermion>  srcConst(srcParams);
      
      QDP::RNG::setrn( quark.seed );

     
      LatticeFermion sour = srcConst(dummy);
			
			
      /*
	multi1d<int> orig(4);
	for (int i = 0 ; i < 4 ; ++i)
	{
	orig=0;
	}

	LatticeColorVector vtr = peekSpin(sour, 2);
	LatticeComplex comp = peekColor( vtr , 0 );
	QDPIO::cout<< "Sourceval = "<< 
	peekSite( comp, orig ) << endl;
      */

      return sour;
	      
    } //dilutedSource		


    LatticeFermion ConstDilutionScheme::dilutedSolution(int t0, int dil) const 
    {
			
      const std::string & filename = quark.time_slices[t0].dilutions[dil].soln_file;

      LatticeFermion soln; 

      XMLReader file_xml, record_xml;

      QDPIO::cout << "reading file = " << filename << endl;
      QDPFileReader from(file_xml, filename, QDPIO_SERIAL);

      read(from, record_xml, soln);

/*
  multi1d<int> orig(4);
  for (int i = 0 ; i < 4 ; ++i)
  {
  orig=0;
  }

  LatticeColorVector vtr = peekSpin(soln, 2);
  LatticeComplex comp = peekColor( vtr , 0 );
  QDPIO::cout<< "Sinkval = "<< 
  peekSite( comp, orig ) << endl;

*/
      return soln;
    }


	
  } // namespace DilutionQuarkSourceConstEnv

}// namespace Chroma 
