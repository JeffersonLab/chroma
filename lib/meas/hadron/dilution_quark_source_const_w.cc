// $Id: dilution_quark_source_const_w.cc,v 1.3 2008-01-07 15:19:43 jbulava Exp $
/*! \file
 * \brief Dilution scheme specified by MAKE_SOURCE and PROPAGATOR calls  
 *
 */

#include "handle.h"
#include "meas/hadron/dilution_quark_source_const_w.h"
#include "meas/hadron/dilution_scheme_factory.h"
#include "meas/inline/io/named_objmap.h"
#include "meas/sources/dilutezN_source_const.h"

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
   * @{
   */
  namespace DilutionQuarkSourceConstEnv
  { 
    //Read Quark timeslice files 
    void read(XMLReader& xml, const string& path, DilutionQuarkSourceConstEnv::Params::QuarkFiles_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "timeslice_files", input.timeslice_files);
    }

		//Read Quark dilution files per timeslice
    void read(XMLReader& xml, const string& path, DilutionQuarkSourceConstEnv::Params::QuarkFiles_t::TimeSliceFiles_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "dilution_files", input.dilution_files);
    }


    //Write Quark time slice files 
    void write(XMLWriter& xml, const string& path, const DilutionQuarkSourceConstEnv::Params::QuarkFiles_t& input)
    {
      push(xml, path);
      write(xml, "timeslice_files", input.timeslice_files);
      pop(xml);
    }

 //Write Quark dilution files 
    void write(XMLWriter& xml, const string& path, const DilutionQuarkSourceConstEnv::Params::QuarkFiles_t::TimeSliceFiles_t& input)
    {
      push(xml, path);
      write(xml, "dilution_files", input.dilution_files);
      pop(xml);
    }

    //! Initialize
    Params::Params()
    {
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

      read(paramtop, "Quark", quark);
    }


    // Writer
    void Params::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);

      int version = 1;
      write(xml, "version", version);
      write(xml, "Quark", quark);
      pop(xml);
    }


  }


  namespace DilutionQuarkSourceConstEnv 
  { 
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



    //-------------------------------------------------------------- 
/*    bool ConstDilutionScheme::hasSupport(int t0, int dil) const
    {
    
			return (quark.[dil].t0 == t0); 
		
		}			
*/	


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

      try
      {

	bool initq = false;

	QDPIO::cout << "Attempt to read solutions " << endl;
		
	quark.dilutions.resize(params.quark.dilution_files.size());

	QDPIO::cout << "dilutions.size = " << quark.dilutions.size() << endl;
	
	for(int dil = 0; dil < quark.dilutions.size(); ++dil)
	{
			quark.dilutions[dil].soln_file = params.quark.dilution_files[dil];

			
	    XMLReader file_xml, record_xml, record_xml_source;

	    QDPIO::cout << "reading file = " << quark.dilutions[dil].soln_file << endl;
	    QDPFileReader from(file_xml, quark.dilutions[dil].soln_file, QDPIO_SERIAL);

			//Read the record xml only
			read(from, record_xml);
	    close(from);

	    read(record_xml, "/Propagator/PropSource", quark.dilutions[dil].source_header);
	    read(record_xml, "/Propagator/ForwardProp", quark.dilutions[dil].prop_header);


	    if (!initq)
	    {
	      read(record_xml, "/Propagator/PropSource/Source/ran_seed",
		   quark.seed);
					
	    	quark.dilutions[dil].decay_dir = quark.dilutions[dil].source_header.j_decay;
	      initq = true;
	    }

	    quark.dilutions[dil].t0 = quark.dilutions[dil].source_header.t_source;
	   
	  }//dil
      } //try
      catch (const string& e) 
      {
	QDPIO::cerr << "Error extracting headers: " << e << endl;
	QDP_abort(1);
      }
      swatch.stop();

      QDPIO::cout << "Source and solution headers successfully read: time= "
		  << swatch.getTimeInSeconds() 
		  << " secs" << endl;


      END_CODE();
    } // init


		
		//create and return the dilution source
		//For gauge invariance testing reasons, this routine could be changed
		//to get the source from the named object map
		LatticeFermion ConstDilutionScheme::dilutedSource( int dil ) const 
		{
			QuarkSourceSolutions_t::Dilutions_t & qq = quark.dilutions[dil];

			//Dummy gauge field to send to source construction;
			multi1d<LatticeColorMatrix> dummy = zero;
			

			// Build source construction
	    QDPIO::cout << "Source_id = " << qq.source_header.source.id << endl;
	    QDPIO::cout << "Source = XX" << qq.source_header.source.xml << "XX" << endl;

	    std::istringstream  xml_s(qq.source_header.source.xml);
	    XMLReader  sourcetop(xml_s);

	    if (qq.source_header.source.id != DiluteZNQuarkSourceConstEnv::name)
	    {
				QDPIO::cerr << "Expected source_type = " << DiluteZNQuarkSourceConstEnv::name << endl;
				QDP_abort(1);
	    }

	    // Manually create the params so I can peek into them and use the source constructor
	    DiluteZNQuarkSourceConstEnv::Params  srcParams(sourcetop, 
							     qq.source_header.source.path);
	    DiluteZNQuarkSourceConstEnv::SourceConst<LatticeFermion>  srcConst(srcParams);
      
			QDP::RNG::setrn( quark[dil].seed );

			 
			return srcConst(dummy);
	      
		} //dilutedSource		


		LatticeFermion ConstDilutionScheme::dilutedSolution(int dil) const 
		{
			
			const std::string & filename = quark.dilutions[dil].soln_file;

			LatticeFermion soln; 

		  XMLReader file_xml, record_xml, record_xml_source;

	    QDPIO::cout << "reading file = " << quark.dilutions[dil].soln_file << endl;
	    QDPFileReader from(file_xml, quark.dilutions[dil].soln_file, QDPIO_SERIAL);

			read(from, record_xml, soln);

			return soln;
		}


	
  } // namespace DilutionQuarkSourceConstEnv


