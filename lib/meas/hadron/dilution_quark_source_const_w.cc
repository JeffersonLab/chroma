// $Id: dilution_quark_source_const_w.cc,v 1.2 2007-12-18 13:40:25 edwards Exp $
/*! \file
 * \brief Inline measurement of stochastic group baryon operator
 *
 */

#include "handle.h"
#include "meas/hadron/dilution_quark_source_const_w.h"
#include "meas/hadron/dilution_operator_factory.h"
#include "meas/inline/io/named_objmap.h"

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
    //Read time dilution files 
    void read(XMLReader& xml, const string& path, DilutionQuarkSourceConstEnv::Params::QuarkFiles_t::TimeDilutions_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "dilution_files", input.dilution_files);
    }


    //Write time dilution files 
    void write(XMLWriter& xml, const string& path, const DilutionQuarkSourceConstEnv::Params::QuarkFiles_t::TimeDilutions_t& input)
    {
      push(xml, path);
      write(xml, "dilution_files", input.dilution_files);
      pop(xml);
    }



    //Read Quark dilution files 
    void read(XMLReader& xml, const string& path, DilutionQuarkSourceConstEnv::Params::QuarkFiles_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "time_files", input.time_files);
    }


    //Write Quark dilution files 
    void write(XMLWriter& xml, const string& path, const DilutionQuarkSourceConstEnv::Params::QuarkFiles_t& input)
    {
      push(xml, path);
      write(xml, "time_files", input.time_files);
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
      DilutionOperator<LatticeFermion>* createOperator(XMLReader& xml_in, 
						       const std::string& path) 
      {
	return new Dilute(Params(xml_in, path));
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
	success &= TheFermDilutionOperatorFactory::Instance().registerObject(name, createOperator);
	registered = true;
      }
      return success;
    }



    //--------------------------------------------------------------
    //! Return the time dilution element to which a particular time belongs
    bool Dilute::hasTimeSupport(const_iterator iter, int time) const
    {
      int ret_t0 = 0;
      bool stop = false; 
				
      for ( int t0 = 0 ; t0 < quark.time_dilutions.size() ;  ++t0 )
      {
	for ( int time0 = 0 ; time0 < quark.time_dilutions[t0].t0.size() ; ++time0 )
	{
	  if ( time == quark.time_dilutions[t0].t0[ time0 ] ) 
	  {	
	    ret_t0 = t0;	
	    stop = true;
	    break;
	  }
	}
					
	if (stop)
	{
	  break;
	}
				
      }
			
      return ret_t0;
    }			
	


    //-------------------------------------------------------------------------------
    // Function call
    void Dilute::init()
    {
      START_CODE();

      StopWatch snoop;
      snoop.reset();
      snoop.start();

      //
      // Read the source and solutions
      //
      StopWatch swatch;
      swatch.reset();
      swatch.start();

      try
      {
	bool initq = false;

	QDPIO::cout << "Attempt to read solutions " << endl;
		
	quark.time_dilutions.resize(params.quark.soln_files.size());

	QDPIO::cout << "time_dilutions.size= " << quark.time_dilutions.size() << endl;
	for(int t=0; t < quark.time_dilutions.size(); ++t)
	{
	  quark.time_dilutions[t].dilutions.resize(params.quark.soln_files[t].dilution_files.size());
	  QDPIO::cout << "dilutions.size= " << quark.time_dilutions[t].dilutions.size() << endl;
	  for(int i=0; i < quark.time_dilutions[t].dilutions.size(); ++i)
	  {
	    // Short-hand
	    const std::string& dilution_file = params.quark.soln_files[t].dilution_files[i];

	    QuarkSourceSolutions_t::TimeSlices_t::Dilutions_t& qq = 
	      quark.time_dilutions[t].dilutions[i];

	    XMLReader file_xml, record_xml, record_xml_source;

	    QDPIO::cout << "reading file= " << dilution_file << endl;
	    QDPFileReader from(file_xml, dilution_file, QDPIO_SERIAL);

	    //For now, read both source and solution
	    read(from, record_xml, qq.soln);
	    //  read(from, record_xml_source, qq.source);
	    close(from);

	
	    read(record_xml, "/Propagator/PropSource", qq.source_header);
	    read(record_xml, "/Propagator/ForwardProp", qq.prop_header);

	    //Get source from the named object map 
	    std::stringstream srcstrm;
	    srcstrm << 	"zN_source_q" << n + 1 << "_t" << 
	      qq.source_header.t_source;

	    std::string source_name = srcstrm.str();

	    qq.source = TheNamedObjMap::Instance().getData< LatticeFermion >(source_name);
	    
			
	    if (!initq)
	    {
	      read(record_xml, "/Propagator/PropSource/Source/ran_seed",
		   quark.seed);
					
	      initq = true;
	    }

	    qq.t0 = qq.source_header.t_source;
	    j_decay = qq.source_header.j_decay;
	   
	  }
	}
      } //try
      catch (const string& e) 
      {
	QDPIO::cerr << "Error extracting headers: " << e << endl;
	QDP_abort(1);
      }
      swatch.stop();

      QDPIO::cout << "Sources and solutions successfully read: time= "
		  << swatch.getTimeInSeconds() 
		  << " secs" << endl;


      END_CODE();
    } // init

  } // namespace DilutionQuarkSourceConstEnv

  /*! @} */  // end of group hadron

} // namespace Chroma
