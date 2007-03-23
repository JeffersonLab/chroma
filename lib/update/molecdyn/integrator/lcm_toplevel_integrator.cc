#include "chromabase.h"
#include "meas/inline/io/named_objmap.h"
#include "update/molecdyn/integrator/lcm_toplevel_integrator.h"
#include "update/molecdyn/integrator/md_integrator_factory.h"
#include "update/molecdyn/integrator/lcm_integrator_leaps.h"

namespace Chroma { 

  //! A reader for an element of the copy list
  void read(XMLReader& xml, const std::string& path, IDPair& p)
  {
    try {
      XMLReader paramtop(xml, path);
      read(paramtop, "./copyFrom", p.source);
      read(paramtop, "./copyTo", p.dest);
    }
    catch( const std::string& e) { 
      QDPIO::cout << "Caught exception reading XML: " << e << endl;
      QDP_abort(1);
    }
  }

  //! A Writer for an element of the copy list
  void write(XMLWriter& xml, const std::string& path, const IDPair& p)
  {
    push(xml, path);
    write(xml, "copyFrom", p.source);
    write(xml, "copyTo", p.dest);
    pop(xml);
  }

  LCMToplevelIntegratorParams::    LCMToplevelIntegratorParams(XMLReader& xml, const std::string& path) {
      try { 
	XMLReader paramtop(xml, path);
	read(paramtop, "tau0", tau0);           // Read Traj Length
	copy_list.resize(0);                    // Initialize 0 length list
	if( paramtop.count("./copyList") > 0 ) { 
	  read(paramtop, "./copyList", copy_list);  // Read Copy List
	}

	//Read the integrator XML
	XMLReader integrator_xml_reader(paramtop, "./Integrator");
	std::ostringstream os;
	integrator_xml_reader.print(os);
	integrator_xml = os.str();

	// Look for anisotropic integration.
	// (Optional)
	if( paramtop.count("anisoP") == 1 ) { 
	  read(paramtop, "anisoP", anisoP);
	  read(paramtop, "t_dir", t_dir);
	  if( t_dir < 0  || t_dir >= Nd ) { 
	    QDPIO::cout << "Value of t_dir must be 0 <= t_dir < Nd. t_dir is " << t_dir << endl;
	    QDP_abort(1);
	  }
	  read(paramtop, "xi_mom", xi_mom);
	}
	else { 
	  // If there is no data then set defaults for isotropy
	  anisoP = false;
	  t_dir = 3;
	  xi_mom = 1;
	}

      }
      catch(const std::string& e) { 
	QDPIO::cout << "Caught Exception Reading XML: " << e << endl;
	QDP_abort(1);
      }

      
  }

  //! Read the Integrator Params
  void read(XMLReader& xml, const std::string& path, LCMToplevelIntegratorParams& p) 
  {
    LCMToplevelIntegratorParams tmp(xml, path);
    p = tmp;
  }
  
  //! Write the Integrator Params
  void write(XMLWriter& xml, const std::string& path, const LCMToplevelIntegratorParams& p) 
  {
    push(xml, path);
    write(xml, "tau0", p.tau0);
    write(xml, "copyList", p.copy_list);
    std::istringstream int_is(p.integrator_xml);
    XMLReader int_reader(int_is);
    xml << int_reader;  
    if( p.anisoP ) {
      write(xml, "anisoP", p.anisoP);
      write(xml, "t_dir", p.t_dir);
      write(xml, "xi_mom", p.xi_mom);
    }
    pop(xml);
  }


  LCMToplevelIntegrator::LCMToplevelIntegrator(
					       const LCMToplevelIntegratorParams& p) :
      params(p) {

      std::istringstream is( p.integrator_xml );
      XMLReader integrator_reader( is );
      std::string root="/Integrator";
      std::string integrator_name;
      try { 
	read(integrator_reader, "/Integrator/Name", integrator_name);
      }
      catch( const std::string& e) {
	QDPIO::cout << "Caught Exception while processing XML: " << e << endl;
	QDP_abort(1);
      }
    
      Handle< AbsComponentIntegrator< multi1d<LatticeColorMatrix>,
	                              multi1d<LatticeColorMatrix> > > h(
									TheMDComponentIntegratorFactory::Instance().createObject(integrator_name, integrator_reader, root));
      
      top_integrator = h;
      
      // Deal with Anisotropic integration
      if ( p.anisoP ) { 
	Real factor = Real(1) / p.xi_mom;
	// Set the step size
	LCMMDIntegratorSteps::theAnisoStepSizeArray::Instance().setAnisoStepSize(p.t_dir, factor);

      }
  }
    
  void LCMToplevelIntegrator::copyFields(void) const { 
      typedef Handle< Monomial< multi1d<LatticeColorMatrix>, 
	multi1d<LatticeColorMatrix> > > MHandle;

      QDPIO::cout << "Working through Copy List of length " << params.copy_list.size() << endl;
      // Loop through all the items in the copy list
      for(int i=0; i < params.copy_list.size(); i++) { 
	std::string source_id = params.copy_list[i].source;
	std::string dest_id = params.copy_list[i].dest;
	MHandle src_mon = 
	  TheNamedObjMap::Instance().getData<MHandle>(source_id);

	MHandle dest_mon = 
	  TheNamedObjMap::Instance().getData<MHandle>(dest_id);

	QDPIO::cout << "Copying monomial: " << source_id << " to " << dest_id<< endl;
	dest_mon->setInternalFields( *src_mon );
      }
    }


}
