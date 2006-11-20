// -*- C++ -*-
// $Id: lcm_toplevel_integrator.h,v 3.1 2006-11-20 19:15:02 bjoo Exp $

#ifndef __LCM_TOPLEVEL_INTEGRATOR_H_
#define __LCM_TOPLEVEL_INTEGRATOR_H_

#include "chromabase.h"
#include <update/molecdyn/integrator/abs_integrator.h>

using namespace QDP;

namespace Chroma {

  //! A Struct for holding pairs of IDs in the copy list
  struct IDPair {
    std::string source;
    std::string dest;
  };
  
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

  //! A Structure to hold the top level parameters
  struct LCMToplevelIntegratorParams { 

    //! Construct from XML
    LCMToplevelIntegratorParams(XMLReader& xml, const std::string& path) {
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
      }
      catch(const std::string& e) { 
	QDPIO::cout << "Caught Exception Reading XML: " << e << endl;
	QDP_abort(1);
      }
    }

    Real tau0;
    multi1d<IDPair> copy_list;
    std::string integrator_xml;
  };

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
    pop(xml);
  }

  class LCMToplevelIntegrator : public 
  AbsMDIntegratorNew< multi1d<LatticeColorMatrix>,
		      multi1d<LatticeColorMatrix> >
  {
  public:

    //! Construct from parameters.
    LCMToplevelIntegrator(const LCMToplevelIntegratorParams& p) :
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
    }

    //! Copy Constructor
    LCMToplevelIntegrator(const LCMToplevelIntegrator& i) : params(i.params), top_integrator(i.top_integrator) {}

    //! Destructor is automagic
    ~LCMToplevelIntegrator() {} 

    //! Get the length of a trajectory
    Real getTrajLength(void) const { 
      return params.tau0;
    }

    //! Copy fields between the monomials of a copy list
    void copyFields(void) const { 
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

    AbsComponentIntegrator< multi1d<LatticeColorMatrix>, 
			    multi1d<LatticeColorMatrix> >& getIntegrator(void) const { 
      return *top_integrator;
    }

  private:
    LCMToplevelIntegratorParams params;
    Handle< AbsComponentIntegrator< multi1d<LatticeColorMatrix>, 
				    multi1d<LatticeColorMatrix> > > top_integrator;
  };


};


#endif
