// $Id: simple_gaugebc.cc,v 3.2 2009-04-17 02:05:35 bjoo Exp $
/*! \file
 *  \brief Simple gauge boundary conditions
 */

#include "actions/gauge/gaugebcs/simple_gaugebc.h"
#include "actions/gauge/gaugebcs/gaugebc_factory.h"


namespace Chroma {

  namespace SimpleGaugeBCEnv 
  { 
    //! Calllback function to register with the factory
    GaugeBC< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >* createGaugeBC(XMLReader& xml, 
										       const string& path)
    {
      QDPIO::cout << "Factory Callback: Creating SimpleGaugeBC " << endl;
      return new SimpleGaugeBC<
        multi1d<LatticeColorMatrix>, 
	multi1d<LatticeColorMatrix> >(SimpleGaugeBCParams(xml, path));
    }

    GaugeBC< multi1d<LatticeColorMatrixF>, multi1d<LatticeColorMatrixF> >* createGaugeBCF(XMLReader& xml, 
										       const string& path)
    {
      QDPIO::cout << "Factory Callback: Creating SimpleGaugeBC " << endl;
      return new SimpleGaugeBC<
        multi1d<LatticeColorMatrixF>, 
	multi1d<LatticeColorMatrixF> >(SimpleGaugeBCParams(xml, path));
    }

    GaugeBC< multi1d<LatticeColorMatrixD>, multi1d<LatticeColorMatrixD> >* createGaugeBCD(XMLReader& xml, 
										       const string& path)
    {
      QDPIO::cout << "Factory Callback: Creating SimpleGaugeBC " << endl;
      return new SimpleGaugeBC<
        multi1d<LatticeColorMatrixD>, 
	multi1d<LatticeColorMatrixD> >(SimpleGaugeBCParams(xml, path));
    }

    const std::string name = "SIMPLE_GAUGEBC";

    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= TheGaugeBCFactory::Instance().registerObject(name, createGaugeBC);
	success &= TheGaugeBCFFactory::Instance().registerObject(name, createGaugeBCF);
	success &= TheGaugeBCDFactory::Instance().registerObject(name, createGaugeBCD);
	registered = true;
      }
      return success;
    }
  }

  SimpleGaugeBCParams::SimpleGaugeBCParams(XMLReader& xml, 
					   const std::string& path) {
      
    XMLReader paramtop(xml, path);
    
    try { 
      read(paramtop, "./boundary", boundary);
    }
    catch( const std::string& e ) { 
      QDPIO::cerr << "Error reading XML: " << e << endl;
      QDP_abort(1);
    }

    QDPIO::cout << "Creating SimpleGaugeBCParams with boundary: ";
    QDPIO::cout << "Boundary.size = " << boundary.size() << endl;
    for(int i = 0; i < boundary.size();i++) { 
      QDPIO::cout << boundary[i] << endl;
    }
    QDPIO::cout << endl;
  }




  void read(XMLReader& xml, const std::string& path, SimpleGaugeBCParams& p) {
    SimpleGaugeBCParams tmp(xml, path);
    p=tmp;
  }
  
  void write(XMLWriter& xml, const std::string& path, const SimpleGaugeBCParams& p) { 
    push(xml, path);
    write(xml, "boundary", p.boundary);
    pop(xml);
  }

} // End namespace Chroma 
