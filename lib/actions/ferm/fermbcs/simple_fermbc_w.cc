// $Id: simple_fermbc_w.cc,v 3.2 2009-04-17 02:05:31 bjoo Exp $
/*! \file
 *  \brief Simple fermionic BC
 */

#include "actions/ferm/fermbcs/simple_fermbc_w.h"
#include "actions/ferm/fermbcs/fermbc_factory_w.h"

namespace Chroma
{

  //! Name and registration
  namespace WilsonTypeSimpleFermBCEnv
  {
    //! Callback function
    // Floating
    FermBC<LatticeFermion,
	   multi1d<LatticeColorMatrix>, 
	   multi1d<LatticeColorMatrix> >* createFermBC(XMLReader& xml_in, 
						       const std::string& path)
    {
      return new SimpleFermBC<LatticeFermion,
	                      multi1d<LatticeColorMatrix>, 
	                      multi1d<LatticeColorMatrix> >(SimpleFermBCParams(xml_in, path));
    }

    FermBC<LatticeFermionF,
	   multi1d<LatticeColorMatrixF>, 
	   multi1d<LatticeColorMatrixF> >* createFermBCF(XMLReader& xml_in, 
						       const std::string& path)
    {
      return new SimpleFermBC<LatticeFermionF,
	                      multi1d<LatticeColorMatrixF>, 
	                      multi1d<LatticeColorMatrixF> >(SimpleFermBCParams(xml_in, path));
    }

    FermBC<LatticeFermionD,
	   multi1d<LatticeColorMatrixD>, 
	   multi1d<LatticeColorMatrixD> >* createFermBCD(XMLReader& xml_in, 
						       const std::string& path)
    {
      return new SimpleFermBC<LatticeFermionD,
	                      multi1d<LatticeColorMatrixD>, 
	                      multi1d<LatticeColorMatrixD> >(SimpleFermBCParams(xml_in, path));
    }


    //! Name to be used
    const std::string name = "SIMPLE_FERMBC";

    static bool registered = false;

    //! Register all objects
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::TheWilsonTypeFermBCFactory::Instance().registerObject(name, createFermBC);

	success &= Chroma::TheWilsonTypeFermBCFFactory::Instance().registerObject(name, createFermBCF);

	success &= Chroma::TheWilsonTypeFermBCDFactory::Instance().registerObject(name, createFermBCD);

	registered = true;
      }
      return success;
    }
  }

}
