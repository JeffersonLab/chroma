// $Id: schr_chromomag_fermbc_w.cc,v 3.1 2006-09-20 20:28:00 edwards Exp $
/*! \file
 *  \brief Schroedinger BC - chromo-magnetic ferm BC
 */

#include "actions/ferm/fermbcs/schr_chromomag_fermbc_w.h"
#include "actions/ferm/fermbcs/fermbc_factory_w.h"

namespace Chroma 
{

  namespace SchrChromoMagFermBCEnv 
  { 
    //! Callback function to register with the factory
    FermBC<LatticeFermion,
	   multi1d<LatticeColorMatrix>, 
	   multi1d<LatticeColorMatrix> >* createFermBC(XMLReader& xml, const string& path)
    {
      return new SchrChromoMagFermBC(SchrChromoMagGaugeBC(SchrGaugeBCParams(xml, path)),
				     SchrFermBCParams(xml, path));
    }

    const std::string name = "SCHROEDINGER_CHROMOMAG_FERMBC";

    static bool registered = false;

    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &=  TheWilsonTypeFermBCFactory::Instance().registerObject(name, createFermBC);
	registered = true;
      }
      return success;
    }
  }

}
