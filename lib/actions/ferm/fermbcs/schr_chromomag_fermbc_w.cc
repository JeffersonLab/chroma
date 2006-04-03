// $Id: schr_chromomag_fermbc_w.cc,v 3.0 2006-04-03 04:58:48 edwards Exp $
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
    const bool registered = TheWilsonTypeFermBCFactory::Instance().registerObject(name,
										  createFermBC);
  }

}
