// $Id: schr_coupling_fermbc_w.cc,v 3.0 2006-04-03 04:58:48 edwards Exp $
/*! \file
 *  \brief Schroedinger BC - use for coupling determinations
 */

#include "actions/ferm/fermbcs/schr_coupling_fermbc_w.h"
#include "actions/ferm/fermbcs/fermbc_factory_w.h"

namespace Chroma 
{

  namespace SchrCouplingFermBCEnv 
  { 
    //! Callback function to register with the factory
    FermBC<LatticeFermion,
	   multi1d<LatticeColorMatrix>, 
	   multi1d<LatticeColorMatrix> >* createFermBC(XMLReader& xml, const string& path)
    {
      return new SchrCouplingFermBC(SchrCouplingGaugeBC(SchrGaugeBCParams(xml, path)),
				    SchrFermBCParams(xml, path));
    }

    const std::string name = "SCHROEDINGER_COUPLING_FERMBC";
    const bool registered = TheWilsonTypeFermBCFactory::Instance().registerObject(name,
										  createFermBC);
  }


}
