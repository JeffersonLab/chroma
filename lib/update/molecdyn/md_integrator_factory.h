// -*- C++ -*-
// $Id: md_integrator_factory.h,v 1.1 2004-12-31 17:17:04 bjoo Exp $
/*! \file
 *  \brief Monomial factories
 */

#ifndef __md_integrator_factory_h__
#define __md_integrator_factory_h__

#include "singleton.h"
#include "objfactory.h"
#include "chromabase.h"

#include "update/molecdyn/abs_hamiltonian.h"
#include "update/molecdyn/abs_integrator.h"

#include <string>

using namespace QDP;
using namespace Chroma;
using namespace std;

namespace Chroma
{


  // Hack -- For some reason TYPELIST3 doesn't like having an
  // Abs Hamiltonian with 2 template params in Typelist. It seems
  // (to be confirmed) to think that the 2 template parameters mean
  // that we should be using typelist 4. This hacks around that using
  // a convenience typedef
  typedef AbsHamiltonian<multi1d<LatticeColorMatrix>,multi1d<LatticeColorMatrix> > LCMHam;
  
  //! A factory for exact non-fermionic monomials
  typedef SingletonHolder< 
  ObjectFactory<
    AbsMDIntegrator< multi1d<LatticeColorMatrix>, 
		     multi1d<LatticeColorMatrix> >,
    std::string,

    TYPELIST_3(XMLReader&, 
	       const std::string&, 
	       Handle< LCMHam >& ),

    AbsMDIntegrator< multi1d<LatticeColorMatrix>, 
		     multi1d<LatticeColorMatrix> >* (*)(
							XMLReader&,
						        const std::string&,
							Handle<LCMHam>&),
    StringFactoryError> >
  TheMDIntegratorFactory;


}; // End namespace Chroma


#endif
