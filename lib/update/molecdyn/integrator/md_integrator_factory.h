// -*- C++ -*-
// $Id: md_integrator_factory.h,v 1.4 2005-04-10 21:46:42 edwards Exp $
/*! \file
 *  \brief Integrator factories
 */

#ifndef __md_integrator_factory_h__
#define __md_integrator_factory_h__

#include "singleton.h"
#include "objfactory.h"
#include "chromabase.h"

#include "update/molecdyn/hamiltonian/abs_hamiltonian.h"
#include "update/molecdyn/integrator/abs_integrator.h"

#include <string>


namespace Chroma
{


  // Hack -- For some reason TYPELIST3 doesn't like having an
  // Abs Hamiltonian with 2 template params in Typelist. It seems
  // (to be confirmed) to think that the 2 template parameters mean
  // that we should be using typelist 4. This hacks around that using
  // a convenience typedef
  /*! @ingroup integrator */
  typedef AbsHamiltonian<multi1d<LatticeColorMatrix>,multi1d<LatticeColorMatrix> > LCMHam;
  
  //! A factory for exact non-fermionic monomials
  /*! @ingroup integrator */
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
