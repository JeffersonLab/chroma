// -*- C++ -*-
// $Id: multi_syssolver_mdagm_factory.h,v 3.1 2006-07-03 15:26:08 edwards Exp $
/*! \file
 *  \brief Factory for producing system solvers for MdagM*psi = chi
 */

#ifndef __multi_syssolver_mdagm_factory_h__
#define __multi_syssolver_mdagm_factory_h__

#include "singleton.h"
#include "objfactory.h"
#include "linearop.h"
#include "state.h"

#include "actions/ferm/invert/multi_syssolver_mdagm.h"

namespace Chroma
{
  namespace { 
    typedef Handle< FermState< LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > FSHandle;
  }

  //! MdagM system solver factory (foundry)
  /*! @ingroup invert */
  typedef SingletonHolder< 
    ObjectFactory<MdagMMultiSystemSolver<LatticeFermion>, 
		  std::string,
		  TYPELIST_4(XMLReader&, const std::string&, FSHandle, Handle< LinearOperator<LatticeFermion> >),
		  MdagMMultiSystemSolver<LatticeFermion>* (*)(XMLReader&,
							      const std::string&,
							      FSHandle,
							      Handle< LinearOperator<LatticeFermion> >), 
		  StringFactoryError> >
  TheMdagMFermMultiSystemSolverFactory;


  //! MdagM system solver factory (foundry)
  /*! @ingroup invert */
  typedef SingletonHolder< 
    ObjectFactory<MdagMMultiSystemSolverArray<LatticeFermion>, 
		  std::string,
		  TYPELIST_4(XMLReader&, const std::string&, FSHandle, Handle< LinearOperatorArray<LatticeFermion> >),
		  MdagMMultiSystemSolverArray<LatticeFermion>* (*)(XMLReader&,
								   const std::string&,
								   FSHandle,
								   Handle< LinearOperatorArray<LatticeFermion> >), 
		  StringFactoryError> >
  TheMdagMFermMultiSystemSolverArrayFactory;


  //! MdagM system solver factory (foundry)
  /*! @ingroup invert */
  typedef SingletonHolder< 
    ObjectFactory<MdagMMultiSystemSolver<LatticeStaggeredFermion>, 
		  std::string,
		  TYPELIST_3(XMLReader&, const std::string&, Handle< LinearOperator<LatticeStaggeredFermion> >),
		  MdagMMultiSystemSolver<LatticeStaggeredFermion>* (*)(XMLReader&,
								       const std::string&,
								       Handle< LinearOperator<LatticeStaggeredFermion> >), 
		  StringFactoryError> >
  TheMdagMStagFermMultiSystemSolverFactory;

}


#endif
