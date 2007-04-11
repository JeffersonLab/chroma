// -*- C++ -*-
// $Id: multi_syssolver_linop_factory.h,v 1.1 2007-04-11 03:41:36 edwards Exp $
/*! \file
 *  \brief Factory for producing system solvers for M*psi = chi
 */

#ifndef __multi_syssolver_linop_factory_h__
#define __multi_syssolver_linop_factory_h__

#include "singleton.h"
#include "objfactory.h"
#include "linearop.h"
#include "actions/ferm/invert/multi_syssolver_linop.h"

namespace Chroma
{

  //! LinOp system solver factory (foundry)
  /*! @ingroup invert */
  typedef SingletonHolder< 
    ObjectFactory<LinOpMultiSystemSolver<LatticeFermion>, 
		  std::string,
		  TYPELIST_3(XMLReader&, const std::string&, Handle< LinearOperator<LatticeFermion> >),
		  LinOpMultiSystemSolver<LatticeFermion>* (*)(XMLReader&,
							      const std::string&,
							      Handle< LinearOperator<LatticeFermion> >), 
		  StringFactoryError> >
  TheLinOpFermMultiSystemSolverFactory;


#if 0
  //! LinOp system solver factory (foundry)
  /*! @ingroup invert */
  typedef SingletonHolder< 
    ObjectFactory<LinOpMultiSystemSolverArray<LatticeFermion>, 
		  std::string,
		  TYPELIST_3(XMLReader&, const std::string&, Handle< LinearOperatorArray<LatticeFermion> >),
		  LinOpMultiSystemSolverArray<LatticeFermion>* (*)(XMLReader&,
								   const std::string&,
								   Handle< LinearOperatorArray<LatticeFermion> >), 
		  StringFactoryError> >
  TheLinOpFermMultiSystemSolverArrayFactory;
#endif


  //! LinOp system solver factory (foundry)
  /*! @ingroup invert */
  typedef SingletonHolder< 
    ObjectFactory<LinOpMultiSystemSolver<LatticeStaggeredFermion>, 
		  std::string,
		  TYPELIST_3(XMLReader&, const std::string&, Handle< LinearOperator<LatticeStaggeredFermion> >),
		  LinOpMultiSystemSolver<LatticeStaggeredFermion>* (*)(XMLReader&,
								       const std::string&,
								       Handle< LinearOperator<LatticeStaggeredFermion> >), 
		  StringFactoryError> >
  TheLinOpStagFermMultiSystemSolverFactory;

}


#endif
