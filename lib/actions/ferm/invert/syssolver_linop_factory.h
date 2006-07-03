// -*- C++ -*-
// $Id: syssolver_linop_factory.h,v 3.1 2006-07-03 15:26:08 edwards Exp $
/*! \file
 *  \brief Factory for solving M*psi=chi where M is not hermitian or pos. def.
 */

#ifndef __syssolver_linop_factory_h__
#define __syssolver_linop_factory_h__

#include "singleton.h"
#include "objfactory.h"
#include "actions/ferm/invert/syssolver_linop.h"

namespace Chroma
{

  //! LinOp system solver factory (foundry)
  /*! @ingroup invert */
  typedef SingletonHolder< 
    ObjectFactory<LinOpSystemSolver<LatticeFermion>, 
		  std::string,
		  TYPELIST_3(XMLReader&, const std::string&, Handle< LinearOperator<LatticeFermion> >),
		  LinOpSystemSolver<LatticeFermion>* (*)(XMLReader&,
							 const std::string&,
							 Handle< LinearOperator<LatticeFermion> >), 
		  StringFactoryError> >
  TheLinOpFermSystemSolverFactory;


  //! LinOp system solver factory (foundry)
  /*! @ingroup invert */
  typedef SingletonHolder< 
    ObjectFactory<LinOpSystemSolverArray<LatticeFermion>, 
		  std::string,
		  TYPELIST_3(XMLReader&, const std::string&, Handle< LinearOperatorArray<LatticeFermion> >),
		  LinOpSystemSolverArray<LatticeFermion>* (*)(XMLReader&,
							      const std::string&,
							      Handle< LinearOperatorArray<LatticeFermion> >), 
		  StringFactoryError> >
  TheLinOpFermSystemSolverArrayFactory;



  //! LinOp system solver factory (foundry)
  /*! @ingroup invert */
  typedef SingletonHolder< 
    ObjectFactory<LinOpSystemSolver<LatticeStaggeredFermion>, 
		  std::string,
		  TYPELIST_3(XMLReader&, const std::string&, Handle< LinearOperator<LatticeStaggeredFermion> >),
		  LinOpSystemSolver<LatticeStaggeredFermion>* (*)(XMLReader&,
								  const std::string&,
								  Handle< LinearOperator<LatticeStaggeredFermion> >), 
		  StringFactoryError> >
  TheLinOpStagFermSystemSolverFactory;


}


#endif
