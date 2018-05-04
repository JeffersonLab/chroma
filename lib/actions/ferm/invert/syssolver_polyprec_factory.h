// -*- C++ -*-
/*! \file
 *  \brief Factory for solving PolyPrec*psi=chi where PolyPrec is hermitian and pos. def.
 */

#ifndef __syssolver_polyprec_factory_h__
#define __syssolver_polyprec_factory_h__

#include "singleton.h"
#include "objfactory.h"
#include "actions/ferm/invert/syssolver_polyprec.h"

namespace Chroma
{

  //! PolyPrec system solver factory (foundry)
  /*! @ingroup invert */
  typedef SingletonHolder< 
    ObjectFactory<PolyPrecSystemSolver<LatticeFermion>, 
		  std::string,
		  TYPELIST_3(XMLReader&, const std::string&, Handle< LinearOperator<LatticeFermion> >),
		  PolyPrecSystemSolver<LatticeFermion>* (*)(XMLReader&,
							    const std::string&,
							    Handle< LinearOperator<LatticeFermion> >), 
		  StringFactoryError> >
  ThePolyPrecFermSystemSolverFactory;

}


#endif
