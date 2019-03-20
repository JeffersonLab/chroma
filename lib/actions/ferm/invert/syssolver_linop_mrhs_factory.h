/*
 * syssolver_linop_mrhs_factory.h
 *
 *  Created on: Mar 11, 2019
 *      Author: bjoo
 */

#ifndef LIB_ACTIONS_FERM_INVERT_SYSSOLVER_LINOP_MRHS_FACTORY_H_
#define LIB_ACTIONS_FERM_INVERT_SYSSOLVER_LINOP_MRHS_FACTORY_H_

#include "chromabase.h"
#include "handle.h"
#include "state.h"
#include "singleton.h"
#include "typelist.h"
#include "objfactory.h"
#include "fermact.h"
#include "actions/ferm/invert/syssolver_linop_mrhs.h"

namespace Chroma
{
  namespace {
  	  using FSHandle =  Handle< FermState< LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > >;
  	  using FA = FermAct4D< LatticeFermion, multi1d<LatticeColorMatrix>,multi1d<LatticeColorMatrix> >;
  }
  //! LinOp system solver factory (foundry)
  /*! @ingroup invert */
  typedef SingletonHolder<
    ObjectFactory<LinOpMRHSSystemSolver<LatticeFermion>,
		  std::string,
		  TYPELIST_4(XMLReader&,const std::string&,const FA&,FSHandle),
		  LinOpMRHSSystemSolver<LatticeFermion>* (*)(XMLReader&,const std::string&,const FA&,FSHandle),
		  StringFactoryError> >
  TheLinOpFermMRHSSystemSolverFactory;


} // Namespace chroma




#endif /* LIB_ACTIONS_FERM_INVERT_SYSSOLVER_LINOP_MRHS_FACTORY_H_ */
