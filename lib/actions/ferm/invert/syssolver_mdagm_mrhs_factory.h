/*
 * syssolver_mdagm_mrhs_factory.h
 *
 *  Created on: Mar 11, 2019
 *      Author: bjoo
 */

#ifndef LIB_ACTIONS_FERM_INVERT_SYSSOLVER_MDAGM_MRHS_FACTORY_H_
#define LIB_ACTIONS_FERM_INVERT_SYSSOLVER_MDAGM_MRHS_FACTORY_H_

#include "chromabase.h"
#include "state.h"
#include "singleton.h"
#include "objfactory.h"
#include "linearop.h"
#include "typelist.h"
#include "fermact.h"
#include "actions/ferm/invert/syssolver_mdagm_mrhs.h"

using namespace QDP;


namespace Chroma
{
  namespace {
  	  using FSHandle =  Handle< FermState< LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > >;
  	  using FA = FermAct4D< LatticeFermion, multi1d<LatticeColorMatrix>,multi1d<LatticeColorMatrix> >;
  }
  //! LinOp system solver factory (foundry)
  /*! @ingroup invert */
  typedef SingletonHolder<
    ObjectFactory<MdagMMRHSSystemSolver<LatticeFermion>,
		  std::string,
		  TYPELIST_4(XMLReader&,const std::string&,const FA&,FSHandle),
		  MdagMMRHSSystemSolver<LatticeFermion>* (*)(XMLReader&,const std::string&,const FA&,FSHandle),
		  StringFactoryError> >
  TheMdagMFermMRHSSystemSolverFactory;


} // Namespace chroma





#endif /* LIB_ACTIONS_FERM_INVERT_SYSSOLVER_MDAGM_MRHS_FACTORY_H_ */
