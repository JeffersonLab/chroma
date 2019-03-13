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
#include "actions/ferm/invert/syssolver_mdagm_mrhs.h"

using namespace QDP;

namespace Chroma
{

  namespace FactoryEnv {
    using FSHandle  = Handle< FermState< LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > >;
    using FSHandleF = Handle< FermState< LatticeFermionF, multi1d<LatticeColorMatrixF>, multi1d<LatticeColorMatrixF> > >;
    using FSHandleD = Handle< FermState< LatticeFermionD, multi1d<LatticeColorMatrixD>, multi1d<LatticeColorMatrixD> > >;
  }

  //! MdagM system solver factory (foundry)
  /*! @ingroup invert */
  typedef SingletonHolder<
    ObjectFactory<MdagMMRHSSystemSolver<LatticeFermion>,
		  std::string,
		  TYPELIST_4(XMLReader&, const std::string&, FactoryEnv::FSHandle, Handle< LinearOperatorArray<LatticeFermion> >),
		  MdagMMRHSSystemSolver<LatticeFermion>* (*)(XMLReader&,
							 const std::string&,
							 FactoryEnv::FSHandle,
							 Handle< LinearOperatorArray<LatticeFermion> >),
		  StringFactoryError> >
  TheMdagMFermMRHSSystemSolverFactory;

  typedef SingletonHolder<
    ObjectFactory<MdagMMRHSSystemSolver<LatticeFermionF>,
		  std::string,
		  TYPELIST_4(XMLReader&, const std::string&, FactoryEnv::FSHandleF, Handle< LinearOperatorArray<LatticeFermionF > >),
		  MdagMMRHSSystemSolver<LatticeFermionF>* (*)(XMLReader&,
							  const std::string&,
							  FactoryEnv::FSHandleF,
							  Handle< LinearOperatorArray<LatticeFermionF> >),
		  StringFactoryError> >
  TheMdagMFFermMRHSSystemSolverFactory;

  typedef SingletonHolder<
    ObjectFactory<MdagMMRHSSystemSolver<LatticeFermionD>,
		  std::string,
		  TYPELIST_4(XMLReader&, const std::string&, FactoryEnv::FSHandleD, Handle< LinearOperatorArray<LatticeFermionD> >),
		  MdagMMRHSSystemSolver<LatticeFermionD>* (*)(XMLReader&,
							  const std::string&,
							  FactoryEnv::FSHandleD,
							  Handle< LinearOperatorArray<LatticeFermionD> >),
		  StringFactoryError> >
  TheMdagMDFermMRHSSystemSolverFactory;

} // namespace



#endif /* LIB_ACTIONS_FERM_INVERT_SYSSOLVER_MDAGM_MRHS_FACTORY_H_ */
