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
#include "actions/ferm/invert/syssolver_linop_mrhs.h"

namespace Chroma
{
  namespace {

    using FSHandle =  Handle< FermState< LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > >;

    using FSHandleF = Handle< FermState< LatticeFermionF, multi1d<LatticeColorMatrixF>, multi1d<LatticeColorMatrixF> > >;

    using FSHandleD = Handle< FermState< LatticeFermionD, multi1d<LatticeColorMatrixD>, multi1d<LatticeColorMatrixD> > >;
  }
  //! LinOp system solver factory (foundry)
  /*! @ingroup invert */
  typedef SingletonHolder<
    ObjectFactory<LinOpMRHSSystemSolver<LatticeFermion>,
		  std::string,
		  TYPELIST_4(XMLReader&, const std::string&, FSHandle,  Handle< LinearOperatorArray<LatticeFermion> >),
		  LinOpMRHSSystemSolver<LatticeFermion>* (*)(XMLReader&,
							 const std::string&,
							 FSHandle,
							 Handle< LinearOperatorArray<LatticeFermion> >),
		  StringFactoryError> >
  TheLinOpFermMRHSSystemSolverFactory;

  typedef SingletonHolder<
    ObjectFactory<LinOpMRHSSystemSolver<LatticeFermionF>,
		  std::string,
		  TYPELIST_4(XMLReader&, const std::string&, FSHandleF,  Handle< LinearOperatorArray<LatticeFermionF> >),
		  LinOpMRHSSystemSolver<LatticeFermionF>* (*)(XMLReader&,
							 const std::string&,
							 FSHandleF,
							 Handle< LinearOperatorArray<LatticeFermionF> >),
		  StringFactoryError> >
  TheLinOpFFermMRHSSystemSolverFactory;

  typedef SingletonHolder<
    ObjectFactory<LinOpMRHSSystemSolver<LatticeFermionD>,
		  std::string,
		  TYPELIST_4(XMLReader&, const std::string&, FSHandleD,  Handle< LinearOperatorArray<LatticeFermionD> >),
		  LinOpMRHSSystemSolver<LatticeFermionD>* (*)(XMLReader&,
							 const std::string&,
							 FSHandleD,
							 Handle< LinearOperatorArray<LatticeFermionD> >),
		  StringFactoryError> >
  TheLinOpDFermMRHSSystemSolverFactory;


} // Namespace chroma




#endif /* LIB_ACTIONS_FERM_INVERT_SYSSOLVER_LINOP_MRHS_FACTORY_H_ */
