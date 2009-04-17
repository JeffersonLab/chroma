// -*- C++ -*-
// $Id: syssolver_linop_factory.h,v 3.2 2009-04-17 02:05:31 bjoo Exp $
/*! \file
 *  \brief Factory for solving M*psi=chi where M is not hermitian or pos. def.
 */

#ifndef __syssolver_linop_factory_h__
#define __syssolver_linop_factory_h__

#include "chromabase.h"
#include "handle.h"
#include "state.h"
#include "singleton.h"
#include "typelist.h"
#include "objfactory.h"
#include "actions/ferm/invert/syssolver_linop.h"

namespace Chroma
{
  namespace { 

    typedef Handle< FermState< LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > FSHandle;

    typedef Handle< FermState< LatticeFermionF, multi1d<LatticeColorMatrixF>, multi1d<LatticeColorMatrixF> > > FSHandleF;

    typedef Handle< FermState< LatticeFermionD, multi1d<LatticeColorMatrixD>, multi1d<LatticeColorMatrixD> > > FSHandleD;
  }
  //! LinOp system solver factory (foundry)
  /*! @ingroup invert */
  typedef SingletonHolder< 
    ObjectFactory<LinOpSystemSolver<LatticeFermion>, 
		  std::string,
		  TYPELIST_4(XMLReader&, const std::string&, FSHandle,  Handle< LinearOperator<LatticeFermion> >),
		  LinOpSystemSolver<LatticeFermion>* (*)(XMLReader&,
							 const std::string&,

							 FSHandle,
							 Handle< LinearOperator<LatticeFermion> >), 
		  StringFactoryError> >
  TheLinOpFermSystemSolverFactory;

  typedef SingletonHolder< 
    ObjectFactory<LinOpSystemSolver<LatticeFermionF>, 
		  std::string,
		  TYPELIST_4(XMLReader&, const std::string&, FSHandleF,  Handle< LinearOperator<LatticeFermionF> >),
		  LinOpSystemSolver<LatticeFermionF>* (*)(XMLReader&,
							 const std::string&,

							 FSHandleF,
							 Handle< LinearOperator<LatticeFermionF> >), 
		  StringFactoryError> >
  TheLinOpFFermSystemSolverFactory;

  typedef SingletonHolder< 
    ObjectFactory<LinOpSystemSolver<LatticeFermionD>, 
		  std::string,
		  TYPELIST_4(XMLReader&, const std::string&, FSHandleD,  Handle< LinearOperator<LatticeFermionD> >),
		  LinOpSystemSolver<LatticeFermionD>* (*)(XMLReader&,
							 const std::string&,

							 FSHandleD,
							 Handle< LinearOperator<LatticeFermionD> >), 
		  StringFactoryError> >
  TheLinOpDFermSystemSolverFactory;


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
