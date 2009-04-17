// -*- C++ -*-
// $Id: syssolver_mdagm_factory.h,v 3.2 2009-04-17 02:05:32 bjoo Exp $
/*! \file
 *  \brief Factory for producing system solvers for MdagM*psi = chi
 */

#ifndef __syssolver_mdagm_factory_h__
#define __syssolver_mdagm_factory_h__

#include "chromabase.h"
#include "state.h"
#include "singleton.h"
#include "objfactory.h"
#include "linearop.h"
#include "typelist.h"
#include "actions/ferm/invert/syssolver_mdagm.h"

using namespace QDP;

namespace Chroma
{

  namespace FactoryEnv { 
    typedef Handle< FermState< LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > FSHandle;
    typedef Handle< FermState< LatticeFermionF, multi1d<LatticeColorMatrixF>, multi1d<LatticeColorMatrixF> > > FSHandleF;
    typedef Handle< FermState< LatticeFermionD, multi1d<LatticeColorMatrixD>, multi1d<LatticeColorMatrixD> > > FSHandleD;
  }

  //! MdagM system solver factory (foundry)
  /*! @ingroup invert */
  typedef SingletonHolder< 
    ObjectFactory<MdagMSystemSolver<LatticeFermion>, 
		  std::string,
		  TYPELIST_4(XMLReader&, const std::string&, FactoryEnv::FSHandle, Handle< LinearOperator<LatticeFermion> >),
		  MdagMSystemSolver<LatticeFermion>* (*)(XMLReader&,
							 const std::string&,
							 FactoryEnv::FSHandle,
							 Handle< LinearOperator<LatticeFermion> >), 
		  StringFactoryError> >
  TheMdagMFermSystemSolverFactory;

  typedef SingletonHolder< 
    ObjectFactory<MdagMSystemSolver<LatticeFermionF>, 
		  std::string,
		  TYPELIST_4(XMLReader&, const std::string&, FactoryEnv::FSHandleF, Handle< LinearOperator<LatticeFermionF > >),
		  MdagMSystemSolver<LatticeFermionF>* (*)(XMLReader&,
							  const std::string&,
							  FactoryEnv::FSHandleF, 
							  Handle< LinearOperator<LatticeFermionF> >), 
		  StringFactoryError> >
  TheMdagMFermFSystemSolverFactory;

  typedef SingletonHolder< 
    ObjectFactory<MdagMSystemSolver<LatticeFermionD>, 
		  std::string,
		  TYPELIST_4(XMLReader&, const std::string&, FactoryEnv::FSHandleD, Handle< LinearOperator<LatticeFermionD> >),
		  MdagMSystemSolver<LatticeFermionD>* (*)(XMLReader&,
							  const std::string&,
							  FactoryEnv::FSHandleD,
							  Handle< LinearOperator<LatticeFermionD> >), 
		  StringFactoryError> >
  TheMdagMFermDSystemSolverFactory;


  //! MdagM system solver factory (foundry)
  /*! @ingroup invert */
  typedef SingletonHolder< 
    ObjectFactory<MdagMSystemSolverArray<LatticeFermion>, 
		  std::string,
		  TYPELIST_3(XMLReader&, const std::string&, Handle< LinearOperatorArray<LatticeFermion> >),
		  MdagMSystemSolverArray<LatticeFermion>* (*)(XMLReader&,
							      const std::string&,
							      Handle< LinearOperatorArray<LatticeFermion> >), 
		  StringFactoryError> >
  TheMdagMFermSystemSolverArrayFactory;


  //! MdagM system solver factory (foundry)
  /*! @ingroup invert */
  typedef SingletonHolder< 
    ObjectFactory<MdagMSystemSolver<LatticeStaggeredFermion>, 
		  std::string,
		  TYPELIST_3(XMLReader&, const std::string&, Handle< LinearOperator<LatticeStaggeredFermion> >),
		  MdagMSystemSolver<LatticeStaggeredFermion>* (*)(XMLReader&,
								  const std::string&,
								  Handle< LinearOperator<LatticeStaggeredFermion> >), 
		  StringFactoryError> >
  TheMdagMStagFermSystemSolverFactory;

}


#endif
