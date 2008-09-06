// -*- C++ -*-
// $Id: multi_syssolver_mdagm_accumulate_factory.h,v 3.2 2008-09-06 18:35:35 bjoo Exp $
/*! \file
 *  \brief Factory for producing system solvers for MdagM*psi = chi
 */

#ifndef __multi_syssolver_mdagm_accumulate_factory_h__
#define __multi_syssolver_mdagm_accumulate_factory_h__

#include "singleton.h"
#include "objfactory.h"
#include "linearop.h"
#include "actions/ferm/invert/multi_syssolver_mdagm_accumulate.h"

namespace Chroma
{

  //! MdagM system solver factory (foundry)
  /*! @ingroup invert */
  typedef SingletonHolder< 
    ObjectFactory<MdagMMultiSystemSolverAccumulate<LatticeFermion>, 
		  std::string,
		  TYPELIST_3(XMLReader&, const std::string&, Handle< LinearOperator<LatticeFermion> >),
		  MdagMMultiSystemSolverAccumulate<LatticeFermion>* (*)(XMLReader&,
							      const std::string&,
							      Handle< LinearOperator<LatticeFermion> >), 
		  StringFactoryError> >
  TheMdagMFermMultiSystemSolverAccumulateFactory;


  //! MdagM system solver factory (foundry)
  /*! @ingroup invert */
  typedef SingletonHolder< 
    ObjectFactory<MdagMMultiSystemSolverAccumulateArray<LatticeFermion>, 
		  std::string,
		  TYPELIST_3(XMLReader&, const std::string&, Handle< LinearOperatorArray<LatticeFermion> >),
		  MdagMMultiSystemSolverAccumulateArray<LatticeFermion>* (*)(XMLReader&,
								   const std::string&,
								   Handle< LinearOperatorArray<LatticeFermion> >), 
		  StringFactoryError> >
  TheMdagMFermMultiSystemSolverAccumulateArrayFactory;


  //! MdagM system solver factory (foundry)
  /*! @ingroup invert */
  typedef SingletonHolder< 
    ObjectFactory<MdagMMultiSystemSolverAccumulate<LatticeStaggeredFermion>, 
		  std::string,
		  TYPELIST_3(XMLReader&, const std::string&, Handle< LinearOperator<LatticeStaggeredFermion> >),
		  MdagMMultiSystemSolverAccumulate<LatticeStaggeredFermion>* (*)(XMLReader&,
								       const std::string&,
								       Handle< LinearOperator<LatticeStaggeredFermion> >), 
		  StringFactoryError> >
  TheMdagMStagFermMultiSystemSolverAccumulateFactory;

}


#endif
