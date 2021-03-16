// -*- C++ -*-
/*! \file
 *  \brief Fermion Boundary Condition factories
 */

#ifndef __fermbc_factory_w_h__
#define __fermbc_factory_w_h__

#include "qdp_precision.h"

#include "singleton.h"
#include "objfactory.h"

#include "fermbc.h"


namespace Chroma
{
  //! FermBC factory
  /*! \ingroup fermbcs */
  typedef Chroma::SingletonHolder< 
    ObjectFactory<FermBC<LatticeFermion, 
			 multi1d<LatticeColorMatrix>, 
			 multi1d<LatticeColorMatrix> >, 
		  std::string,
		  TYPELIST_2(XMLReader&, const std::string&),
		  FermBC<LatticeFermion,
			 multi1d<LatticeColorMatrix>, 
			 multi1d<LatticeColorMatrix> >* (*)(XMLReader&, 
							    const std::string&), 
		  StringFactoryError> >
  TheWilsonTypeFermBCFactory;

  typedef Chroma::SingletonHolder< 
    ObjectFactory<FermBC<LatticeFermionF, 
			 multi1d<LatticeColorMatrixF>, 
			 multi1d<LatticeColorMatrixF> >, 
		  std::string,
		  TYPELIST_2(XMLReader&, const std::string&),
		  FermBC<LatticeFermionF,
			 multi1d<LatticeColorMatrixF>, 
			 multi1d<LatticeColorMatrixF> >* (*)(XMLReader&, 
							    const std::string&), 
		  StringFactoryError> >
  TheWilsonTypeFermBCFFactory;

  typedef Chroma::SingletonHolder< 
    ObjectFactory<FermBC<LatticeFermionD, 
			 multi1d<LatticeColorMatrixD>, 
			 multi1d<LatticeColorMatrixD> >, 
		  std::string,
		  TYPELIST_2(XMLReader&, const std::string&),
		  FermBC<LatticeFermionD,
			 multi1d<LatticeColorMatrixD>, 
			 multi1d<LatticeColorMatrixD> >* (*)(XMLReader&, 
							    const std::string&), 
		  StringFactoryError> >
  TheWilsonTypeFermBCDFactory;


}


#endif
