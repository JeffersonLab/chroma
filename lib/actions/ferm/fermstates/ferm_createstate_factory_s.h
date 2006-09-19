// -*- C++ -*-
// $Id: ferm_createstate_factory_s.h,v 1.1 2006-09-19 17:53:36 edwards Exp $
/*! \file
 *  \brief Fermion create state factory
 */

#ifndef __ferm_createstate_factory_s_h__
#define __ferm_createstate_factory_s_h__

#include "singleton.h"
#include "objfactory.h"
#include "chromabase.h"

#include "create_state.h"

namespace Chroma
{

  //! CreateFermState Factory 
  /*! @ingroup fermstates */
  typedef SingletonHolder< 
  ObjectFactory<CreateFermState<LatticeStaggeredFermion,
				multi1d<LatticeColorMatrix>, 
				multi1d<LatticeColorMatrix> >, 
    std::string,
    TYPELIST_2(XMLReader&, const std::string&),
    CreateFermState<LatticeStaggeredFermion,
		    multi1d<LatticeColorMatrix>, 
		    multi1d<LatticeColorMatrix> >* (*)(XMLReader&, 
						       const std::string&), 
		StringFactoryError> >
  TheStaggeredCreateFermStateFactory;

} // end namespace Chroma


#endif
