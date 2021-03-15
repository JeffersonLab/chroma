// -*- C++ -*-
/*! \file
 *  \brief Fermion create state factory
 */

#ifndef __ferm_createstate_factory_w_h__
#define __ferm_createstate_factory_w_h__

#include "singleton.h"
#include "objfactory.h"
#include "chromabase.h"

#include "create_state.h"

namespace Chroma
{

  //! CreateFermState Factory 
  /*! @ingroup fermstates */
  typedef Chroma::SingletonHolder< 
  ObjectFactory<CreateFermState<LatticeFermion,
				multi1d<LatticeColorMatrix>, 
				multi1d<LatticeColorMatrix> >, 
    std::string,
    TYPELIST_2(XMLReader&, const std::string&),
    CreateFermState<LatticeFermion,
		    multi1d<LatticeColorMatrix>, 
		    multi1d<LatticeColorMatrix> >* (*)(XMLReader&, 
						       const std::string&), 
		StringFactoryError> >
  TheCreateFermStateFactory;

} // end namespace Chroma


#endif
