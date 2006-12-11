// -*- C++ -*-
// $Id: baryon_spinmat_funcmap_w.h,v 3.1 2006-12-11 17:20:34 edwards Exp $
/*! \file
 *  \brief Factory for producing baryon spin matrix contraction objects
 */

#ifndef __baryon_spinmat_funcmap_w_h__
#define __baryon_spinmat_funcmap_w_h__

#include "singleton.h"
#include "funcmap.h"
#include "chromabase.h"

namespace Chroma
{
  //! Registration aggregator
  /*! \ingroup hadron */
  namespace BaryonSpinMatrixEnv
  {
    struct DumbDisambiguator {};

    //! Spin matrix factory
    /*! @ingroup hadron */
    typedef SingletonHolder< 
      FunctionMap<DumbDisambiguator,
		  SpinMatrix, 
		  std::string,
		  TYPELIST_2(XMLReader&, const std::string&),
		  SpinMatrix (*)(XMLReader&,
				 const std::string&), StringFunctionMapError> >
    TheBarSpinMatFuncMap;

    bool registerAll();
  }
}


#endif
