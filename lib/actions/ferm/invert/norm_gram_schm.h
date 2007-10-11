// -*- C++ -*-
// $Id: norm_gram_schm.h,v 1.1 2007-10-11 19:00:09 edwards Exp $
/*! \file
 *  \brief Gram-Schmidt with normalization
 */

#ifndef __norm_gram_schmidt_h__
#define __norm_gram_schmidt_h__

#include "chromabase.h"

namespace Chroma 
{

  //! Gram-Schmidt with normalization
  /*! \ingroup invert
   * @{
   */
  void normGramSchmidt(multi1d<LatticeFermionF>& vec, 
		       int f,
		       int t,
		       const Subset& sub);

  void normGramSchmidt(multi1d<LatticeFermionD>& vec, 
		       int f,
		       int t,
		       const Subset& sub);

  void normGramSchmidt(multi1d<LatticeStaggeredFermionF>& vec, 
		       int f,
		       int t,
		       const Subset& sub);

  void normGramSchmidt(multi1d<LatticeStaggeredFermionD>& vec, 
		       int f,
		       int t,
		       const Subset& sub);

  /*! @} */  // end of group invert

  
}// End Namespace Chroma

#endif 
