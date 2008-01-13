// -*- C++ -*-
// $Id: norm_gram_schm.h,v 1.2 2008-01-13 21:08:14 edwards Exp $
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
  // 4D versions
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

  // 5D versions
  void normGramSchmidt(multi2d<LatticeFermionF>& vec, 
		       int f,
		       int t,
		       const Subset& sub);

  void normGramSchmidt(multi2d<LatticeFermionD>& vec, 
		       int f,
		       int t,
		       const Subset& sub);

  /*! @} */  // end of group invert

  
}// End Namespace Chroma

#endif 
