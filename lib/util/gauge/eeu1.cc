// $Id: eeu1.cc,v 3.1 2007-08-29 13:15:37 edwards Exp $
/*! \file
 *  \brief Exactly exponentiate a U(1) lie algebra element
 */

#include "chromabase.h"
#include "util/gauge/eeu1.h"

namespace Chroma 
{
  
  //! Exactly exponentiate a U(1) lie algebra element
  /*!
   * \ingroup gauge
   *
   *  \param lambda      LatticeColorMatrix          (Modify)
   */
  
  void eeu1(LatticeColorMatrix& lambda)
  {
    START_CODE();
    
    if ( Nc != 1 )
    {
      QDPIO::cerr << __func__ << ": can only handle U(1)" << endl;
      QDP_abort(1);
    }

    LatticeReal phi = imag(trace(lambda));	   
    LatticeComplex a = cmplx(cos(phi),sin(phi));
    pokeColor(lambda, a, 0,0);
    
    END_CODE();
  }

}  // end namespace Chroma

