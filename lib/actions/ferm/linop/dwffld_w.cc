// $Id: dwffld_w.cc,v 1.6 2004-12-12 21:22:15 edwards Exp $
/*! \file
 *  \brief DWF parity/rotation operator
 *
 * Apply a DWF style `parity' or rotation operator to move from 
 * boundary field basis to single field basis or back again 
 */

#include "chromabase.h"
#include "actions/ferm/linop/dwffld_w.h"

using namespace QDP;

namespace Chroma 
{ 
//! DWF parity/rotation operator
/*! \ingroup linop
 *
 *  Chi  :=  P^{isign} . Psi    where  P is the rotation operator 
 *
 *  \param psi        Pseudofermion field                     (Read)
 *  \param chi        Pseudofermion field                     (Write)
 *  \param isign      Sign (Plus/Minus)    	              (Read)
 */

void DwfFld(multi1d<LatticeFermion>& chi, const multi1d<LatticeFermion>& psi, enum PlusMinus isign)
{
  START_CODE();
    
  const int N5 = psi.size();

  switch (isign)
  {
  case PLUS:
    for(int n=0; n < N5-1; ++n)
      chi[n] = chiralProjectMinus(psi[n]) + chiralProjectPlus(psi[n+1]);

    chi[N5-1] = chiralProjectMinus(psi[N5-1]) + chiralProjectPlus(psi[0]);
  break;

  case MINUS:
    chi[0] = chiralProjectMinus(psi[0]) + chiralProjectPlus(psi[N5-1]);
    for(int n=1; n < N5; ++n)
      chi[n] = chiralProjectMinus(psi[n]) + chiralProjectPlus(psi[n-1]);
    break;

  default:
    QDP_error_exit("invalid option", isign);
  }

  END_CODE();
}


}; // End Namespace Chroma



