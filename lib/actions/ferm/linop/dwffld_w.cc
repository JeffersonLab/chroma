// $Id: dwffld_w.cc,v 3.1 2007-09-01 23:44:10 uid3790 Exp $
/*! \file
 *  \brief DWF parity/rotation operator
 *
 * Apply a DWF style `parity' or rotation operator to move from 
 * boundary field basis to single field basis or back again 
 */

#include "chromabase.h"
#include "actions/ferm/linop/dwffld_w.h"


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
    chi.resize(N5);

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


} // End Namespace Chroma
