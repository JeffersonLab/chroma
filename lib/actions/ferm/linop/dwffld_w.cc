// $Id: dwffld_w.cc,v 1.2 2003-11-12 22:11:53 edwards Exp $
/*! \file
 *  \brief DWF parity/rotation operator
 *
 * Apply a DWF style `parity' or rotation operator to move from 
 * boundary field basis to single field basis or back again 
 */

#include "chromabase.h"
#include "actions/ferm/linop/dwffld_w.h"

using namespace QDP;

// HACK
/* Terrible implementation just to get the silly thing going */
static inline LatticeFermion 
chiralProjectPlus(const LatticeFermion& l)
{
  return 0.5*(l + Gamma(15)*l);
}

static inline LatticeFermion 
chiralProjectMinus(const LatticeFermion& l)
{
  return 0.5*(l - Gamma(15)*l);
}


//! DWF parity/rotation operator
/*! \ingroup linop
 *
 *  Chi  :=  P^{isign} . Psi    where  P is the rotation operator 
 *
 *  \param psi        Pseudofermion field                     (Read)
 *  \param chi        Pseudofermion field                     (Write)
 *  \param isign      Sign (Plus/Minus)    	              (Read)
 */

void DwfFld(LatticeDWFermion& chi, const LatticeDWFermion& psi, enum PlusMinus isign)
{
  START_CODE("DwfFld");
    
  switch (isign)
  {
  case PLUS:
    for(int i=0; i < Ls-1; ++i)
      pokeDW(chi, chiralProjectMinus(peekDW(psi,i)) + chiralProjectPlus(peekDW(psi,i+1)), i);

    pokeDW(chi, chiralProjectMinus(peekDW(psi,Ls-1)) + chiralProjectPlus(peekDW(psi,0)), Ls-1);
  break;

  case MINUS:
    pokeDW(chi, chiralProjectMinus(peekDW(psi,0)) + chiralProjectPlus(peekDW(psi,Ls-1)), 0);
    for(int i=1; i < Ls; ++i)
      pokeDW(chi, chiralProjectMinus(peekDW(psi,i)) + chiralProjectPlus(peekDW(psi,i-1)), i);
    break;

  default:
    QDP_error_exit("invalid option", isign);
  }

  END_CODE("DwfFld");
}


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
  START_CODE("DwfFld");
    
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

  END_CODE("DwfFld");
}



