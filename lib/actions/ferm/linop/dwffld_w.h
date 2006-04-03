// -*- C++ -*-
// $Id: dwffld_w.h,v 3.0 2006-04-03 04:58:50 edwards Exp $
/*! \file
 *  \brief DWF parity/rotation operator
 *
 * Apply a DWF style `parity' or rotation operator to move from 
 * boundary field basis to single field basis or back again 
 */

#ifndef __dwffld_w_h__
#define __dwffld_w_h__

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

void DwfFld(multi1d<LatticeFermion>& chi, const multi1d<LatticeFermion>& psi, enum PlusMinus isign);

} // End Namespace Chroma


#endif
