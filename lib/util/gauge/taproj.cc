// $Id: taproj.cc,v 1.2 2003-01-26 04:30:41 edwards Exp $
// TAPROJ

/*! \file
 *  \brief Take the traceless antihermitian projection of a color matrix
 */

/*!
 *  a = (1/2)[a - a_dag] - Tr[(1/2)*(a - a_dag)]/Nc
 *
 * that is the anti-hermitian traceless part of a 
 *
 * Arguments:
 *
 *  \param a        LatticeColorMatrix          (Modify)
 */

#include "qdp.h"
#include "szin.h"

using namespace QDP;

void taproj(LatticeColorMatrix& a)
{
  START_CODE("taproj");

  // a = a - a_dagger  --- a -> antihermitian matrix
  LatticeColorMatrix aux_1 = a;
  a -= adj(aux_1);
 
  if (Nc > 1)
  {
    // tmp = Tr[ a ] * (1/Nc)
    LatticeReal tmp = imag(trace(a))*Real(1/Real(Nc));

    // a = a - (1/Nc) * Tr[ a] --- a -> traceless matrix
    a -= cmplx(0,tmp);
  }

  // Normalisation to make taproj idempotent
  a *= 0.5;
  
  END_CODE("taproj");
}
