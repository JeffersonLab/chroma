// $Id: expmat.cc,v 1.7 2005-01-14 20:13:08 edwards Exp $
/*! \file
 *  \brief Exponentiate a SU(n) lie algebra element by some method,
 */

#include "chromabase.h"
#include "util/gauge/expmat.h"
#include "util/gauge/expm12.h"
#include "util/gauge/expsu3.h"
#include "util/gauge/eeu1.h"
#include "util/gauge/eesu2.h"
#include "util/gauge/eesu3.h"
#include "util/gauge/reunit.h"
//#include "util/gauge/eesu3.h"


namespace Chroma { 
//! Exponentiate a SU(n) lie algebra element by some method.
/*!
 * \ingroup gauge
 *
 * This routine is a driver for other routines. For example, expsu3
 * will exponentiate EFFICIENTLY an SU(3) matrix and, if desired,
 * check for errors in the truncation of the Taylor series.
 * The routine eesu3 EXACTLY exponentiates a SU(3) matrix. 
 *
 *  \param a        LatticeColorMatrix          (Modify)
 *  \param opt      Method of exponentiation    (Read)
 */

void expmat(LatticeColorMatrix& a,
	    enum ExpMat_t opt)
{
  START_CODE();

  switch (Nc)
  {
  case 1:
    eeu1(a);
    break;

#if 0
  case 2:
    /* Not tested yet */
    eesu2(a);
    break;
#endif

  default:
    switch (opt)
    {
    case EXP_TWELFTH_ORDER:
      switch (Nc)
      {
      case 3:
	expsu3(a, REUNITARIZE_LABEL);
	break;
	
      default:
	expm12(a);
	break;
      }
      break;

    case EXP_EXACT:
      switch (Nc)
      {
      case 3:
	eesu3(a);
	break;

      default:
	QDP_error_exit("exact exponentation not implemented for this Nc", Nc);
      }
      break;

    default:
      QDP_error_exit("unknown option", opt);
    }
  }

  END_CODE();
}

}; // End namespace Chroma
