// $Id: expmat.cc,v 1.2 2004-01-05 00:47:20 edwards Exp $
/*! \file
 *  \brief Exponentiate a SU(n) lie algebra element by some method,
 */

#include "chromabase.h"
#include "util/gauge/expmat.h"
#include "util/gauge/expm12.h"
#include "util/gauge/eeu1.h"
#include "util/gauge/eesu2.h"
//#include "util/gauge/eesu3.h"

using namespace QDP;

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
	    enum Expmat_t opt)
{
  START_CODE("expmat");

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
    case EXP_TWELTH_ORDER:
      switch (Nc)
      {
      case 3:
	expsu3(a, EXP_TWELTH_ORDER);
	break;
	
      default:
	expm12(a);
	break;
      }
      break;

    case EXP_EXACT:
      switch (Nc)
      {
#if 0
      case 3:
	/* Should be rechecked */
	eesu3 (a);
	break;
#endif

      default:
	QDP_error_exit("exact exponentation not implemented for this Nc", Nc);
      }
      break;

    default:
      QDP_error_exit("unknown option", opt);
    }
  }

  END_CODE("expmat");
}
