// $Id: expmat.cc,v 3.3 2009/11/14 20:01:46 eneil Exp $
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

  namespace ExpMatEnv {
    static double time_spent = 0;

    double getTime() { return time_spent ; }

  };


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
  QDP::StopWatch swatch;
  swatch.reset();
  swatch.start();

  switch (Nc)
  {
  case 1:
    eeu1(a);
    break;

  case 2:
    eesu2(a);
    break;

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
	QDPIO::cerr << __func__ << ": exact exponentation not implemented for this Nc=" << Nc << endl;
	QDP_abort(1);
      }
      break;

    default:
      QDPIO::cerr << __func__ << ": unknown option = " << opt << endl;
      QDP_abort(1);
    }
  }

  swatch.stop();
  ExpMatEnv::time_spent += swatch.getTimeInSeconds();
  END_CODE();
}

}; // End namespace Chroma
