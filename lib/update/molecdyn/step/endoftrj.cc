// $Id: endoftrj.cc,v 1.1 2003-12-30 19:52:28 edwards Exp $

#error "NOT FULLY CONVERTED - NEED TO MOVE AlgETrj into params of Integ. functor"

#include "chromabase.h"

using namespace QDP;

//! Test for end of a trajectory
/*!
 * \ingroup molecdyn
 *
 * \param t      integration time so far ( Read )
 * \return EndP  flag indicating end of trajectory ( Write ) 
 */

bool EndOfTrj(const Real& t)
{
  Real r;
  
  START_CODE("EndOfTrj");
  
  switch (AlgETrj)
  {
  case FIXED_LENGTH:
    if ( t > tau0+dt/TO_REAL(2) )
      EndP = true;
    else
      EndP = false;

    break;
  case EXPONENTIAL_LENGTH:
    random(r);
    if ( r < dt/max(tau0-dt,dt) )
      EndP = true;
    else
      EndP = false;

    break;
  default:
    QDP_error_exit("unknown algorithm termination", AlgETrj);
  }
  
  END_CODE("EndOfTrj");
}
