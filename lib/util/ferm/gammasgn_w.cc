// $Id: gammasgn_w.cc,v 1.1 2004-06-11 20:31:25 edwards Exp $
/*! \file
 *  \brief Compute gamma matrix multiplication table factors
 */

#include "chromabase.h"
#include "util/ferm/gammasgn.h"

using namespace QDP;

static multi2d<int> meson_eta2;
static bool initP = false;


//! Init gamma matrix multiplication table factors
/*!
 * \ingroup ferm
 *
 * Initialize signs needed for  Gamma(n)*Gamma(m)=sgn(n,m)*Gamma(n ^ m)
 */

static void gammaSgn_init()
{
  START_CODE("gammasgn_init");

  sgnn.resize(Ns*Ns,Ns*Ns);

  multi1d sgnn(2);
  
  sgnn[0] =  1;
  sgnn[1] = -1;
  
  for(int n = 0; n < Ns*Ns; ++n)
  {
    int nn = n;
    int n1 = nn & 1;
    nn >>= 1;
    int n2 = nn & 1;
    nn >>= 1;
    int n3 = nn & 1;
    int n4 = nn >> 1;

    int ssum = (n1 + n2 + n3 + n4) & 1;

    for(int m = 0; m < Ns*Ns; ++m)
    {
      int mm;
      int m1 = mm & 1;
      mm >>= 1;
      m2 = mm & 1;
      mm >>= 1;
      m3 = mm & 1;
      m4 = mm >> 1;

      ssum = ((n2+n3+n4)*m1 + (n3+n4)*m2 + n4*m3) & 1;
      meson_eta2(n,m) = sgnn[ssum];
    }
  }
  
  initP = true;

  END_CODE("gammasgn_init");
}

//! Return gamma matrix multiplication table factors
/*!
 * \ingroup ferm
 *
 * Initialize signs needed for  Gamma(n)*Gamma(m)=sgn(n,m)*Gamma(n ^ m)
 */

int gammaSgn(int n, int m)
{
  if (! initP)
    gammaSgn_init();

  return meson_eta2(n,m);
}

