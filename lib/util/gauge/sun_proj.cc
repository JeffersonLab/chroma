//  $Id: sun_proj.cc,v 1.1 2003-03-30 15:35:22 edwards Exp $
/*! \file
 *  \brief Project a complex Nc x Nc matrix W onto SU(Nc) by maximizing Tr(VW)
 */

#include "chromabase.h"
#include "util/gauge/sun_proj.h"
#include "util/gauge/su3proj.h"
#include "util/gauge/reunit.h"

using namespace QDP;

//! Project a complex Nc x Nc matrix W onto SU(Nc) by maximizing Tr(VW)
/*!
 * Arguments:
 *
 *  \param w            complex Nc x Nc matrix (Read)
 *  \param v            the projected SU(Nc) Matrix (Write)
 *  \param BlkAccu      accuracy in SU(Nc) projection (Read)
 *  \param BlkMax       max number of iterations in SU(Nc) projection (Read)
 */

void sun_proj(const LatticeColorMatrix& w, LatticeColorMatrix& v,
	      const Real& BlkAccu, int BlkMax)
{
  Double new_tr;
  Double ddummy;

  START_CODE("sun_proj");

  /*
   * Project back to SU(3) by maximizing tr(v w).
   * This is done by looping proj_iter times over the 3 SU(2) subgroups.
   */

  /* The initial trace */
  Double old_tr = sum(real(trace(v * w))) / double(Layout::vol()*Nc);

  int iter = 0;
  int wrswitch = 1;			/* Write out iterations? */
//  Double conver = 1.0;
  Real conver = 1.0;

  while ( toBool(conver > BlkAccu)  &&  iter < BlkMax )
  {
    iter = iter + 1;

    // Loop over SU(2) subgroup index
    for(int su2_index = 0; su2_index < Nc*(Nc-1)/2; ++su2_index)
      su3proj(v, w, su2_index);

    // Reunitarize
    reunit(v);

    // Calculate the trace
    new_tr = sum(real(trace(v * w))) / double(Layout::vol()*Nc);

    if( wrswitch == 1 )
    {
      cerr << "iter = " << iter << endl;
      cerr << "old_tr = " << old_tr << endl;
      cerr << "new_tr = " << new_tr << endl;
    }

    // Normalized convergence criterion:
    conver = fabs((new_tr - old_tr) / old_tr);
    old_tr = new_tr;
  }

  if ( wrswitch == 1 )
  {
//    Push(nml,"Final_sun_proj");
//    Write(nml, iter);
//    Write(nml, new_tr);
//    Pop(nml);
    cerr << "iter = " << iter << endl;
    cerr << "new_tr = " << new_tr << endl;
  }

  END_CODE("sun_proj");
}
