// $Id: axgauge.cc,v 1.2 2004-07-28 02:38:03 edwards Exp $
/*! \file
 *  \brief Axial gauge fixing 
 */

#include "chromabase.h"
#include "meas/gfix/axgauge.h"

using namespace QDP;

//! Axial gauge fixing
/*!
 * \ingroup gfix
 *
 * Transfroms a gauge configuration, in place, into axial
 *            gauge with special direction j_decay.
 *
 * Note: The non-unity time-like gauge fields from the last time slice
 *       will be copied to all other time slices.
 *
 * \param ug         gauge field and its axial gauge transform (Modify)
 * \param j_decay    time direction (Read) 
 */

void axGauge(multi1d<LatticeColorMatrix>& ug, int j_decay)
{
  START_CODE();

  int lsizet = Layout::lattSize()[j_decay];

  LatticeColorMatrix v = 1;
  LatticeInteger t_coord = Layout::latticeCoordinate(j_decay);

  /* Transform the "time-like" links to unity, slice by slice except for the */
  /* last slice and thereby construct the gauge transformation V. */
  for(int t = 1; t < lsizet; ++t)
  {
    LatticeBoolean btmp = t_coord == t;

    LatticeColorMatrix tmp_1 = shift(ug[j_decay], BACKWARD, j_decay);
    copymask(v, btmp, tmp_1);
    LatticeColorMatrix tmp_2 = tmp_1 * ug[j_decay];
    copymask(ug[j_decay], btmp, tmp_2);
  }

  /* Now do the gauge transformation on the space-like links */
  for(int mu = 0; mu < Nd; ++mu)
  {
    if ( mu != j_decay )
    {
      LatticeColorMatrix tmp_2 = ug[mu] * shift(adj(v), FORWARD, mu);
      ug[mu] = v * tmp_2;
    }
  }

  /* Finally "broadcast" the t-link from the last time slice to all others */
  for(int t = lsizet-2; t >= 0; --t)
  {
    LatticeBoolean btmp = t_coord == t;
    
    copymask(ug[j_decay], btmp, LatticeColorMatrix(shift(ug[j_decay], FORWARD, j_decay)));
  }


  END_CODE();
}
