// $Id: temporal_gauge.cc,v 3.1 2009-10-09 15:33:43 bjoo Exp $
/*! \file
 *  \brief Axial gauge fixing 
 */

#include "chromabase.h"
#include "meas/gfix/temporal_gauge.h"

namespace Chroma 
{

  //! Temporal gauge fixing
  /*!
   * \ingroup gfix
   *
   * Transfroms a gauge configuration, in place, into axial
   *            gauge with special direction decay_dir.
   *
   * Note: The non-unity time-like gauge fields from the last time slice
   *       will be copied to all other time slices.
   *
   * \param ug         gauge field and its axial gauge transform (Modify)
   * \param g          gauge rotation matrix (Write)
   * \param decay_dir  time direction (Read) 
   */

  void temporalGauge(multi1d<LatticeColorMatrix>& ug, LatticeColorMatrix& g, int decay_dir)
  {
    START_CODE();

    // Sanity check
    if (decay_dir < 0 || decay_dir >= ug.size())
    {
      QDPIO::cerr << __func__ << ": invalid decay_dir= " << decay_dir << std::endl;
      QDP_abort(1);
    }

    // Initialize
    int N_t = Layout::lattSize()[decay_dir];
    LatticeInteger t_coord = Layout::latticeCoordinate(decay_dir);

    {
      LatticeColorMatrixD u = ug[decay_dir];
      LatticeColorMatrixD gt = 1;

      // The first shift
      LatticeColorMatrixD U_prev = shift(u, BACKWARD, decay_dir);
      LatticeColorMatrixD G_prev = gt; // Not shifting becuse it is unit at this point

      for(int t = 1; t < N_t; ++t) {
	// Mask timeslice t
	LatticeBoolean btmp = (t_coord == t); 
	LatticeColorMatrixD t1 = G_prev*U_prev ; // DO whole lattice... That Sucks
	copymask(gt, btmp, t1);
	
	// Reuse t1:
	t1 = shift(u,BACKWARD, decay_dir);
	U_prev = t1;
	t1 = shift(gt, BACKWARD, decay_dir);
	G_prev = t1;
      }
      g=gt;
    }

    /* Now do the gauge transformation on all the links */
    for(int mu = 0; mu < Nd; ++mu) {

      LatticeColorMatrix tmp = ug[mu] * shift(adj(g), FORWARD, mu);
      ug[mu] = g * tmp;
    }

#if 0
    // Check temporal links (except last timeslice are unit
    {
      // Take a copy
      LatticeColorMatrix u = ug[decay_dir];
      LatticeColorMatrix g_unit = 1; // For copying into the last timeslice
      
      // Blast last non-unit timeslice with a unit one...
      LatticeBoolean btmp  = (t_coord == (N_t-1) );
      copymask(u, btmp, g_unit); // Replace last timeslice with units...

      // U should now be indistinguishable from unit gauge
      QDPIO::cout << "Norm of Unit-Gauge violation / link = " << sqrt( norm2(u-g_unit) )/ Layout::vol() << endl;

    }
#endif

    END_CODE();
  }

}  // end namespace Chroma
