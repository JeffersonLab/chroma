//  $Id: mesons_w.cc,v 1.1 2003-01-04 04:10:56 edwards Exp $


#include <szin.h>

using namespace QDP;

//! Function used for constructing the time-slice set
static const int j_decay = Nd-1;
static int set_timeslice_func(const multi1d<int>& coordinate) {return coordinate[j_decay];}
 
//! Meson 2-pt functions
/* This routine is specific to Wilson fermions!
 *
 * Construct meson propagators
 * The two propagators can be identical or different.
 *
 * quark_prop_1 -- first quark propagator ( Read )
 * quark_prop_2 -- second (anti-) quark propagator ( Read )
 * meson_propagator -- Ns^2 mesons ( Modify )
 * t_source -- cartesian coordinates of the source ( Read )
 * j_decay -- direction of the exponential decay ( Read )
 *
 *        ____
 *        \
 * m(t) =  >  < m(t_source, 0) m(t + t_source, x) >
 *        /
 *        ----
 *          x
 */

void mesons(const LatticePropagator& quark_prop_1, const LatticePropagator& quark_prop_2, 
	    multi2d<Real>& meson_propagator, 
	    const multi1d<int>& t_source)
{
  int length = Layout::lattSize()[j_decay];
  int t0 = t_source[j_decay];
  int G5 = Ns*Ns-1;
  multi1d<Double> hsum(length);

  // Create the time-slice set
  Set timeslice(set_timeslice_func, length);

  // Initialize the propagator so that we just add to it below
  meson_propagator = 0.0;

  // Contruct the antiquark prop
  LatticePropagator anti_quark_prop =  Gamma(G5) * quark_prop_2 * Gamma(G5);

  for(int n = 0; n < (Ns*Ns); ++n)
  {
    LatticeReal psi_sq = real(trace(conj(anti_quark_prop) * Gamma(n) * quark_prop_1 * Gamma(n)));

    // Do a slice-wise sum.
    hsum = sumMulti(psi_sq, timeslice);

    for(int t = 0; t < length; ++t)
    {
      int t_eff = (t - t0 + length) % length;

      meson_propagator[n][t_eff] += Real(hsum[t]);
    }
  }
}
