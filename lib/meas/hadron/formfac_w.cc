// $Id: formfac_w.cc,v 1.1 2003-02-25 20:25:28 edwards Exp $
/*! \file
 *  \brief Form-factors 
 *
 * Form factors constructed from a quark and a sequential quark propagator
 */

#include "chromabase.h"
#include "meas/hadron/formfac_w.h"

using namespace QDP;

//! Function object used for constructing the time-slice set
class TimeSliceFunc : public SetFunc
{
public:
  TimeSliceFunc(int dir): dir_decay(dir) {}

  int operator() (const multi1d<int>& coordinate) const {return coordinate[dir_decay];}
  int numSubsets() const {return Layout::lattSize()[dir_decay];}

  int dir_decay;

private:
  TimeSliceFunc() {}  // hide default constructor
};

 
//! Compute contractions for current insertion 3-point functions.
/*!
 *
 * This routine is specific to Wilson fermions!
 *
 * \param u    -- gauge fields (used for non-local currents) ( Read )
 * \param quark_propagator -- quark propagator ( Read )
 * \param seq_quark_prop -- sequential quark propagator ( Read )
 * \param t_source -- cartesian coordinates of the source ( Read )
 * \param t_sink -- time coordinate of the sink ( Read )
 * \param j_decay -- direction of the exponential decay ( Read ) 
 * \param nml   -- namelist file object ( Read )
 */

void FormFac(const multi1d<LatticeColorMatrix>& u, 
	     const LatticePropagator& quark_propagator,
	     const LatticePropagator& seq_quark_prop, 
	     const multi1d<int>& t_source, 
	     int t_sink, int j_decay,
	     NmlWriter& nml)
{
  // Create the time-slice set
  Set timeslice;
  timeslice.make(TimeSliceFunc(j_decay));

  // Length of lattice in j_decay direction and 3pt correlations fcns
  int length = timeslice.numSubsets();
  multi1d<Complex> local_cur3ptfn(length);
  multi1d<Complex> nonlocal_cur3ptfn(length);
  
//  START_CODE("formfac");

  int t0 = t_source[j_decay];
  int G5 = Ns*Ns-1;
  
  /*
   * Coordinates for insertion momenta
   */
  multi1d<LatticeInteger> my_coord(Nd);
  for(int mu = 0; mu < Nd; ++mu)
    my_coord[mu] = Layout::latticeCoordinate(mu);

  
  /*
   * Construct the anti-quark propagator from the seq. quark prop.
   */ 
  LatticePropagator anti_quark_prop = Gamma(G5) * seq_quark_prop * Gamma(G5);


  /*
   * Loop over mu of insertion current
   */
  for(int mu = 0; mu < Nd; ++mu)
  {
    /* 
     * The local non-conserved vector-current matrix element 
     */
    int n = 1 << mu;
    LatticeComplex corr_local_fn = trace(adj(anti_quark_prop) * Gamma(n) * quark_propagator);

    /* 
     * Construct the conserved non-local vector-current matrix element 
     * NOTE: this is only really conserved for the Wilson fermion quark action
     *
     * The form of J_mu = (1/2)*[psibar(x+mu)*U^dag_mu*(1+gamma_mu)*psi(x) -
     *                           psibar(x)*U_mu*(1-gamma_mu)*psi(x+mu)]
     * NOTE: the 1/2  is included down below in the sumMulti stuff
     */
    LatticeComplex corr_nonlocal_fn = 
      trace(adj(u[mu] * shift(anti_quark_prop, FORWARD, mu)) * 
	    (quark_propagator + Gamma(n) * quark_propagator));

    LatticePropagator tmp_prop1 = u[mu] * shift(quark_propagator, FORWARD, mu);
    corr_nonlocal_fn -= trace(adj(anti_quark_prop) * (tmp_prop1 - Gamma(n) * tmp_prop1));


    /*
     * Loop over non-zero insertion momenta
     * Do this by constructing a 5^(Nd-1) grid in momenta centered about the
     * origin. Loop lexicographically over all the "sites" (momenta value)
     * and toss out ones considered too large to give any reasonable statistics
     *
     * NOTE: spatial anisotropy is no allowed here
     */
    multi1d<int> mom_size(Nd-1);
    int Ndm1 = Nd-1;
    int L = 5;
    int mom_vol = 1;

    for(int nu=0; nu < Ndm1; ++nu)
    {
      mom_vol *= L;
      mom_size[nu] = L;
    }

    for(int n = 0; n < mom_vol; ++n)
    {
      multi1d<int> inser_mom = crtesn(n, mom_size);

      int q_sq = 0;
      for(int nu = 0; nu < Ndm1; ++nu)
      {
	inser_mom[nu] -= (L-1)/2;
	q_sq += inser_mom[nu]*inser_mom[nu];
      }

      // Arbitrarily set the cutoff on max allowed momentum to be [2,1,0]
      if (q_sq > 4) continue;

      LatticeReal p_dot_x(float(0.0));

      int j = 0;
      for(int nu = 0; nu < Nd; ++nu)
      {
	const Real twopi = 6.283185307179586476925286;

	if (nu == j_decay)
	  continue;

	p_dot_x += LatticeReal(my_coord[nu]) * twopi
	  * Real(inser_mom[j]) / Layout::lattSize()[nu];
	j++;
      }

      // The complex phases  exp(i p.x )
      LatticeComplex  phasefac = cmplx(cos(p_dot_x), sin(p_dot_x));

      // Local current
      multi1d<DComplex> hsum(length);
      hsum = sumMulti(phasefac * corr_local_fn, timeslice);

      for(int t = 0; t < length; ++t)
      {
	int t_eff = (t - t0 + length) % length;

	local_cur3ptfn[t_eff] = Complex(hsum[t]);
      }


      // Non-local current
      hsum = sumMulti(phasefac * corr_nonlocal_fn, timeslice);

      for(int t = 0; t < length; ++t)
      {
	int t_eff = (t - t0 + length) % length;

	nonlocal_cur3ptfn[t_eff] = 0.5 * Complex(hsum[t]);
      }

      // Print out the results
      push(nml,"Wilson_Current_3Pt_fn");
      Write(nml,mu);
      Write(nml,j_decay);
      Write(nml,inser_mom);
      Write(nml,local_cur3ptfn);
      Write(nml,nonlocal_cur3ptfn);
      pop(nml);
    }
  }
                            
//  END_CODE();
}
