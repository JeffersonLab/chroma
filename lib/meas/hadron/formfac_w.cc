// $Id: formfac_w.cc,v 1.4 2003-03-06 00:25:23 flemingg Exp $
/*! \file
 *  \brief Form-factors 
 *
 * Form factors constructed from a quark and a sequential quark propagator
 *
 * $Log: formfac_w.cc,v $
 * Revision 1.4  2003-03-06 00:25:23  flemingg
 * Added some comments about things that need to be looked at more carefully.
 *
 * Revision 1.3  2003/03/05 01:25:47  flemingg
 * Modified to insert all Nd^2 Gamma matrices in the local current.  For
 * vector and axial vector insertions, the non-local currents are computed
 * as well.
 *
 */

#include "chromabase.h"
#include "meas/hadron/formfac_w.h"
#include "proto.h"                  // part of QDP++, for crtesn()

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
 * \param source_mom2_max -- max source hadron mom squared ( Read )
 * \param t_sink -- time coordinate of the sink ( Read )
 * \param sink_mom -- sink hadron momentum ( Read )
 * \param j_decay -- direction of the exponential decay ( Read ) 
 * \param nml   -- namelist file object ( Read )
 */

void FormFac(const multi1d<LatticeColorMatrix>& u, 
	     const LatticePropagator& quark_propagator,
	     const LatticePropagator& seq_quark_prop, 
	     const multi1d<int>& t_source, 
	     int source_mom2_max,
	     int t_sink,
	     const multi1d<int>& sink_mom, 
	     int j_decay,
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
   * Coordinates for source momenta
   */
  multi1d<LatticeInteger> my_coord(Nd);
  for(int mu = 0; mu < Nd; ++mu)
    my_coord[mu] = Layout::latticeCoordinate(mu);

  
  /*
   * Construct the anti-quark propagator from the seq. quark prop.
   */ 
  LatticePropagator anti_quark_prop = Gamma(G5) * seq_quark_prop * Gamma(G5);


  /*
   * Loop over gamma matrices of the insertion current of insertion current
   */
  for(int gamma_value = 0; gamma_value < Nd*Nd; ++gamma_value)
  {
    /*
     *  For the case where the gamma value indicates we are evaluating either
     *  the vector or axial vector currents, we will also evaluate
     *  the non-local currents.  The non-local vector current is the conserved
     *  current.  The non-local axial vector current would be partially
     *  conserved but for the Wilson term.  In these cases we will set
     *  mu = corresponding direction.  In all other cases, we will set mu = -1.
     */

    int compute_nonlocal;
    int mu;

    switch(gamma_value){
    case  1:
    case 14:
      mu = 0;
      compute_nonlocal = 1;
      break;
    case  2:
    case 13:
      mu = 1;
      compute_nonlocal = 1;
      break;
    case  4:
    case 11:
      mu = 2;
      compute_nonlocal = 1;
      break;
    case  8:
    case  7:
      mu = 3;
      compute_nonlocal = 1;
      break;
    default:
      mu = -1;
      compute_nonlocal = 0;
    }

    /* 
     * The local non-conserved vector-current matrix element 
     */
    LatticeComplex corr_local_fn =
      trace(adj(anti_quark_prop) * Gamma(gamma_value) * quark_propagator);

    /*
     * Construct the non-local current matrix element 
     *
     * The form of J_mu = (1/2)*[psibar(x+mu)*U^dag_mu*(1+gamma_mu)*psi(x) -
     *                           psibar(x)*U_mu*(1-gamma_mu)*psi(x+mu)]
     * NOTE: the 1/2  is included down below in the sumMulti stuff
     */
    LatticeComplex corr_nonlocal_fn;
    if(compute_nonlocal){
      corr_nonlocal_fn =
        trace(adj(u[mu] * shift(anti_quark_prop, FORWARD, mu)) *
          (quark_propagator + Gamma(gamma_value) * quark_propagator));
      LatticePropagator tmp_prop1 = u[mu] *
        shift(quark_propagator, FORWARD, mu);
      corr_nonlocal_fn -= trace(adj(anti_quark_prop) *
                            (tmp_prop1 - Gamma(gamma_value) * tmp_prop1));
    }

    /*
     * Loop over allowed source momenta: (source_mom)^2 <= source_mom2_max.
     * Do this by constructing a L^(Nd-1) grid in momenta centered about the
     * origin. Loop lexicographically over all the "sites" (momenta value)
     * and toss out ones considered too large to give any reasonable statistics
     *
     * NOTE: spatial anisotropy is no allowed here
     */
    int Ndm1 = Nd-1;
    multi1d<int> mom_size(Ndm1);
    int L ;
    int mom_vol = 1;

    for (L=0; (L+1)*(L+1) <= source_mom2_max; ++L) ;

    for(int nu=0; nu < Ndm1; ++nu)
    {
      mom_vol *= (2*L)+1;
      mom_size[nu] = (2*L)+1;
    }

    for(int n = 0; n < mom_vol; ++n)
    {
      multi1d<int> source_mom = crtesn(n, mom_size);
      multi1d<int> inser_mom(Ndm1) ;

      int source_mom2 = 0;
      for(int nu = 0; nu < Ndm1; ++nu)
      {
	source_mom[nu] -= L;
	source_mom2 += source_mom[nu]*source_mom[nu];
      }

      // Arbitrarily set the cutoff on max allowed momentum to be [2,1,0]
      if (source_mom2 > source_mom2_max) continue;

      for(int nu = 0; nu < Ndm1; ++nu)
        inser_mom[nu] = sink_mom[nu] - source_mom[nu] ;

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

      // The complex phases exp(i p.x).  GTF CHECK p vs. -p
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
      if(compute_nonlocal){
        hsum = sumMulti(phasefac * corr_nonlocal_fn, timeslice);

        for(int t = 0; t < length; ++t)
        {
          int t_eff = (t - t0 + length) % length;

          nonlocal_cur3ptfn[t_eff] = 0.5 * Complex(hsum[t]);
        }
      }

      // Print out the results
      push(nml,"Wilson_Current_3Pt_fn");
      Write(nml,gamma_value);
      // GTF comment: we should decide here whether to always write mu
      //   or only when mu is meaningful, i.e. 0 <= mu < Nd
      Write(nml,mu);
      Write(nml,j_decay);
      Write(nml,inser_mom);
      Write(nml,local_cur3ptfn);
      if(compute_nonlocal){ Write(nml,nonlocal_cur3ptfn); }
      pop(nml);
    }
  }
                            
//  END_CODE();
}
