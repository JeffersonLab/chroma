// $Id: curcor2_w.cc,v 3.0 2006-04-03 04:58:59 edwards Exp $
/*! \file
 *  \brief Mesonic current correlators
 */

#include "chromabase.h"
#include "util/ft/sftmom.h"
#include "meas/hadron/mesons_w.h"
#include "qdp_util.h"                 // part of QDP++, for crtesn()

namespace Chroma {

//! Construct current correlators 
/*!
 * \ingroup hadron
 *
 * This routine is specific to Wilson fermions!
 *
 *  The two propagators can be identical or different.

 * This includes the "rho_1--rho_2" correlators used for O(a) improvement

 * For use with "rotated" propagators we added the possibility of also
 * computing the local vector current, when no_vec_cur = 4. In this
 * case the 3 local currents come last.

 * \param u               gauge field ( Read )
 * \param quark_prop_1    first quark propagator ( Read )
 * \param quark_prop_2    second (anti-) quark propagator ( Read )
 * \param phases          fourier transform phase factors ( Read )
 * \param t0              timeslice coordinate of the source ( Read )
 * \param no_vec_cur      number of vector current types, 3 or 4 ( Read )
 * \param xml             namelist file object ( Read )
 * \param xml_group       string used for writing xml data ( Read )
 *
 *         ____
 *         \
 * cc(t) =  >  < m(t_source, 0) c(t + t_source, x) >
 *         /                    
 *         ----
 *           x
 */


void curcor2(const multi1d<LatticeColorMatrix>& u, 
	     const LatticePropagator& quark_prop_1, 
	     const LatticePropagator& quark_prop_2, 
	     const SftMom& phases,
	     int t0,
	     int no_vec_cur,
	     XMLWriter& xml,
	     const string& xml_group)
{
  START_CODE();

  if ( no_vec_cur < 2 || no_vec_cur > 4 )
    QDP_error_exit("no_vec_cur must be 2 or 3 or 4", no_vec_cur);

  // Initial group
  push(xml, xml_group);

  write(xml, "num_vec_cur", no_vec_cur);

  // Length of lattice in decay direction
  int length  = phases.numSubsets();
  int j_decay = phases.getDir();

  LatticePropagator tmp_prop1;
  LatticePropagator tmp_prop2;

  LatticeReal psi_sq;
  LatticeReal chi_sq;

  multi1d<Double> hsum(length);

  // Construct the anti-quark propagator from quark_prop_2
  int G5 = Ns*Ns-1;
  LatticePropagator anti_quark_prop =  Gamma(G5) * quark_prop_2 * Gamma(G5);

  // Vector currents
  {
    multi2d<Real> vector_current(no_vec_cur*(Nd-1), length);

    /* Construct the 2*(Nd-1) non-local vector-current to rho correlators */
    int kv = -1;
    int kcv = Nd-2;

    for(int k = 0; k < Nd; ++k)
    {
      if( k != j_decay )
      {
	int n = 1 << k;
	kv = kv + 1;
	kcv = kcv + 1;

	tmp_prop2 = u[k] * shift(quark_prop_1, FORWARD, k) * Gamma(n);
	chi_sq = - real(trace(adj(anti_quark_prop) * tmp_prop2));

	tmp_prop1 = Gamma(n) * tmp_prop2;
	psi_sq = real(trace(adj(anti_quark_prop) * tmp_prop1));

	tmp_prop2 = u[k] * shift(anti_quark_prop, FORWARD, k) * Gamma(n);
	chi_sq += real(trace(adj(tmp_prop2) * quark_prop_1));

	tmp_prop1 = Gamma(n) * tmp_prop2;
	psi_sq += real(trace(adj(tmp_prop1) * quark_prop_1));

	chi_sq += psi_sq;

	/* Do a slice-wise sum. */

//    Real dummy1 = Real(meson_eta[n]) / Real(2);
	Real dummy1 = 0.5;

	/* The nonconserved vector current first */
	hsum = sumMulti(psi_sq, phases.getSet());

	for(int t = 0; t < length; ++t)
	{
	  int t_eff = (t - t0 + length) % length;
	  vector_current[kv][t_eff] = dummy1 * Real(hsum[t]);
	}

	/* The conserved vector current next */
	hsum = sumMulti(chi_sq, phases.getSet());

	for(int t = 0; t < length; ++t)
	{
	  int t_eff = (t - t0 + length) % length;
	  vector_current[kcv][t_eff] = dummy1 * Real(hsum[t]);
	}
      }
    }


    /* Construct the O(a) improved vector-current to rho correlators,
       if desired */
    if ( no_vec_cur >= 3 )
    {
      kv = 2*Nd-3;
      int jd = 1 << j_decay;

      for(int k = 0; k < Nd; ++k)
      {
	if( k != j_decay )
	{
	  int n = 1 << k;
	  kv = kv + 1;
	  int n1 = n ^ jd;

	  psi_sq = real(trace(adj(anti_quark_prop) * Gamma(n1) * quark_prop_1 * Gamma(n)));

//	dummy1 = - Real(meson_eta[n]);
	  Real dummy1 = - 1;

	  /* Do a slice-wise sum. */
	  hsum = sumMulti(psi_sq, phases.getSet());

	  for(int t = 0; t < length; ++t)
	  {
	    int t_eff = (t - t0 + length) % length;
	    vector_current[kv][t_eff] = dummy1 * Real(hsum[t]);
	  }
	}
      }
    }


    /* Construct the local vector-current to rho correlators, if desired */
    if ( no_vec_cur >= 4 )
    {
      kv = 3*Nd-4;

      for(int k = 0; k < Nd; ++k)
      {
	if( k != j_decay )
	{
	  int n = 1 << k;
	  kv = kv + 1;

	  psi_sq = real(trace(adj(anti_quark_prop) * Gamma(n) * quark_prop_1 * Gamma(n)));

//	dummy1 = Real(meson_eta[n]);
	  Real dummy1 = 1;

	  /* Do a slice-wise sum. */
	  hsum = sumMulti(psi_sq, phases.getSet());

	  for(int t = 0; t < length; ++t)
	  {
	    int t_eff = (t - t0 + length) % length;
	    vector_current[kv][t_eff] = dummy1 * Real(hsum[t]);
	  }
	}
      }
    }


    // Loop over currents to print
    XMLArrayWriter xml_cur(xml,vector_current.size2());
    push(xml_cur, "Vector_currents");

    for (int current_value=0; current_value < vector_current.size2(); ++current_value)
    {
      push(xml_cur);     // next array element

      write(xml_cur, "current_value", current_value);
      write(xml_cur, "vector_current", vector_current[current_value]);

      pop(xml_cur);
    }

    pop(xml_cur);
  }


  //
  // Axial currents
  //
  {
    multi2d<Real> axial_current(2, length);

    /* Construct the 2 axial-current to pion correlators */
    int n = G5 ^ (1 << j_decay);

    /* The local axial current first */
    psi_sq = real(trace(adj(anti_quark_prop) * Gamma(n) * quark_prop_1 * Gamma(G5)));

    /* The nonlocal axial current next */
    chi_sq  = real(trace(adj(anti_quark_prop) * Gamma(n) * 
			 u[j_decay] * shift(quark_prop_1, FORWARD, j_decay) * Gamma(G5)));

    // The () forces precedence
    chi_sq -= real(trace(adj(Gamma(n) * (u[j_decay] * shift(anti_quark_prop, FORWARD, j_decay)) * 
			     Gamma(G5)) * quark_prop_1));

    /* Do a slice-wise sum. */

    Real dummy1 = Real(-1) / Real(2);

    /* The local axial current first */
    hsum = sumMulti(psi_sq, phases.getSet());

    for(int t = 0; t < length; ++t)
    {
      int t_eff = (t - t0 + length) % length;
      axial_current[1][t_eff] = - Real(hsum[t]);
    }

    /* The nonlocal axial current next */
    hsum = sumMulti(chi_sq, phases.getSet());

    for(int t = 0; t < length; ++t)
    {
      int t_eff = (t - t0 + length) % length;
      axial_current[0][t_eff] = dummy1 * Real(hsum[t]);
    }


    // Loop over currents to print
    XMLArrayWriter xml_cur(xml,axial_current.size2());
    push(xml_cur, "Axial_currents");

    for (int current_value=0; current_value < axial_current.size2(); ++current_value)
    {
      push(xml_cur);     // next array element

      write(xml_cur, "current_value", current_value);
      write(xml_cur, "axial_current", axial_current[current_value]);

      pop(xml_cur);
    }

    pop(xml_cur);
  }


  pop(xml);  // xml_group
              
  END_CODE();
}

}  // end namespace Chroma
