// $Id: curcor2_w.cc,v 1.2 2003-09-30 21:01:04 edwards Exp $
/*! \file
 *  \brief Mesonic current correlators
 */

#include "chromabase.h"
#include "util/ft/sftmom.h"
#include "meas/hadron/mesons_w.h"
#include "qdp_util.h"                 // part of QDP++, for crtesn()

using namespace QDP;

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
 * \param t0              timeslice coordinate of the source ( Read )
 * \param j_decay         direction of the exponential decay ( Read )
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
	     int j_decay, int no_vec_cur,
	     XMLWriter& xml,
	     const string& xml_group)
{
  START_CODE("curcor2");

  if ( no_vec_cur < 2 || no_vec_cur > 4 )
    QDP_error_exit("no_vec_cur must be 2 or 3 or 4", no_vec_cur);

  // Beginning group
  push(xml, xml_group);

  // Length of lattice in decay direction
  int length = phases.numSubsets() ;

  LatticePropagator tmp_prop1;
  LatticePropagator tmp_prop2;

  multi1d<Double> hsum(length);

  // Construct the anti-quark propagator from quark_prop_2
  int G5 = Ns*Ns-1;
  LatticePropagator anti_quark_prop =  Gamma(G5) * quark_prop_2 * Gamma(G5);
        
  /* Construct the 2*(Nd-1) non-local vector-current to rho correlators */
  LatticeReal psi_sq;
  LatticeReal chi_sq;

  int kv = -1;
  int kcv = Nd-2;

  {
    XMLArrayWriter xml_dir(xml,Nd-1);
    push(xml_dir, "Vector_current");

    for(int k = 0; k < Nd; ++k)
    {
      if( k != j_decay )
      {
	push(xml_dir);     // next array element
	Write(xml_dir, k);

	int n = 1 << k;
	kv = kv + 1;
	kcv = kcv + 1;
	psi_sq = 0;
	chi_sq = 0;

	tmp_prop2 = u[k] * shift(quark_prop_1, FORWARD, k) * Gamma(n);
	chi_sq -= real(trace(adj(anti_quark_prop) * tmp_prop2));

	tmp_prop1 = Gamma(n) * tmp_prop2;
	psi_sq += real(trace(adj(anti_quark_prop) * tmp_prop1));

	tmp_prop2 = u[k] * shift(anti_quark_prop, FORWARD, k) * Gamma(n);
	chi_sq += real(trace(adj(tmp_prop2) * quark_prop_1));

	tmp_prop1 = Gamma(n) * tmp_prop2;
	psi_sq += real(trace(adj(tmp_prop1) * quark_prop_1));

	/* Do a slice-wise sum. */

//	Real dummy1 = Real(meson_eta[n]) / Real(2);
	Real dummy1 = Real(1) / Real(2);

	/* The nonconserved vector current first */

	hsum = sumMulti(psi_sq, phases.getSubset());

	multi1d<Real> nonconserved_vector_current(length);
	for(int t = 0; t < length; ++t)
	{
	  int t_eff = (t - t0 + length) % length;
	  nonconserved_vector_current[t_eff] += dummy1 * Real(hsum[t]);
	}
	Write(xml_dir, nonconserved_vector_current);

	/* The conserved vector current next */

	hsum = sumMulti(chi_sq, phases.getSubset());

	multi1d<Real> conserved_vector_current(length);
	for(int t = 0; t < length; ++t)
	{
	  int t_eff = (t - t0 + length) % length;
	  conserved_vector_current[t_eff] += dummy1 * Real(hsum[t]);
	}
	Write(xml_dir, conserved_vector_current);
      
	pop(xml_dir);
      }
    }

    pop(xml_dir);
  }

  /* Construct the O(a) improved vector-current to rho correlators,
     if desired */

  if ( no_vec_cur >= 3 )
  {
    XMLArrayWriter xml_dir(xml,Nd-1);
    push(xml_dir, "Oa_improved_vector_current");

    kv = 2*Nd-3;
    int jd = 1 << j_decay;

    for(int k = 0; k < Nd; ++k)
    {
      if( k != j_decay )
      {
	int n = 1 << k;
	kv = kv + 1;
	psi_sq = 0;
	int n1 = n ^ jd;

	psi_sq = real(trace(adj(anti_quark_prop) * Gamma(n1) * quark_prop_1 * Gamma(n)));

	/* Do a slice-wise sum. */
	hsum = sumMulti(psi_sq, phases.getSubset());
//	dummy1 = - Real(meson_eta[n]);
	Real dummy1 = - 1;

	multi1d<Real> improved_vector_current(length);
	for(int t = 0; t < length; ++t)
	{
	  int t_eff = (t - t0 + length) % length;
	  improved_vector_current[t_eff] += dummy1 * Real(hsum[t]);
	}
	Write(xml_dir, improved_vector_current);
      
	pop(xml_dir);
      }
    }

    pop(xml_dir);
  }

  /* Construct the local vector-current to rho correlators, if desired */

  if ( no_vec_cur == 4 )
  {
    XMLArrayWriter xml_dir(xml,Nd-1);
    push(xml_dir, "Local_vector_current");

    kv = 3*Nd-4;

    for(int k = 0; k < Nd; ++k)
    {
      if( k != j_decay )
      {
	int n = 1 << k;
	kv = kv + 1;
	psi_sq = 0;

	psi_sq = real(trace(adj(anti_quark_prop) * Gamma(n) * quark_prop_1 * Gamma(n)));

	/* Do a slice-wise sum. */
	hsum = sumMulti(psi_sq, phases.getSubset());
//	dummy1 = Real(meson_eta[n]);
	Real dummy1 = 1;

	multi1d<Real> local_vector_current(length);
	for(int t = 0; t < length; ++t)
	{
	  int t_eff = (t - t0 + length) % length;
	  local_vector_current[t_eff] += dummy1 * Real(hsum[t]);
	}
	Write(xml_dir, local_vector_current);
      
	pop(xml_dir);
      }
    }

    pop(xml_dir);
  }


  {
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

    hsum = sumMulti(psi_sq, phases.getSubset());

    multi1d<Real> local_axial_current(length);
    for(int t = 0; t < length; ++t)
    {
      int t_eff = (t - t0 + length) % length;
      local_axial_current[t_eff] -= Real(hsum[t]);
    }
    Write(xml, local_axial_current);


    /* The nonlocal axial current next */

    hsum = sumMulti(chi_sq, phases.getSubset());

    multi1d<Real> nonlocal_axial_current(length);
    for(int t = 0; t < length; ++t)
    {
      int t_eff = (t - t0 + length) % length;
      nonlocal_axial_current[t_eff] += dummy1 * Real(hsum[t]);
    }
    Write(xml, nonlocal_axial_current);
  }

  pop(xml);
              
  END_CODE("curcor2");
}
