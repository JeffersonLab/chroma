// -*- C++ -*-
// $Id: sfcurrents_w.cc,v 3.2 2007-10-16 20:11:17 edwards Exp $
/*! \file
 *  \brief Schroedinger functional application of PCAC
 */

#include "meas/schrfun/sfcurrents_w.h"
#include "util/ferm/transf.h"
#include "actions/ferm/fermbcs/schroedinger_fermbc_w.h"

namespace Chroma 
{
  //! Compute the kprop used in PCAC
  /*! @ingroup schrfun */
  Propagator SFKprop(const LatticePropagator& quark_prop_f,
		     Handle< FermState<LatticeFermion, multi1d<LatticeColorMatrix>,
		     multi1d<LatticeColorMatrix> > > state,
		     const SftMom& phases)
  {
    START_CODE();
  
    QDPIO::cout << __func__ << ": entering" << endl;

    // Decay direction
    int j_decay = phases.getDir();
    int jd = 1 << j_decay;

    // Need to downcast to the appropriate BC
    const SchrFermBC& fermbc = dynamic_cast<const SchrFermBC&>(state->getBC());

    // Location of upper wall source
    int tmax = fermbc.getDecayMax();

    // Grab the links from the state
    const multi1d<LatticeColorMatrix>& u = state->getLinks();

    // Spin projectors
    SpinMatrix P_plus = 0.5*(SpinMatrix(1.0) + (Gamma(jd) * SpinMatrix(1.0)));

    // Common to all currents
    // Sum the forward propagator over time slices to get Kprop
    Propagator kprop = sum(P_plus * (adj(u[j_decay]) * quark_prop_f), phases.getSet()[tmax]);

    QDPIO::cout << __func__ << ": exiting" << endl;

    return kprop;
  }      


  //! Compute Z_V
  /*! @ingroup schrfun */
  void SFCurrentZV(XMLWriter& xml_out, 
		   const string& xml_group,
		   const LatticePropagator& quark_prop_f,
		   const LatticePropagator& quark_prop_b,
		   Handle< SystemSolver<LatticeFermion> > qprop,
		   Handle< FermState<LatticeFermion, multi1d<LatticeColorMatrix>,
		   multi1d<LatticeColorMatrix> > > state,
		   const SftMom& phases)
  {
    START_CODE();
  
    QDPIO::cout << __func__ << ": entering" << endl;

    // Length of lattice in decay direction
    int length = phases.numSubsets();
    int j_decay = phases.getDir();

    // Other useful stuff
    int G5 = Ns*Ns-1;
    int jd = 1 << j_decay;

    // Grab the links from the state
    const multi1d<LatticeColorMatrix>& u = state->getLinks();

    // 3-space volume normalization
    Real norm = 1.0 / Real(QDP::Layout::vol());
    norm *= Real(QDP::Layout::lattSize()[j_decay]);

    // Common to all currents
    // Sum the forward propagator over time slices to get Kprop
    Propagator kprop = SFKprop(quark_prop_f, state, phases);
    kprop *= Real(2)*norm;
      
    /* Construct f_1 */
    Real f_1 = 0.5 * norm2(kprop);
      
    // Construct H'' = H' gamma_5 K, where H' is the propagator
    // from the upper boundary. 
    LatticePropagator quark_prop_bg5k = quark_prop_b * Gamma(G5) * kprop;

    // Construct f_V
    int n = G5 ^ jd;
    LatticeReal r_tmp1 = real(trace(adj(quark_prop_bg5k) * (Gamma(n) * quark_prop_f)));
    multi1d<Double> hrsum = sumMulti(r_tmp1, phases.getSet());

    multi1d<Real> vector_corr(length);
    for(int t = 0; t < length; t++)
      vector_corr[t] = norm * real(hrsum[t]);

    push(xml_out, xml_group);
    write(xml_out, "f_1", f_1);
    write(xml_out, "vector_corr", vector_corr);
    pop(xml_out);

    QDPIO::cout << __func__ << ": exiting" << endl;

    END_CODE();
  }


  //! Compute Z_V
  /*! 
   * @ingroup schrfun 
   *
   * @return number of inverter iterations
   */
  int SFCurrentZA(XMLWriter& xml_out, 
		  const string& xml_group,
		  const multi1d<Real>& pseudo_prop_f,
		  const multi1d<Real>& axial_prop_f,
		  const multi1d<Real>& pseudo_prop_b,
		  const multi1d<Real>& axial_prop_b,
		  const LatticePropagator& quark_prop_f,
		  const LatticePropagator& quark_prop_b,
		  Handle< SystemSolver<LatticeFermion> > qprop,
		  Handle< FermState<LatticeFermion, multi1d<LatticeColorMatrix>,
		  multi1d<LatticeColorMatrix> > > state,
		  const SftMom& phases,
		  int x0, int y0)
  {
    START_CODE();
       
    QDPIO::cout << __func__ << ": entering" << endl;

    // Sanity checks
    if ( x0 < y0 )
    {
      QDPIO::cerr << __func__ << ": Z_A computation requires x0 > y0: x0,y0=" 
		  << x0 << " " << y0 << endl;
      QDP_abort(1);
    }

    // Number of cg iterations
    int ncg_had = 0;

    // Length of lattice in decay direction
    int j_decay = phases.getDir();

    // Other useful stuff
    int G5 = Ns*Ns-1;
    int jd = 1 << j_decay;

    // 3-space volume normalization
    Real norm = 1.0 / Real(QDP::Layout::vol());
    norm *= Real(QDP::Layout::lattSize()[j_decay]);

    // Common to all currents
    // Sum the forward propagator over time slices to get Kprop
    Propagator kprop = SFKprop(quark_prop_f, state, phases);
    kprop *= Real(2)*norm;
      
    /* Construct f_1 */
    Real f_1 = 0.5 * norm2(kprop);
      
    // Construct H'' = H' gamma_5 K, where H' is the propagator
    // from the upper boundary. 
    LatticePropagator quark_prop_bg5k = quark_prop_b * Gamma(G5) * kprop;

    /* Construct the ingredients for f^I_AA */
    LatticeReal faa_tmp = zero;
    LatticeReal fap_tmp = zero;
    LatticeReal fpp_tmp = zero;
    LatticeInteger t_coord = QDP::Layout::latticeCoordinate(j_decay);

    /* "right" A_0 insertion at x */
    LatticeBoolean tmask = (t_coord == x0);
    LatticePropagator quark_prop_s = zero;    // sequential propagator

    for(int color_source = 0; color_source < Nc; ++color_source)
    {
      for(int spin_source = 0; spin_source < Ns; ++spin_source)
      {
	LatticeFermion chi;
	PropToFerm(quark_prop_f, chi, color_source, spin_source);

	int n = jd ^ G5;
	LatticeFermion psi = Gamma(n) * chi;
	/* This gives multiplication with gamma_5 * gamma_0. */
	/* Include the minus sign for multiplication with */
	/* gamma_0 * gamma_5 below. */

	chi = where(tmask, LatticeFermion(-psi), LatticeFermion(zero));
	psi = zero;
	SystemSolverResults_t res = (*qprop)(psi, chi);
	ncg_had += res.n_count;

	FermToProp(psi, quark_prop_s, color_source, spin_source);
      }
    }

    // "left" P insertion at y
    LatticeReal r_tmp1 = -real(trace(adj(quark_prop_bg5k) * quark_prop_s));

    fap_tmp += where(t_coord == (y0+1), r_tmp1, LatticeReal(zero));
    fap_tmp -= where(t_coord == (y0-1), r_tmp1, LatticeReal(zero));

    /* "left" A_0 insertion at y */
    r_tmp1 = real(trace(adj(quark_prop_bg5k) * (Gamma(jd) * quark_prop_s)));
    faa_tmp += where(t_coord == y0, r_tmp1, LatticeReal(zero));

    /* "right" A_0 insertion at y */
    tmask = (t_coord == y0);
    int n = jd ^ G5;
    quark_prop_s = zero;
    for(int color_source = 0; color_source < Nc; ++color_source)
    {
      for(int spin_source = 0; spin_source < Ns; ++spin_source)
      {
	LatticeFermion psi, chi;
	PropToFerm(quark_prop_f, chi, color_source, spin_source);
	psi = Gamma(n) * chi;
	/* This gives multiplication with gamma_5 * gamma_0. */
	/* Include the minus sign for multiplication with */
	/* gamma_0 * gamma_5 below. */

	chi = where(tmask, LatticeFermion(-psi), LatticeFermion(zero));
	psi = zero;
	SystemSolverResults_t res = (*qprop)(psi, chi);
	ncg_had += res.n_count;

	FermToProp(psi, quark_prop_s, color_source, spin_source);
      }
    }

    /* "left" P insertion at x */
    r_tmp1 = real(trace(adj(quark_prop_bg5k) * quark_prop_s));
    fap_tmp += where(t_coord == (x0+1), r_tmp1, LatticeReal(zero));
    fap_tmp -= where(t_coord == (x0-1), r_tmp1, LatticeReal(zero));

    /* "left" A_0 insertion at x */
    r_tmp1 = -real(trace(adj(quark_prop_bg5k) * (Gamma(jd) * quark_prop_s)));
    faa_tmp += where(t_coord == x0, r_tmp1, LatticeReal(zero));

    /* "right" P insertion at x */
    quark_prop_s = zero;
    for(int color_source = 0; color_source < Nc; ++color_source)
    {
      for(int spin_source = 0; spin_source < Ns; ++spin_source)
      {
	LatticeFermion psi, chi;
	PropToFerm(quark_prop_f, chi, color_source, spin_source);
	psi = Gamma(G5) * chi;

	chi  = where(t_coord == (x0+1), psi, LatticeFermion(zero));
	chi -= where(t_coord == (x0-1), psi, LatticeFermion(zero));

	psi = zero;
	SystemSolverResults_t res = (*qprop)(psi, chi);
	ncg_had += res.n_count;

	FermToProp(psi, quark_prop_s, color_source, spin_source);
      }
    }

    /* "left" P insertion at y */
    r_tmp1 = -real(trace(adj(quark_prop_bg5k) * quark_prop_s));
    fpp_tmp += where(t_coord == (y0+1), r_tmp1, LatticeReal(zero));
    fpp_tmp -= where(t_coord == (y0-1), r_tmp1, LatticeReal(zero));

    /* "left" A_0 insertion at y */
    r_tmp1 = real(trace(adj(quark_prop_bg5k) * (Gamma(jd) * quark_prop_s)));
    fap_tmp += where(t_coord == y0, r_tmp1, LatticeReal(zero));

    /* "right" P insertion at y */
    quark_prop_s = zero;
    for(int color_source = 0; color_source < Nc; ++color_source)
    {
      for(int spin_source = 0; spin_source < Ns; ++spin_source)
      {
	LatticeFermion psi, chi;
	PropToFerm(quark_prop_f, chi, color_source, spin_source);
	psi = Gamma(G5) * chi;

	chi  = where(t_coord == (y0+1), psi, LatticeFermion(zero));
	chi -= where(t_coord == (y0-1), psi, LatticeFermion(zero));
	psi = zero;
	SystemSolverResults_t res = (*qprop)(psi, chi);
	ncg_had += res.n_count;

	FermToProp(psi, quark_prop_s, color_source, spin_source);
      }
    }

    /* "left" P insertion at x */
    r_tmp1 = real(trace(adj(quark_prop_bg5k) * quark_prop_s));
    fpp_tmp += where(t_coord == (x0+1), r_tmp1, LatticeReal(zero));
    fpp_tmp -= where(t_coord == (x0-1), r_tmp1, LatticeReal(zero));

    /* "left" A_0 insertion at x */
    r_tmp1 = -real(trace(adj(quark_prop_bg5k) * (Gamma(jd) * quark_prop_s)));
    fap_tmp += where(t_coord == x0, r_tmp1, LatticeReal(zero));

    multi1d<Double> hrsum = sumMulti(faa_tmp, phases.getSet());
    Real f_AA = (real(hrsum[x0]) + real(hrsum[y0])) * norm;

    hrsum = sumMulti(fap_tmp, phases.getSet());
    Real f_AP_PA  = real(hrsum[x0]);
    f_AP_PA += real(hrsum[y0]);
    f_AP_PA += real(hrsum[x0+1]);
    f_AP_PA += real(hrsum[x0-1]);
    f_AP_PA += real(hrsum[y0+1]);
    f_AP_PA += real(hrsum[y0-1]);
    f_AP_PA *= 0.5 * norm;

    hrsum = sumMulti(fpp_tmp, phases.getSet());
    Real f_PP  = real(hrsum[x0+1]);
    f_PP += real(hrsum[x0-1]);
    f_PP += real(hrsum[y0+1]);
    f_PP += real(hrsum[y0-1]);
    f_PP *= 0.25 * norm;

    Real fd_AA;
    fd_AA  = axial_prop_b[y0] * axial_prop_f[x0];
    fd_AA -= axial_prop_b[x0] * axial_prop_f[y0];

    Real fd_AP_PA;
    fd_AP_PA  = pseudo_prop_b[y0+1] * axial_prop_f[x0];
    fd_AP_PA -= pseudo_prop_b[y0-1] * axial_prop_f[x0];
    fd_AP_PA += axial_prop_b[y0] * pseudo_prop_f[x0+1];
    fd_AP_PA -= axial_prop_b[y0] * pseudo_prop_f[x0-1];
    fd_AP_PA -= axial_prop_b[x0] * pseudo_prop_f[y0+1];
    fd_AP_PA += axial_prop_b[x0] * pseudo_prop_f[y0-1];
    fd_AP_PA -= pseudo_prop_b[x0+1] * axial_prop_f[y0];
    fd_AP_PA += pseudo_prop_b[x0-1] * axial_prop_f[y0];
    fd_AP_PA *= 0.5;

    Real fd_PP;
    fd_PP  = pseudo_prop_b[y0+1] * pseudo_prop_f[x0+1];
    fd_PP -= pseudo_prop_b[y0+1] * pseudo_prop_f[x0-1];
    fd_PP -= pseudo_prop_b[y0-1] * pseudo_prop_f[x0+1];
    fd_PP += pseudo_prop_b[y0-1] * pseudo_prop_f[x0-1];
    fd_PP -= pseudo_prop_b[x0+1] * pseudo_prop_f[y0+1];
    fd_PP += pseudo_prop_b[x0+1] * pseudo_prop_f[y0-1];
    fd_PP += pseudo_prop_b[x0-1] * pseudo_prop_f[y0+1];
    fd_PP -= pseudo_prop_b[x0-1] * pseudo_prop_f[y0-1];
    fd_PP *= 0.25;

    push(xml_out, xml_group);
    write(xml_out, "f_1", f_1);
    write(xml_out, "f_AA", f_AA);
    write(xml_out, "f_AP_PA", f_AP_PA);
    write(xml_out, "f_PP", f_PP);
    write(xml_out, "fd_AA", fd_AA);
    write(xml_out, "fd_AP_PA", fd_AP_PA);
    write(xml_out, "fd_PP", fd_PP);
    pop(xml_out);

    QDPIO::cout << __func__ << ": exiting" << endl;

    END_CODE();

    return ncg_had;
  }

} // namespace Chroma
