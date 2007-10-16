// -*- C++ -*-
// $Id: sfpcac_w.cc,v 3.14 2007-10-16 20:12:11 edwards Exp $
/*! \file
 *  \brief Schroedinger functional application of PCAC
 */

#include "meas/schrfun/sfpcac_w.h"
#include "meas/schrfun/sfcorr_w.h"
#include "meas/schrfun/sfcurrents_w.h"
#include "meas/sources/walfil_w.h"
#include "util/ferm/transf.h"

#include "actions/ferm/fermbcs/schroedinger_fermbc_w.h"

namespace Chroma 
{
  //! Schroedinger functional stuff
  /*!
   * @ingroup schrfun
   *
   * Compute correlation functions between axial current or pseudescalar
   * density and boundary fields using Schroedinger BC.
   *
   * Also computed, on demand, are correlation functions between both
   * boundaries with zero, one (vector current) and two (axial current or
   * pseudoscalar density) insertions in the bulk. These currents are
   * controlled by the ZVfactP and ZAfactP boolean flags.
   *
   * Compute quark propagators by using the qprop SystemSolver.
   * The initial guess for the inverter is zero.
   *
   * The results are written to the xml file.
   *
   * For further details see the comments in the dependent subroutines.
   *
   * \param state         gauge field state ( Read )
   * \param qprop         propagator solver ( Read )
   * \param phases        object holds list of momenta and Fourier phases ( Read )
   * \param ZVfactP       flag for doing Z_V measurements ( Read )
   * \param ZAfactP       flag for doing Z_A measurements ( Read )
   * \param x0            time slices with axial current insertions ( Read ) 
   * \param y0            time slices with axial current insertions ( Read ) 
   * \param xml           xml file object ( Write )
   * \param xml_group     string used for writing xml data ( Read )
   */
  void SFpcac(Handle< SystemSolver<LatticeFermion> > qprop,
	      Handle< FermState<LatticeFermion, multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> > > state,
	      const SftMom& phases,
	      bool ZVfactP, bool ZAfactP, 
	      int x0, int y0,
	      XMLWriter& xml_out,
	      const string& xml_group)
  {
    START_CODE();
  
    QDPIO::cout << __func__ << ": entering" << endl;

    if ( Ns != 4 )
    {
      QDPIO::cerr << __func__ << ": only supports 4 spin components" << endl;
      QDP_abort(1);
    }
  
    // Need to downcast to the appropriate BC
    const SchrFermBC& fermbc = dynamic_cast<const SchrFermBC&>(state->getBC());

    // Outside group
    push(xml_out, xml_group);

    // Length of lattice in decay direction
    int length = phases.numSubsets();
    int j_decay = phases.getDir();

    // Other useful stuff
    int G5 = Ns*Ns-1;
    int jd = 1 << j_decay;

    multi1d<Real> pseudo_prop_f(length);
    multi1d<Real> axial_prop_f(length);
    multi1d<Real> pseudo_prop_b(length);
    multi1d<Real> axial_prop_b(length);

    // 3-space volume normalization
    Real norm = 1.0 / Real(QDP::Layout::vol());
    norm *= Real(QDP::Layout::lattSize()[j_decay]);

    // Spin projectors
    SpinMatrix g_one = 1.0;
    SpinMatrix P_plus  = 0.5*(g_one + (Gamma(jd) * g_one));
    SpinMatrix P_minus = 0.5*(g_one - (Gamma(jd) * g_one));

    /* Location of upper wall source */
    int tmin = fermbc.getDecayMin();
    int tmax = fermbc.getDecayMax();

    // Grab the links from the state
    const multi1d<LatticeColorMatrix>& u = state->getLinks();

    // Total number of inversions
    int ncg_had = 0;

    LatticePropagator quark_prop_f = zero;
    LatticePropagator quark_prop_b = zero;

    // Temporaries
    multi1d<Real> pseudo_prop_tmp;
    multi1d<Real> axial_prop_tmp;

    push(xml_out, "PCAC_measurements");
    for(int direction = -1; direction <= 1; direction+=2)
    {
      int t0 = (direction == -1) ? tmax : tmin;

      for(int color_source = 0; color_source < Nc; ++color_source)
      {
	for(int spin_source = 0; spin_source < Ns; ++spin_source)
	{
	  /* Compute quark propagator "psi" using source "chi" with type specified */
	  /* by WALL_SOURCE, and color and spin equal to color_source and spin_source. */
	  LatticeFermion psi, chi;
	  {
	    LatticeFermion tmp1;
	    walfil(tmp1, t0, j_decay, color_source, spin_source);

	    if (direction == -1)
	    {
	      chi = P_minus * (u[j_decay] * tmp1);
	    }
	    else
	    {
	      chi = shift(P_plus * (adj(u[j_decay]) * tmp1), BACKWARD, j_decay);
	    }
  	  
	    // Solve for the propagator
	    psi = zero;
	    SystemSolverResults_t res = (*qprop)(psi, chi);
	    ncg_had += res.n_count;
	  }


	  /* Store in quark_prop_f/b */
	  if (direction == -1)
	  {
	    FermToProp(psi, quark_prop_b, color_source, spin_source);
	  }
	  else
	  {
	    FermToProp(psi, quark_prop_f, color_source, spin_source);
	  }

	}
      }

      /* Time reverse backwards source */
      if (direction == -1)
      {
	// Construct the pion axial current divergence and the pion correlator
	SFcorr(pseudo_prop_b, axial_prop_b, quark_prop_b, phases);

	// Normalize to compare to everybody else
	pseudo_prop_b *= norm;
	axial_prop_b  *= norm;

	pseudo_prop_tmp.resize(length);
	axial_prop_tmp.resize(length);

	// Time reverse
	for(int t = 0; t < (length-1)/2 + 1; ++t)
	{
	  int t_eff = length - t - 1;

	  pseudo_prop_tmp[t]     = pseudo_prop_b[t_eff];
	  pseudo_prop_tmp[t_eff] = pseudo_prop_b[t];

	  axial_prop_tmp[t]      = -axial_prop_b[t_eff];
	  axial_prop_tmp[t_eff]  = -axial_prop_b[t];
	}

      }
      else
      {
	// Construct the pion axial current divergence and the pion correlator
	SFcorr(pseudo_prop_f, axial_prop_f, quark_prop_f, phases);

	// Normalize to compare to everybody else
	pseudo_prop_f *= norm;
	axial_prop_f  *= norm;

	pseudo_prop_tmp = pseudo_prop_f;
	axial_prop_tmp  = axial_prop_f;
      }

      // Write out results
      push(xml_out, "elem");
      write(xml_out, "direction", direction);
      write(xml_out, "pseudo_prop", pseudo_prop_tmp);
      write(xml_out, "axial_prop", axial_prop_tmp);
      pop(xml_out);

    }  // end for direction
    pop(xml_out);  // PCAC_measurements

#if 0
    {
      multi1d<Double> prop_corr = sumMulti(localNorm2(quark_prop_f), 
					   phases.getSet());

      push(xml_out, "Forward_prop_test");
      write(xml_out, "quark_prop_f", prop_corr);
      write(xml_out, "pseudo_prop_f", pseudo_prop_f);
      pop(xml_out);
    }
#endif

    //
    // Currents
    //
    // Compute Z_V
    if (ZVfactP)
      SFCurrentZV(xml_out, "ZV_measurements",
		  quark_prop_f, quark_prop_b, qprop, state, phases);

    // Compute Z_A
    if (ZAfactP)
      ncg_had += SFCurrentZA(xml_out, "ZA_measurements", 
			     pseudo_prop_f, axial_prop_f, pseudo_prop_b, axial_prop_b,
			     quark_prop_f, quark_prop_b, qprop, state, phases, x0, y0);
  
    QDPIO::cout << __func__ << ": print iterations" << endl;

    push(xml_out,"Relaxation_Iterations");
    write(xml_out, "ncg_had", ncg_had);
    pop(xml_out);

    pop(xml_out);   // xml_group

    QDPIO::cout << __func__ << ": exiting" << endl;

    END_CODE();
  }

}
