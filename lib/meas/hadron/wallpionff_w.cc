// $Id: wallpionff_w.cc,v 1.7 2004-04-07 04:38:36 edwards Exp $
/*! \file
 *  \brief Wall-sink pion form-factors 
 *
 *  Form factors constructed from a quark and a backward quark propagator
 */

#include "chromabase.h"
#include "meas/hadron/wallpionff_w.h"

using namespace QDP;

//! Compute contractions for current insertion 3-point functions.
/*!
 * \ingroup hadron
 *
 * This routine is specific to Wilson fermions!
 *
 * \param xml                xml output ( Modify )
 * \param u                  gauge fields (used for non-local currents) ( Read )
 * \param forw_prop          forward quark propagator ( Read )
 * \param back_prop          backward quark propagator ( Read )
 * \param phases             fourier transform phase factors ( Read )
 * \param t0                 cartesian coordinates of the source ( Read )
 */

void wallPionFormFac(XMLWriter& xml,
		     const multi1d<LatticeColorMatrix>& u, 
		     const LatticePropagator& forw_prop,
		     const LatticePropagator& back_prop, 
		     const SftMom& phases,
		     int t0, int t_sink)
{
  START_CODE("wallPionFormFac");

  // Length of lattice in j_decay direction and 3pt correlations fcns
  int length = phases.numSubsets();

  const int G5 = Ns*Ns-1;
  
  // Project propagator onto zero momentum: Do a slice-wise sum.
  Propagator q_x2 = sum(forw_prop, phases.getSet()[t_sink]);

  // Start new array group
  XMLArrayWriter xml_array(xml, Nd);
  push(xml_array, "WallPionFormFac");

  // Loop over gamma matrices of the insertion current of insertion current
  for(int mu = 0; mu < Nd; ++mu)
  {
    push(xml_array);
    write(xml_array, "mu", mu);

    int gamma_value = 1 << mu;

    // The local non-conserved vector-current matrix element 
    LatticeComplex corr_local_fn =
      trace(back_prop*Gamma(G5)*q_x2*adj(forw_prop)*Gamma(G5)*Gamma(gamma_value));

    multi2d<DComplex> hsum_local = phases.sft(corr_local_fn);

    // Construct the non-local current matrix element 
    //
    // The form of J_mu = (1/2)*[psibar(x+mu)*U^dag_mu*(1+gamma_mu)*psi(x) -
    //                           psibar(x)*U_mu*(1-gamma_mu)*psi(x+mu)]
    // NOTE: the 1/2  is included down below in the sumMulti stuff
    LatticePropagator tmp_prop1 = back_prop * Gamma(15) * q_x2;
    LatticePropagator tmp_prop2 = shift(tmp_prop1, FORWARD, mu);
    LatticeComplex corr_nonlocal_fn =
      trace(adj(u[mu] * shift(forw_prop, FORWARD, mu)) * Gamma(15) * 
	    (tmp_prop1 + Gamma(gamma_value) * tmp_prop1)) -
      trace(adj(forw_prop) * u[mu] * Gamma(15) *
	    (tmp_prop2 - Gamma(gamma_value) * tmp_prop2));

    multi2d<DComplex> hsum_nonlocal = phases.sft(corr_nonlocal_fn);

  
    XMLArrayWriter xml_inser_mom(xml_array, phases.numMom());
    push(xml_inser_mom, "Momenta");

    // Loop over insertion momenta and print out results
    for(int inser_mom_num=0; inser_mom_num<phases.numMom(); ++inser_mom_num) 
    {
      push(xml_inser_mom);
      write(xml_inser_mom, "inser_mom_num", inser_mom_num);
      write(xml_inser_mom, "inser_mom", phases.numToMom(inser_mom_num)) ;

//      form.formFac[gamma_value].momenta[inser_mom_num].inser_mom = phases.numToMom(inser_mom_num);

      multi1d<Complex> local_cur3ptfn(length);
      multi1d<Complex> nonlocal_cur3ptfn(length);

      for (int t=0; t < length; ++t) 
      {
        int t_eff = (t - t0 + length) % length;

        local_cur3ptfn[t_eff] = hsum_local[inser_mom_num][t];
	nonlocal_cur3ptfn[t_eff] = 0.5 * hsum_nonlocal[inser_mom_num][t];
      } // end for(t)

      // Print out the results
      write(xml_inser_mom, "local_cur3ptfn", local_cur3ptfn);
      write(xml_inser_mom, "nonlocal_cur3ptfn", nonlocal_cur3ptfn);

//      form.formFac[gamma_value].momenta[inser_mom_num].local_current    = local_cur3ptfn;
//      form.formFac[gamma_value].momenta[inser_mom_num].nonlocal_current = nonlocal_cur3ptfn;

      pop(xml_inser_mom);  // elem
    } // end for(inser_mom_num)

    pop(xml_inser_mom);    // Momenta
    pop(xml_array);        // elem
  } // end for(gamma_value)
                            
  pop(xml_array);          // WallPionFormFac

  END_CODE("FormFac");
}
