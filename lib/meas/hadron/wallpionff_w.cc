// $Id: wallpionff_w.cc,v 1.2 2004-01-13 03:57:54 edwards Exp $
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

  int G5 = Ns*Ns-1;
  
  // Project propagator onto zero momentum: Do a slice-wise sum.
  Propagator q_x2 = sum(forw_prop, phases.getSubset()[t_sink]);

  // Start new array group
  XMLArrayWriter xml_array(xml, Nd);
  push(xml_array, "WallPionFormFac");

  // Loop over gamma matrices of the insertion current of insertion current
  for(int mu = 0; mu < Nd; ++mu)
  {
    push(xml_array);
    Write(xml_array, mu);

    int gamma_value = 1 << mu;

    // The local non-conserved vector-current matrix element 
    LatticeComplex corr_local_fn =
      trace(back_prop*Gamma(G5)*q_x2*adj(forw_prop)*Gamma(15)*Gamma(gamma_value));

    multi2d<DComplex> hsum;
//    multi2d<DComplex> hsum_nonlocal;
    hsum = phases.sft(corr_local_fn);

    XMLArrayWriter xml_inser_mom(xml_array, phases.numMom());
    push(xml_inser_mom, "Momenta");

    // Loop over insertion momenta and print out results
    for(int inser_mom_num=0; inser_mom_num<phases.numMom(); ++inser_mom_num) 
    {
      push(xml_inser_mom);
      Write(xml_inser_mom, inser_mom_num);
      write(xml_inser_mom, "inser_mom", phases.numToMom(inser_mom_num)) ;

//      form.formFac[gamma_value].momenta[inser_mom_num].inser_mom = phases.numToMom(inser_mom_num);

      multi1d<Complex> local_cur3ptfn(length); // always compute
//      multi1d<Complex> nonlocal_cur3ptfn;
//      if (compute_nonlocal)
//	nonlocal_cur3ptfn.resize(length);      // possibly compute

      for (int t=0; t < length; ++t) 
      {
        int t_eff = (t - t0 + length) % length;

        local_cur3ptfn[t_eff] = Complex(hsum[inser_mom_num][t]);
//        if (compute_nonlocal)
//          nonlocal_cur3ptfn[t_eff] = 0.5 * Complex(hsum_nonlocal[inser_mom_num][t]);

      } // end for(t)

      // Print out the results
      push(xml_inser_mom, "Wilson_Local_Current_3Pt_fn");
      Write(xml_inser_mom, local_cur3ptfn);
      pop(xml_inser_mom);

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
