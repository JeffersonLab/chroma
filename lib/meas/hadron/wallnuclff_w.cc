// $Id: wallnuclff_w.cc,v 1.5 2004-04-14 04:42:06 edwards Exp $
/*! \file
 *  \brief Wall-sink nucleon form-factors 
 *
 *  Form factors constructed from a quark and a backward quark propagator
 */

#include "chromabase.h"
#include "util/ft/sftmom.h"
#include "meas/hadron/wallnuclff_w.h"

using namespace QDP;

/*
 * Structures for hadron parts
 *
 * \ingroup hadron
 *
 * @{
 */


/*! @} */  // end of group hadron


//! Compute contractions for current insertion 3-point functions.
/*!
 * \ingroup hadron
 *
 * This routine is specific to Wilson fermions!
 *
 * \param xml                buffer for writing the data ( Write )
 * \param u                  gauge fields (used for non-local currents) ( Read )
 * \param forw_u_prop        forward U quark propagator ( Read )
 * \param back_u_prop        backward D quark propagator ( Read )
 * \param forw_d_prop        forward U quark propagator ( Read )
 * \param back_d_prop        backward D quark propagator ( Read )
 * \param phases             fourier transform phase factors ( Read )
 * \param t0                 time coordinates of the source ( Read )
 * \param t_sink             time coordinates of the sink ( Read )
 */

void wallNuclFormFac(XMLWriter& xml,
		     const multi1d<LatticeColorMatrix>& u, 
		     const LatticePropagator& forw_u_prop,
		     const LatticePropagator& back_u_prop, 
		     const LatticePropagator& forw_d_prop,
		     const LatticePropagator& back_d_prop, 
		     const SftMom& phases,
		     int t0, int t_sink)
{
  START_CODE("wallNuclFormFac");

  // Start new array group
  XMLArrayWriter xml_array(xml, Nd);
  push(xml_array, "WallNuclFormFac");

  // Length of lattice in j_decay direction and 3pt correlations fcns
  int length = phases.numSubsets();

  int G5 = Ns*Ns-1;
  
  Real e_u = 2.0/3.0;
  Real e_d = -1.0/3.0;

  LatticePropagator q1_tmp;

  // Project propagator onto zero momentum: Do a slice-wise sum.
  Propagator u_x2 = sum(forw_u_prop, phases.getSet()[t_sink]);
  Propagator d_x2 = sum(forw_d_prop, phases.getSet()[t_sink]);

//  form.formFac.resize(Nd*Nd);

  // Loop over gamma matrices of the insertion current of insertion current
  for(int gamma_value = 0; gamma_value < Nd*Nd; ++gamma_value)
  {
    push(xml_array);
    write(xml_array, "gamma_value", gamma_value);

    //  For the case where the gamma value indicates we are evaluating either
    //  the vector or axial vector currents, we will also evaluate
    //  the non-local currents.  The non-local vector current is the conserved
    //  current.  The non-local axial vector current would be partially
    //  conserved but for the Wilson term.  In these cases we will set
    //  mu = corresponding direction.  In all other cases, we will set mu = -1.

    bool compute_nonlocal;
    int mu;

    switch(gamma_value){
    case  1:
    case 14:
      mu = 0;
      compute_nonlocal = true;
      break;
    case  2:
    case 13:
      mu = 1;
      compute_nonlocal = true;
      break;
    case  4:
    case 11:
      mu = 2;
      compute_nonlocal = true;
      break;
    case  8:
    case  7:
      mu = 3;
      compute_nonlocal = true;
      break;
    default:
      mu = -1;
      compute_nonlocal = false;
    }

    /* "\bar u O u" insertion in proton, ie. "(u C gamma_5 d) u" */
    /* T = (1 + gamma_4) / 2 = (1 + Gamma(8)) / 2 */
    /* C gamma_5 = Gamma(5) = - (C gamma_5)^T */

    q1_tmp = 1;

    // The local non-conserved vector-current matrix element 
    LatticeComplex corr_local_fn;
    LatticePropagator anti_u_prop = Gamma(15)*back_u_prop*Gamma(15);
    LatticePropagator anti_d_prop = Gamma(15)*back_d_prop*Gamma(15);

    // Term 1
    corr_local_fn = -e_u*trace(anti_u_prop*Gamma(gamma_value)*forw_u_prop*Gamma(5)*
			       quarkContract13(Gamma(5)*d_x2, u_x2*(q1_tmp+Gamma(8)*q1_tmp)));

    // Term 2
    corr_local_fn += -e_u*trace(traceSpin(u_x2+Gamma(8)*u_x2) *
				quarkContract13(anti_u_prop*Gamma(gamma_value)*forw_u_prop*Gamma(5),
						Gamma(5)*d_x2));

    // Term 3
    corr_local_fn += e_u*trace(traceSpin(anti_u_prop*Gamma(gamma_value)*(forw_u_prop+forw_u_prop*Gamma(8)))*
			       quarkContract13(Gamma(5)*d_x2, u_x2*Gamma(5)));

    // Term 4
    corr_local_fn += e_u*trace(anti_u_prop*Gamma(gamma_value)*(forw_u_prop+forw_u_prop*Gamma(8))*
			       quarkContract13(u_x2*Gamma(5), Gamma(5)*d_x2));

    // Term 5
    corr_local_fn += e_d*trace(Gamma(5)*anti_d_prop*Gamma(gamma_value)*forw_d_prop * 
			       quarkContract14(u_x2*Gamma(5), u_x2+u_x2*Gamma(8)));

    // Term 6
    corr_local_fn += e_d*trace(traceSpin(u_x2+Gamma(8)*u_x2) * 
			       quarkContract14(u_x2*Gamma(5), 
					       Gamma(5)*anti_d_prop*Gamma(gamma_value)*forw_d_prop));
    corr_local_fn *= 0.5;


    multi2d<DComplex> hsum, hsum_nonlocal;
    hsum = phases.sft(corr_local_fn);

#if 0
    // Construct the non-local current matrix element 
    //
    // The form of J_mu = (1/2)*[psibar(x+mu)*U^dag_mu*(1+gamma_mu)*psi(x) -
    //                           psibar(x)*U_mu*(1-gamma_mu)*psi(x+mu)]
    // NOTE: the 1/2  is included down below in the sumMulti stuff
    LatticeComplex corr_nonlocal_fn;
    if(compute_nonlocal){
      corr_nonlocal_fn =
        trace(adj(u[mu] * shift(anti_quark_prop, FORWARD, mu)) *
          (quark_propagator + Gamma(gamma_value) * quark_propagator));
      LatticePropagator tmp_prop1 = u[mu] *
        shift(quark_propagator, FORWARD, mu);
      corr_nonlocal_fn -= trace(adj(anti_quark_prop) *
                            (tmp_prop1 - Gamma(gamma_value) * tmp_prop1));

      hsum_nonlocal = phases.sft(corr_nonlocal_fn);
    }
#endif
  
//    form.formFac[gamma_value].gamma_value = gamma_value;
//    form.formFac[gamma_value].momenta.resize(phases.numMom());  // hold momenta output
    
    XMLArrayWriter xml_inser_mom(xml_array, phases.numMom());
    push(xml_inser_mom, "Momenta");

    // Loop over insertion momenta and print out results
    for(int inser_mom_num=0; inser_mom_num<phases.numMom(); ++inser_mom_num) 
    {
      push(xml_inser_mom);
      write(xml_inser_mom, "inser_mom_num", inser_mom_num);
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
      write(xml_inser_mom, "local_cur3ptfn", local_cur3ptfn);
//o      write(xml_inser_mom, "nonlocal_cur3ptfn", nonlocal_cur3ptfn);

//      form.formFac[gamma_value].momenta[inser_mom_num].local_current    = local_cur3ptfn;
//      form.formFac[gamma_value].momenta[inser_mom_num].nonlocal_current = nonlocal_cur3ptfn;

      pop(xml_inser_mom);  // elem
    } // end for(inser_mom_num)

    pop(xml_inser_mom);    // Momenta
    pop(xml_array);        // elem
  } // end for(gamma_value)
                            
  pop(xml_array);          // WallNuclFormFac

  END_CODE("wallNuclFormFac");
}
