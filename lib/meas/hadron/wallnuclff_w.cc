// $Id: wallnuclff_w.cc,v 1.15 2004-05-04 21:28:50 edwards Exp $
/*! \file
 *  \brief Wall-sink nucleon form-factors 
 *
 *  Form factors constructed from a quark and a backward quark propagator
 */

#include "chromabase.h"
#include "util/ft/sftmom.h"
#include "meas/hadron/wallnuclff_w.h"

using namespace QDP;


//! Compute nonlocal current propagator
/*!
 * \ingroup hadron
 *
 * The form of J_mu = (1/2)*[psibar(x+mu)*U^dag_mu*(1+gamma_mu)*psi(x) -
 *                           psibar(x)*U_mu*(1-gamma_mu)*psi(x+mu)]
 *
 * \param u                  gauge fields ( Read )
 * \param mu                 direction ( Read )
 * \param forw_prop          forward propagator ( Read )
 * \param anti_prop          anti-quark version of forward propagator ( Read )
 *
 * \return nonlocal current propagator
 */
LatticePropagator nonlocalCurrentProp(const multi1d<LatticeColorMatrix>& u, 
				      int mu, 
				      const LatticePropagator& forw_prop,
				      const LatticePropagator& anti_prop)
{
  int gamma_value = 1 << mu;

  LatticePropagator S = shift(anti_prop, FORWARD, mu) * adj(u[mu])
    * (forw_prop + Gamma(gamma_value)*forw_prop)
    - anti_prop * u[mu] * shift(forw_prop - Gamma(gamma_value)*forw_prop, FORWARD, mu);

  return S;
}


//! Compute dbar-d current insertion in nucleon
/*!
 * \ingroup hadron
 *
 * quark contraction within a baryon
 *
 * \param q1        first quark ( Read )
 * \param q2        second quark ( Read )
 * \param q3        third quark ( Read )
 *
 * \return color-contracted spin object
 */
template<class T1, class T2, class T3>
static
LatticeSpinMatrix baryonContract(const T1& q1,
				 const T2& q2, 
				 const T3& q3)
{
  LatticeSpinMatrix  S; 

  S = traceColor(q1 * traceSpin(quarkContract13(q3*Gamma(5), Gamma(5)*q2)))
    + traceColor(q1 * quarkContract13(q3*Gamma(5), Gamma(5)*q2));

  return S;
}


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

  // Length of lattice in j_decay direction and 3pt correlations fcns
  int length = phases.numSubsets();

  multi1d<Complex> local_cur3ptfn(length);
  multi1d<Complex> nonlocal_cur3ptfn(length);

  int G5 = Ns*Ns-1;
  
  // Project propagator onto zero momentum: Do a slice-wise sum.
  Propagator u_x2 = sum(forw_u_prop, phases.getSet()[t_sink]);
  Propagator d_x2 = sum(forw_d_prop, phases.getSet()[t_sink]);
  LatticePropagator anti_u_prop = adj(Gamma(G5)*back_u_prop*Gamma(G5));
  LatticePropagator anti_d_prop = adj(Gamma(G5)*back_d_prop*Gamma(G5));

  // Loop over appropriate form-factor contractions for this system
  XMLArrayWriter xml_seq_src(xml, 4);
  push(xml_seq_src, "FormFac");

  for (int seq_src = 0; seq_src < 4; ++seq_src) 
  {
    push(xml_seq_src);
    write(xml_seq_src, "seq_src", seq_src);

    // Loop over gamma matrices of the insertion current of insertion current
    XMLArrayWriter xml_array(xml_seq_src, Nd);
    push(xml_array, "Insertions");

    for(int mu = 0; mu < Nd; ++mu)
    {
      int gamma_value = 1 << mu;

      push(xml_array);
      write(xml_array, "mu", mu);
      write(xml_array, "gamma_value", gamma_value);

      LatticeComplex corr_local_fn;
      LatticeComplex corr_nonlocal_fn;

      switch (seq_src)
      {
      case 0:
      case 2:
      {

	// "\bar u O u" insertion in proton, ie. "(u C gamma_5 d) u"
	// The local non-conserved current contraction
	LatticePropagator local_insert_prop = anti_u_prop*Gamma(gamma_value)*forw_u_prop;
	LatticeSpinMatrix local_contract = 
	  baryonContract(local_insert_prop, u_x2, d_x2) + baryonContract(u_x2, local_insert_prop, d_x2);

	// Construct the non-local (possibly conserved) current contraction
	LatticePropagator nonlocal_insert_prop = nonlocalCurrentProp(u, mu, forw_u_prop, anti_u_prop);
	LatticeSpinMatrix nonlocal_contract = 
	  baryonContract(nonlocal_insert_prop, u_x2, d_x2) + baryonContract(u_x2, nonlocal_insert_prop, d_x2);
	  
	if (seq_src == 0)
	{
	  /* "\bar u O u" insertion in proton, ie. "(u C gamma_5 d) u" */
	  /* T = (1 + gamma_4) / 2 = (1 + Gamma(8)) / 2 */
	  /* C gamma_5 = Gamma(5) = - (C gamma_5)^T */
	
	  // The local non-conserved vector-current matrix element 
	  corr_local_fn = 0.5 * traceSpin(local_contract + Gamma(8)*local_contract);

	  // The nonlocal (possibly conserved) current matrix element 
	  corr_nonlocal_fn = 0.25 * traceSpin(nonlocal_contract + Gamma(8)*nonlocal_contract);
	}
	else
	{
	  /* "\bar u O u" insertion in proton, ie. "(u C gamma_5 d) u" */
	  /* T = \Sigma_3 (1 + gamma_4) / 2 = -i (Gamma(3) + Gamma(11)) / 2 */
	  /* C gamma_5 = Gamma(5) = - (C gamma_5)^T */

	  // The local non-conserved vector-current matrix element 
	  corr_local_fn = 0.5 * timesMinusI(traceSpin(Gamma(3)*local_contract + Gamma(11)*local_contract));

	  // The nonlocal (possibly conserved) current matrix element 
	  corr_nonlocal_fn = 0.25 * timesMinusI(traceSpin(Gamma(3)*nonlocal_contract + Gamma(11)*nonlocal_contract));
	}
      }
      break;

      case 1:
      case 3:
      {
	// "\bar d O d" insertion in proton, ie. "(u C gamma_5 d) u"
	// The local non-conserved current contraction
	LatticePropagator local_insert_prop = anti_d_prop*Gamma(gamma_value)*forw_d_prop;
	LatticeSpinMatrix local_contract = 
	  baryonContract(u_x2, u_x2, local_insert_prop);

	// Construct the non-local (possibly conserved) current contraction
	LatticePropagator nonlocal_insert_prop = nonlocalCurrentProp(u, mu, forw_d_prop, anti_d_prop);
	LatticeSpinMatrix nonlocal_contract = 
	  baryonContract(u_x2, u_x2, nonlocal_insert_prop);

	if (seq_src == 1)
	{
	  /* "\bar d O d" insertion in proton, ie. "(u C gamma_5 d) u" */
	  /* T = (1 + gamma_4) / 2 = (1 + Gamma(8)) / 2 */

	  // The local non-conserved vector-current matrix element 
	  corr_local_fn = 0.5 * traceSpin(local_contract + Gamma(8)*local_contract);

	  // The nonlocal (possibly conserved) current matrix element 
	  corr_nonlocal_fn = 0.25 * traceSpin(nonlocal_contract + Gamma(8)*nonlocal_contract);
	}
	else
	{
	  /* "\bar d O d" insertion in proton, ie. "(u C gamma_5 d) u" */
	  /* T = \Sigma_3 (1 + gamma_4) / 2 = -i (Gamma(3) + Gamma(11)) / 2 */
	  /* C gamma_5 = Gamma(5) = - (C gamma_5)^T */

	  // The local non-conserved vector-current matrix element 
	  corr_local_fn = 0.5 * timesMinusI(traceSpin(Gamma(3)*local_contract + Gamma(11)*local_contract));

	  // The nonlocal (possibly conserved) current matrix element 
	  corr_nonlocal_fn = 0.25 * timesMinusI(traceSpin(Gamma(3)*nonlocal_contract + Gamma(11)*nonlocal_contract));
	}
      }
      break;

      default:
	QDP_error_exit("Unknown sequential source type", seq_src);
      }

      multi2d<DComplex> hsum_local = phases.sft(corr_local_fn);

      multi2d<DComplex> hsum_nonlocal = phases.sft(corr_nonlocal_fn);
  
      XMLArrayWriter xml_inser_mom(xml_array, phases.numMom());
      push(xml_inser_mom, "Momenta");

      // Loop over insertion momenta and print out results
      for(int inser_mom_num=0; inser_mom_num<phases.numMom(); ++inser_mom_num) 
      {
	push(xml_inser_mom);
	write(xml_inser_mom, "inser_mom_num", inser_mom_num);
	write(xml_inser_mom, "inser_mom", phases.numToMom(inser_mom_num)) ;

	for (int t=0; t < length; ++t) 
	{
	  int t_eff = (t - t0 + length) % length;

	  local_cur3ptfn[t_eff] = Complex(hsum_local[inser_mom_num][t]);
	  nonlocal_cur3ptfn[t_eff] = Complex(hsum_nonlocal[inser_mom_num][t]);
	} // end for(t)

	// Print out the results
	write(xml_inser_mom, "local_cur3ptfn", local_cur3ptfn);
	write(xml_inser_mom, "nonlocal_cur3ptfn", nonlocal_cur3ptfn);

	pop(xml_inser_mom);  // elem
      } // end for(inser_mom_num)

      pop(xml_inser_mom);    // Momenta
      pop(xml_array);        // elem
    } // end for(gamma_value)
                            
    pop(xml_array);          // FormFac
    pop(xml_seq_src);        // elem
  } // end for(seq_src)
                            
  pop(xml_seq_src);          // WallNuclFormFac


  END_CODE("wallNuclFormFac");
}
