// $Id: wallnuclff_w.cc,v 1.8 2004-04-18 20:39:33 edwards Exp $
/*! \file
 *  \brief Wall-sink nucleon form-factors 
 *
 *  Form factors constructed from a quark and a backward quark propagator
 */

#include "chromabase.h"
#include "util/ft/sftmom.h"
#include "meas/hadron/wallnuclff_w.h"

using namespace QDP;


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
  LatticePropagator anti_u_prop = Gamma(G5)*back_u_prop*Gamma(G5);
  LatticePropagator anti_d_prop = Gamma(G5)*back_d_prop*Gamma(G5);

  // Loop over appropriate form-factor contractions for this system
  XMLArrayWriter xml_seq_src(xml, 4);
  push(xml_seq_src, "FormFac");

  for (int seq_src = 0; seq_src < 4; ++seq_src) 
  {
    push(xml_seq_src);
    write(xml_seq_src, "seq_src", seq_src);

    QDPIO::cout << "WallNuclFormFac: seq_src " << seq_src << endl;

    // Loop over gamma matrices of the insertion current of insertion current
    XMLArrayWriter xml_array(xml_seq_src, Nd);
    push(xml_array, "Insertions");

    for(int mu = 0; mu < Nd; ++mu)
    {
      int gamma_value = 1 << mu;

      push(xml_array);
      write(xml_array, "mu", mu);
      write(xml_array, "gamma_value", gamma_value);

      QDPIO::cout << "WallNuclFormFac: gamma_value " << gamma_value << endl;
    
      LatticeComplex corr_local_fn;
      LatticeComplex corr_nonlocal_fn;

      // For conserved current
      LatticePropagator tmp_prop1 = u[mu] * shift(forw_u_prop, FORWARD, mu);
      LatticePropagator tmp_prop2 = adj(u[mu] * shift(anti_u_prop, FORWARD, mu))
	* (forw_u_prop + Gamma(gamma_value)*forw_u_prop)
	- adj(anti_u_prop) * (tmp_prop1 - Gamma(gamma_value)*tmp_prop1);

      switch (seq_src)
      {
      case 0:
      {
	/* "\bar u O u" insertion in proton, ie. "(u C gamma_5 d) u" */
	/* T = (1 + gamma_4) / 2 = (1 + Gamma(8)) / 2 */
	/* C gamma_5 = Gamma(5) = - (C gamma_5)^T */
	
	// The local non-conserved vector-current matrix element 
	// Term 1
	corr_local_fn = -trace(anti_u_prop*Gamma(gamma_value)*forw_u_prop*Gamma(5)*
			       quarkContract13(Gamma(5)*d_x2, u_x2+u_x2*Gamma(8)));
	// Term 2
	corr_local_fn += -trace(traceSpin(u_x2+Gamma(8)*u_x2) *
				quarkContract13(anti_u_prop*Gamma(gamma_value)*forw_u_prop*Gamma(5),
						Gamma(5)*d_x2));
	// Term 3
	corr_local_fn += trace(traceSpin(anti_u_prop*Gamma(gamma_value)*(forw_u_prop+forw_u_prop*Gamma(8)))*
			       quarkContract13(Gamma(5)*d_x2, u_x2*Gamma(5)));
	// Term 4
	corr_local_fn += trace(anti_u_prop*Gamma(gamma_value)*(forw_u_prop+forw_u_prop*Gamma(8))*
			       quarkContract13(u_x2*Gamma(5), Gamma(5)*d_x2));

	/*
	 * Construct the non-local current matrix element 
	 *
	 * The form of J_mu = (1/2)*[psibar(x+mu)*U^dag_mu*(1+gamma_mu)*psi(x) -
	 *                           psibar(x)*U_mu*(1-gamma_mu)*psi(x+mu)]
	 */
	// Term 1
	corr_nonlocal_fn = -trace(tmp_prop2 * Gamma(5) * 
				  quarkContract13(Gamma(5)*d_x2, u_x2+u_x2*Gamma(8)));
	// Term 2
	corr_nonlocal_fn += -trace(traceSpin(u_x2+Gamma(8)*u_x2) *
				   quarkContract13(tmp_prop2*Gamma(5), 
						Gamma(5)*d_x2));
	// Term 3
	corr_nonlocal_fn += trace(traceSpin(tmp_prop2 + tmp_prop2*Gamma(8))*
				  quarkContract13(Gamma(5)*d_x2, u_x2*Gamma(5)));
	// Term 4
	corr_nonlocal_fn += trace((tmp_prop2 + tmp_prop2*Gamma(8))*
				  quarkContract13(u_x2*Gamma(5), Gamma(5)*d_x2));
      }
      break;
      
      case 1:
      {
	/* "\bar d O d" insertion in proton, ie. "(u C gamma_5 d) u" */
	/* T = (1 + gamma_4) / 2 = (1 + Gamma(8)) / 2 */

	// The local non-conserved vector-current matrix element 
	// Term 5
	corr_local_fn = trace(Gamma(5)*anti_d_prop*Gamma(gamma_value)*forw_d_prop * 
			      quarkContract14(u_x2*Gamma(5), u_x2+u_x2*Gamma(8)));
	// Term 6
	corr_local_fn += trace(traceSpin(u_x2+Gamma(8)*u_x2) * 
			       quarkContract13(u_x2*Gamma(5), 
					       Gamma(5)*anti_d_prop*Gamma(gamma_value)*forw_d_prop));

	// Construct the non-local current matrix element 
	// Term 5
	corr_local_fn = trace(Gamma(5)*tmp_prop2*
			      quarkContract14(u_x2*Gamma(5), u_x2+u_x2*Gamma(8)));
	// Term 6
	corr_local_fn += trace(traceSpin(u_x2+Gamma(8)*u_x2) * 
			       quarkContract13(u_x2*Gamma(5), Gamma(5)*tmp_prop2));
      }
      break;
	
      case 2:
      {
	/* "\bar u O u" insertion in proton, ie. "(u C gamma_5 d) u" */
	/* T = \Sigma_3 (1 + gamma_4) / 2 = -i (Gamma(3) + Gamma(11)) / 2 */
	/* C gamma_5 = Gamma(5) = - (C gamma_5)^T */

	// The local non-conserved vector-current matrix element 
	// NOTE: extract the common  "-i" piece
	LatticeComplex local_tmp;
	// Term 1
	local_tmp = -trace(anti_u_prop*Gamma(gamma_value)*forw_u_prop*Gamma(5)*
			   quarkContract13(Gamma(5)*d_x2, u_x2*Gamma(3)+u_x2*Gamma(11)));
	// Term 2
	local_tmp += -trace(traceSpin(Gamma(3)*u_x2+Gamma(11)*u_x2) *
			    quarkContract13(anti_u_prop*Gamma(gamma_value)*forw_u_prop*Gamma(5),
					    Gamma(5)*d_x2));
	// Term 3
	local_tmp += trace(traceSpin(anti_u_prop*Gamma(gamma_value)*(forw_u_prop*Gamma(3)+forw_u_prop*Gamma(11)))*
			   quarkContract13(Gamma(5)*d_x2, u_x2*Gamma(5)));
	// Term 4
	local_tmp += trace(anti_u_prop*Gamma(gamma_value)*(forw_u_prop*Gamma(3)+forw_u_prop*Gamma(11))*
			   quarkContract13(u_x2*Gamma(5), Gamma(5)*d_x2));
	corr_local_fn = timesMinusI(local_tmp);

	/*
	 * Construct the non-local current matrix element 
	 *
	 * The form of J_mu = (1/2)*[psibar(x+mu)*U^dag_mu*(1+gamma_mu)*psi(x) -
	 *                           psibar(x)*U_mu*(1-gamma_mu)*psi(x+mu)]
	 */
	LatticeComplex nonlocal_tmp;
	// Term 1
	nonlocal_tmp = -trace(tmp_prop2 * Gamma(5) * 
			      quarkContract13(Gamma(5)*d_x2, u_x2*Gamma(3)+u_x2*Gamma(11)));
	// Term 2
	nonlocal_tmp += -trace(traceSpin(Gamma(3)*u_x2+Gamma(11)*u_x2) *
			       quarkContract13(tmp_prop2*Gamma(5), 
					       Gamma(5)*d_x2));
	// Term 3
	nonlocal_tmp += trace(traceSpin(tmp_prop2*Gamma(3) + tmp_prop2*Gamma(11))*
			      quarkContract13(Gamma(5)*d_x2, u_x2*Gamma(5)));
	// Term 4
	nonlocal_tmp += trace((tmp_prop2*Gamma(3) + tmp_prop2*Gamma(11))*
			      quarkContract13(u_x2*Gamma(5), Gamma(5)*d_x2));
	corr_nonlocal_fn = timesMinusI(nonlocal_tmp);
      }
      break;

      case 3:
      {
	/* "\bar d O d" insertion in proton, ie. "(u C gamma_5 d) u" */
	/* T = \Sigma_3 (1 + gamma_4) / 2 = -i (Gamma(3) + Gamma(11)) / 2 */
	/* C gamma_5 = Gamma(5) = - (C gamma_5)^T */

	// The local non-conserved vector-current matrix element 
	// NOTE: extract the common  "-i" piece
	LatticeComplex local_tmp;
	// Term 5
	local_tmp = trace(Gamma(5)*anti_d_prop*Gamma(gamma_value)*forw_d_prop * 
			  quarkContract14(u_x2*Gamma(5), u_x2*Gamma(3)+u_x2*Gamma(11)));
	// Term 6
	local_tmp += trace(traceSpin(Gamma(3)*u_x2+Gamma(11)*u_x2) * 
			   quarkContract13(u_x2*Gamma(5), 
					   Gamma(5)*anti_d_prop*Gamma(gamma_value)*forw_d_prop));
	corr_local_fn = timesMinusI(local_tmp); 

	// Construct the non-local current matrix element 
	// NOTE: extract the common  "-i" piece
	LatticeComplex nonlocal_tmp;
	// Term 5
	nonlocal_tmp = trace(Gamma(5)*tmp_prop2*
			     quarkContract14(u_x2*Gamma(5), u_x2*Gamma(3)+u_x2*Gamma(11)));
	// Term 6
	nonlocal_tmp += trace(traceSpin(Gamma(3)*u_x2+Gamma(11)*u_x2) * 
			      quarkContract13(u_x2*Gamma(5), Gamma(5)*tmp_prop2));
	corr_nonlocal_fn = timesMinusI(nonlocal_tmp);
      }
      break;

      default:
	QDP_error_exit("Unknown sequential source type", seq_src);
      }

      QDPIO::cout << "WallNuclFormFac: here" << endl;
    
      corr_local_fn *= 0.5;
      multi2d<DComplex> hsum_local = phases.sft(corr_local_fn);

      corr_local_fn *= 0.25;
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
