// $Id: wallrhopiff_w.cc,v 1.1 2004-06-04 04:02:57 edwards Exp $
/*! \file
 *  \brief Wall-sink rho^+ -> gamma+pi^+ form-factors 
 *
 *  Form factors constructed from a quark and a backward quark propagator
 */

#include "chromabase.h"
#include "util/ft/sftmom.h"
#include "meas/hadron/wallff_w.h"
#include "meas/hadron/walldeltapff_w.h"

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

void wallRhoPiFormFac(XMLWriter& xml,
		       const multi1d<LatticeColorMatrix>& u, 
		       const LatticePropagator& forw_u_prop,
		       const LatticePropagator& back_u_prop, 
		       const LatticePropagator& forw_d_prop,
		       const LatticePropagator& back_d_prop, 
		       const SftMom& phases,
		       int t0, int t_sink)
{
  START_CODE("wallRhoPiFormFac");

  if ( Ns != 4 || Nc != 3 || Nd != 4 )	// Code is specific to Ns=4, Nc=3, Nd=4
    return;

  // Length of lattice in j_decay direction and 3pt correlations fcns
  int length = phases.numSubsets();

  // Mega-structure holding form-factors
  WallFormFac_formfacs_t form;

  int G5 = Ns*Ns-1;
  
  // Spin projectors
  multi1d<SpinMatrix> S_proj(Nd);
  SpinMatrix  g_one = 1.0;

  // T = (1 + gamma_4) / 2 = (1 + Gamma(8)) / 2
  S_proj[0] = 0.5 * (g_one + (g_one * Gamma(8)));

  // T = \Sigma_k (1 + gamma_4) / 2 = -i gamma_i gamma_j S_proj[0] i,j cyclic perms
  S_proj[1] = timesMinusI(Gamma(1 << 1) * (Gamma(1 << 2) * S_proj[0]));
  S_proj[2] = timesMinusI(Gamma(1 << 2) * (Gamma(1 << 0) * S_proj[0]));
  S_proj[3] = timesMinusI(Gamma(1 << 0) * (Gamma(1 << 1) * S_proj[0]));


  // Project propagator onto zero momentum: Do a slice-wise sum.
  Propagator u_x2 = sum(forw_u_prop, phases.getSet()[t_sink]);
  Propagator d_x2 = sum(forw_d_prop, phases.getSet()[t_sink]);
  LatticePropagator anti_u_prop = adj(Gamma(G5)*back_u_prop*Gamma(G5));
  LatticePropagator anti_d_prop = adj(Gamma(G5)*back_d_prop*Gamma(G5));


  // Resize some things
  form.formFacs.resize(2*Nd*2);
  for (int seq_src = 0; seq_src < form.formFacs.size(); ++seq_src) 
    form.formFacs[seq_src].insertions.resize(Nd);


  // For calculational purpose, loop over insertions first.
  // This is out-of-order from storage within the data structure
  for(int mu = 0; mu < Nd; ++mu)
  {
    int gamma_value = 1 << mu;

    // !!!!!!!!!!!!!!!! HACK !!!!!!!!!!!!!!!
    int sigma = 0;    // NEED TO FIX THIS!!!!!
    // !!!!!!!!!!!!!!!! END-HACK !!!!!!!!!!!!!!!

    multi1d<LatticeSpinMatrix> local_contract(2);
    multi1d<LatticeSpinMatrix> nonlocal_contract(2);

    // Loop over "u"=0 or "d"=1 pieces
    for(int ud = 0; ud < 2; ++ud)
    {
      switch (ud)
      {
      case 0:
      {
	// "\bar u O u" insertion in delta and proton, ie. "(u C gamma_5 d) u"
	// The local non-conserved current contraction
	LatticePropagator local_insert_prop = 
	  anti_u_prop*Gamma(gamma_value)*forw_u_prop;

	local_contract[0] = 
	  rhopiContract(local_insert_prop, u_x2, d_x2, sigma) + 
	  rhopiContract(u_x2, local_insert_prop, d_x2, sigma);

	local_contract[1] = 
	  pDeltaContract(local_insert_prop, u_x2, d_x2, sigma) + 
	  pDeltaContract(u_x2, local_insert_prop, d_x2, sigma);

	// Construct the non-local (possibly conserved) current contraction
	LatticePropagator nonlocal_insert_prop = 
	  nonlocalCurrentProp(u, mu, forw_u_prop, anti_u_prop);

	nonlocal_contract[0] = 
	  rhopiContract(nonlocal_insert_prop, u_x2, d_x2, sigma) + 
	  rhopiContract(u_x2, nonlocal_insert_prop, d_x2, sigma);

	nonlocal_contract[1] = 
	  pDeltaContract(nonlocal_insert_prop, u_x2, d_x2, sigma) + 
	  pDeltaContract(u_x2, nonlocal_insert_prop, d_x2, sigma);
      }
      break;

      case 1:
      {
	// "\bar d O d" insertion in delta and proton, ie. "(u C gamma_5 d) u"
	// The local non-conserved current contraction
	LatticePropagator local_insert_prop = 
	  anti_d_prop*Gamma(gamma_value)*forw_d_prop;

	local_contract[0] =
	  rhopiContract(u_x2, u_x2, local_insert_prop, sigma);

	local_contract[1] =
	  pDeltaContract(u_x2, u_x2, local_insert_prop, sigma);

	// Construct the non-local (possibly conserved) current contraction
	LatticePropagator nonlocal_insert_prop = 
	  nonlocalCurrentProp(u, mu, forw_d_prop, anti_d_prop);

	nonlocal_contract[0] =
	  rhopiContract(u_x2, u_x2, nonlocal_insert_prop, sigma);

	nonlocal_contract[1] =
	  pDeltaContract(u_x2, u_x2, nonlocal_insert_prop, sigma);
      }
      break;

      default:
	QDP_error_exit("Unknown ud type", ud);
      }


      // Loop over "delta->p"=0 or "p->delta"=1 types of form-factors
      for(int dp = 0; dp < 2; ++dp)
      {
	// Loop over insertions types - these are the spin projectors
	for (int proj = 0; proj < Nd; ++proj) 
	{
	  int seq_src = ud + 2*(proj + Nd*(dp));   // encode which form-factor
	  if (2*Nd*2 != form.formFacs.size())
	  {
	    QDPIO::cerr << "wallRhoPiff: internal sizing error" << endl;
	    QDP_abort(1);
	  }

	  form.formFacs[seq_src].seq_src = seq_src;

	  // The local non-conserved vector-current matrix element 
	  LatticeComplex corr_local_fn = traceSpin(S_proj[proj] * local_contract[dp]);

	  // The nonlocal (possibly conserved) current matrix element 
	  LatticeComplex corr_nonlocal_fn = traceSpin(S_proj[proj] * nonlocal_contract[dp]);
	
	  multi2d<DComplex> hsum_local = phases.sft(corr_local_fn);
	  multi2d<DComplex> hsum_nonlocal = phases.sft(corr_nonlocal_fn);
 
 
	  form.formFacs[seq_src].insertions[mu].gamma_value = gamma_value;
	  form.formFacs[seq_src].insertions[mu].momenta.resize(phases.numMom());  // hold momenta output

	  // Loop over insertion momenta
	  for(int inser_mom_num=0; inser_mom_num < phases.numMom(); ++inser_mom_num) 
	  {
	    form.formFacs[seq_src].insertions[mu].momenta[inser_mom_num].inser_mom_num =
	      inser_mom_num;
	    form.formFacs[seq_src].insertions[mu].momenta[inser_mom_num].inser_mom = 
	      phases.numToMom(inser_mom_num);

	    multi1d<Complex> local_cur3ptfn(length); // always compute
	    multi1d<Complex> nonlocal_cur3ptfn(length); // always compute
	    
	    for (int t=0; t < length; ++t) 
	    {
	      int t_eff = (t - t0 + length) % length;

	      local_cur3ptfn[t_eff] = Complex(hsum_local[inser_mom_num][t]);
	      nonlocal_cur3ptfn[t_eff] = Complex(hsum_nonlocal[inser_mom_num][t]);
	    } // end for(t)

	    form.formFacs[seq_src].insertions[mu].momenta[inser_mom_num].local_current    = local_cur3ptfn;
	    form.formFacs[seq_src].insertions[mu].momenta[inser_mom_num].nonlocal_current = nonlocal_cur3ptfn;
	  } // end for(inser_mom_num)
	} // end for(proj)
      }  // end for(dp)
    } // end for(ud)
  } // end for(mu)

  // Finally, dump the structure
  write(xml, "FormFac", form);

  END_CODE("wallRhoPiFormFac");
}

