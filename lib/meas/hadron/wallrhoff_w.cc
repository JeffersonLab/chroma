// $Id: wallrhoff_w.cc,v 1.1 2004-06-04 04:02:38 edwards Exp $
/*! \file
 *  \brief Wall-sink rho-> gamma+rho form-factors 
 *
 *  Form factors constructed from a quark and a backward quark propagator
 */

#include "chromabase.h"
#include "util/ft/sftmom.h"
#include "meas/hadron/wallff_w.h"
#include "meas/hadron/wallrhoff_w.h"

using namespace QDP;


//! Wall-sink rho-> gamma+rho form-factors 
/*!
 * \ingroup hadron
 *
 * This routine is specific to Wilson fermions!
 *
 * \param form               Mega-structure holding form-factors ( Write )
 * \param u                  gauge fields (used for non-local currents) ( Read )
 * \param forw_u_prop        forward U quark propagator ( Read )
 * \param back_u_prop        backward D quark propagator ( Read )
 * \param forw_d_prop        forward U quark propagator ( Read )
 * \param back_d_prop        backward D quark propagator ( Read )
 * \param phases             fourier transform phase factors ( Read )
 * \param t0                 time coordinates of the source ( Read )
 * \param t_sink             time coordinates of the sink ( Read )
 */

void wallRhoFormFac(WallFormFac_formfacs_t& form,
		    const multi1d<LatticeColorMatrix>& u, 
		    const LatticePropagator& forw_u_prop,
		    const LatticePropagator& back_u_prop, 
		    const LatticePropagator& forw_d_prop,
		    const LatticePropagator& back_d_prop, 
		    const SftMom& phases,
		    int t0, int t_sink)
{
  START_CODE("wallRhoFormFac");

  if ( Ns != 4 || Nc != 3 || Nd != 4 )	// Code is specific to Ns=4, Nc=3, Nd=4
    return;

  // Length of lattice in j_decay direction and 3pt correlations fcns
  int length = phases.numSubsets();

  int G5 = Ns*Ns-1;
  
  // The list of meaningful insertions
  //   The gamma matrices specified correspond to 
  // V_mu and A_mu = gamma_mu gamma_5, specificially
  // where in the Chroma code gamma_5 = g_3 g_2 g_1 g_0
  //
  //               GAMMA              CURRENT
  //                 1                  V_0
  //                 2                  V_1
  //                 4                  V_2
  //                 8                  V_3
  //                 14                -A_0
  //                 13                 A_1
  //                 11                -A_2
  //                 7                  A_3
  multi1d<int> gamma_list(2*Nd);
  gamma_list[0] = 1;
  gamma_list[1] = 2;
  gamma_list[2] = 4;
  gamma_list[3] = 8;
  gamma_list[4] = 14;
  gamma_list[5] = 13;
  gamma_list[6] = 11;
  gamma_list[7] = 7;

  // Project propagator onto zero momentum: Do a slice-wise sum.
  Propagator u_x2 = sum(forw_u_prop, phases.getSet()[t_sink]);
  Propagator d_x2 = sum(forw_d_prop, phases.getSet()[t_sink]);
  LatticePropagator anti_u_prop = adj(Gamma(G5)*back_u_prop*Gamma(G5));
  LatticePropagator anti_d_prop = adj(Gamma(G5)*back_d_prop*Gamma(G5));


  // Resize some things - this is needed upfront because I traverse the 
  // structure in a non-recursive scheme
  form.quark.resize(2);
  for (int ud=0; ud < form.quark.size(); ++ud) 
  {
    form.quark[ud].lorentz.resize(Nd*Nd);
    for(int lorz = 0; lorz < form.quark[ud].lorentz.size(); ++lorz)
    {
      form.quark[ud].lorentz[lorz].formfac.resize(1);
      for(int dp = 0; dp < form.quark[ud].lorentz[lorz].formfac.size(); ++dp)
      {
	form.quark[ud].lorentz[lorz].formfac[dp].projector.resize(1);
	for (int proj = 0; proj < form.quark[ud].lorentz[lorz].formfac[dp].projector.size(); ++proj) 
	{
	  form.quark[ud].lorentz[lorz].formfac[dp].projector[proj].insertion.resize(gamma_list.size());
	}
      }
    }
  }


  // For calculational purpose, loop over insertions first.
  // This is out-of-order from storage within the data structure
  // Loop over gamma matrices of the insertion current of insertion current
  for(int gamma_ctr = 0; gamma_ctr < gamma_list.size(); ++gamma_ctr)
  {
    int gamma_value = gamma_list[gamma_ctr];
    int mu = gamma_value % Nd;
    bool compute_nonlocal = (gamma_ctr < Nd) ? true : false;

    LatticePropagator local_insert_prop, nonlocal_insert_prop;

    // Loop over "u"=0 or "d"=1 pieces
    for(int ud = 0; ud < form.quark.size(); ++ud)
    {
      WallFormFac_quark_t& quark = form.quark[ud];

      switch (ud)
      {
      case 0:
      {
	// "\bar u O u" insertion in rho
	// The local non-conserved current contraction
	local_insert_prop = anti_u_prop*Gamma(gamma_value)*forw_u_prop;

	if (compute_nonlocal)
	{
	  // Construct the non-local (possibly conserved) current contraction
	  nonlocal_insert_prop = nonlocalCurrentProp(u, mu, forw_u_prop, anti_u_prop);
	}
      }
      break;

      case 1:
      {
	// "\bar d O d" insertion in rho
	// The local non-conserved current contraction
	LatticePropagator local_insert_prop = 
	  anti_d_prop*Gamma(gamma_value)*forw_d_prop;

	if (compute_nonlocal)
	{
	  // Construct the non-local (possibly conserved) current contraction
	  LatticePropagator nonlocal_insert_prop = 
	    nonlocalCurrentProp(u, mu, forw_d_prop, anti_d_prop);
	}
      }
      break;

      default:
	QDP_error_exit("Unknown ud type", ud);
      }


      // Loop over Lorentz indices of source and sink hadron operators
      for(int lorz = 0; lorz < quark.lorentz.size(); ++lorz)
      {
	WallFormFac_lorentz_t& lorentz = quark.lorentz[lorz];

	int sigma;  // Lorentz index at sink
	int tau;    // Lorentz index at source

	// Encoding of lorentz indices
	sigma = lorz % Nd;      
	tau   = int(lorz / Nd);

	int gamma_value1 = 1 << sigma;
	int gamma_value2 = 1 << tau;

	lorentz.snk_gamma = gamma_value1;
	lorentz.src_gamma = gamma_value2;

	multi1d<LatticeComplex> local_contract(lorentz.formfac.size());
	multi1d<LatticeComplex> nonlocal_contract(lorentz.formfac.size());

	// Contractions depend on "ud" (u or d quark contribution)
	// and source/sink operators
	switch (ud)
	{
	case 0:
	{
	  // "\bar u O u" insertion in rho
	  // The local non-conserved current contraction
	  local_contract[0] = 
	    trace(Gamma(gamma_value1)*local_insert_prop*Gamma(gamma_value2)*
		  Gamma(G5)*adj(d_x2)*Gamma(G5));

	  if (compute_nonlocal)
	  {
	    // Construct the non-local (possibly conserved) current contraction
	    nonlocal_contract[0] = 
	      trace(Gamma(gamma_value1)*nonlocal_insert_prop*Gamma(gamma_value2)*
		    Gamma(G5)*adj(d_x2)*Gamma(G5));
	  }
	}
	break;

	case 1:
	{
	  // "\bar d O d" insertion in rho
	  // The local non-conserved current contraction
	  local_contract[0] =
	    trace(Gamma(gamma_value1)*u_x2*Gamma(gamma_value2)*
		  Gamma(G5)*adj(local_insert_prop)*Gamma(G5));

	  if (compute_nonlocal)
	  {
	    // Construct the non-local (possibly conserved) current contraction
	    nonlocal_contract[0] =
	      trace(Gamma(gamma_value1)*u_x2*Gamma(gamma_value2)*
		    Gamma(G5)*adj(nonlocal_insert_prop)*Gamma(G5));
	  }
	}
	break;

	default:
	  QDP_error_exit("Unknown ud type", ud);
	}


	// Loop over "rho->rho"=0 types of form-factors
	for(int dp = 0; dp < lorentz.formfac.size(); ++dp)
	{
	  WallFormFac_formfac_t& formfac = lorentz.formfac[dp];

	  // Loop over the spin projectors - none here
	  for (int proj = 0; proj < formfac.projector.size(); ++proj) 
	  {
	    WallFormFac_projector_t& projector = formfac.projector[proj];
	    WallFormFac_insertion_t& insertion = projector.insertion[gamma_ctr];

	    insertion.gamma_ctr   = gamma_ctr;
	    insertion.mu          = mu;
	    insertion.gamma_value = gamma_value;
	    insertion.momenta.resize(phases.numMom());  // hold momenta output

	    // The local non-conserved vector-current matrix element 
	    LatticeComplex corr_local_fn = local_contract[0];

	    // The nonlocal (possibly conserved) current matrix element 
	    LatticeComplex corr_nonlocal_fn = nonlocal_contract[0];
	
	    multi1d<WallFormFac_momenta_t>& momenta = insertion.momenta;

	    wallFormFacSft(momenta, corr_local_fn, corr_nonlocal_fn, phases,
			   compute_nonlocal, t0, t_sink);

	  } // end for(proj)
	}  // end for(dp)
      } // end for(lorz)
    } // end for(ud)
  } // end for(gamma_ctr)

  END_CODE("wallRhoFormFac");
}

