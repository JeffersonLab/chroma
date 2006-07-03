// $Id: unprec_ppdwf4d_linop_w.cc,v 3.1 2006-07-03 15:26:09 edwards Exp $
/*! \file
 *  \brief Unpreconditioned projected DWF operator to 4D using prec 5D bits
 */

#include "chromabase.h"
#include "actions/ferm/invert/invcg2_array.h"
#include "actions/ferm/linop/unprec_ppdwf4d_linop_w.h"
#include "actions/ferm/linop/dwffld_w.h"


namespace Chroma
{

  // Apply unpreconditioned linear operator
  template<>
  void 
  UnprecPPDWF4DLinOp<LatticeFermion, 
		     multi1d<LatticeColorMatrix>,
		     multi1d<LatticeColorMatrix> >::operator()(LatticeFermion& chi, 
							       const LatticeFermion& psi, 
							       enum PlusMinus isign) const
  {
    START_CODE();

    const int  N5 = size();   // array size better match
    SystemSolverResults_t res;
  
    // Initialize the 5D fields
    multi1d<LatticeFermion> psi5(N5);
    //  psi5 = (psi,0,0,0,..,0)^T
    psi5 = zero;
    psi5[0] = psi;

    // tmp5 = P . psi5
    multi1d<LatticeFermion> tmp5(N5);
    DwfFld(tmp5, psi5, PLUS);

    multi1d<LatticeFermion> chi5(N5);

    switch(isign)
    {
    case PLUS:
    {
      // chi5 = D5(m_q) . tmp5 =  D5(m_q) . P . (chi,0,0,..,0)^T 
      D->unprecLinOp(chi5, tmp5, PLUS);

      // Solve  D5(1) . psi5 = chi5
      psi5 = chi5;
  
      /* Step (i) */
      /* chi_tmp_o =  chi_o - D_oe * A_ee^-1 * chi_e */
      multi1d<LatticeFermion> chi_tmp(N5);
      {
	multi1d<LatticeFermion> tmp1(N5);
	multi1d<LatticeFermion> tmp2(N5);

	PV->evenEvenInvLinOp(tmp1, chi5, PLUS);
	PV->oddEvenLinOp(tmp2, tmp1, PLUS);
	for(int n=0; n < N5; ++n)
	  chi_tmp[n][rb[1]] = chi5[n] - tmp2[n];
      }

//      case CG_INVERTER: 
      {
	/* tmp5 = M_dag(u) * chi_tmp */
	(*PV)(tmp5, chi_tmp, MINUS);

	/* psi5 = (M^dag * M)^(-1) chi_tmp */
	res = InvCG2 (*PV, tmp5, psi5, invParam.RsdCG, invParam.MaxCG);
      }

      /* Step (ii) */
      /* psi_e = A_ee^-1 * [chi_e  -  D_eo * psi_o] */
      {
	multi1d<LatticeFermion> tmp1(N5);
	multi1d<LatticeFermion> tmp2(N5);

	PV->evenOddLinOp(tmp1, psi5, PLUS);
	for(int n=0; n < N5; ++n)
	  tmp2[n][rb[0]] = chi5[n] - tmp1[n];
	PV->evenEvenInvLinOp(psi5, tmp2, PLUS);
      }
    }
    break;

    case MINUS:
    {
      // Solve  D5(1) . tmp5 = chi5
      chi5 = tmp5;

      /* Step (i) */
      /* chi_tmp_o =  chi5_o - D_oe^dag * (A_ee^-1)^dag * chi5_e */
      multi1d<LatticeFermion> chi_tmp(N5);
      {
	multi1d<LatticeFermion> tmp1(N5);
	multi1d<LatticeFermion> tmp2(N5);

	PV->evenEvenInvLinOp(tmp1, chi5, MINUS);
	PV->oddEvenLinOp(tmp2, tmp1, MINUS);
	for(int n=0; n < N5; ++n)
	  chi_tmp[n][rb[1]] = chi5[n] - tmp2[n];
      }

//      case CG_INVERTER: 
      {
	multi1d<LatticeFermion> tmp1(N5);
	for(int n=0; n < N5; ++n)
	  tmp1[n][rb[1]] = zero;
  
	/* tmp1 = (M^dag * M)^(-1) chi_tmp */
	res = InvCG2(*PV, chi_tmp, tmp1, invParam.RsdCG, invParam.MaxCG);

	/* tmp5_o = M(u) * tmp1_o = (M^dag)^(-1) * chi_tmp_o */
	(*PV)(tmp5, tmp1, PLUS);
      }

      /* Step (ii) */
      /* tmp5_e = (A_ee^-1)^dag * [chi5_e  -  D_eo^dag * tmp5_o] */
      {
	multi1d<LatticeFermion> tmp1(N5);
	multi1d<LatticeFermion> tmp2(N5);

	PV->evenOddLinOp(tmp1, tmp5, MINUS);
	for(int n=0; n < N5; ++n)
	  tmp2[n][rb[0]] = chi5[n] - tmp1[n];
	PV->evenEvenInvLinOp(tmp5, tmp2, MINUS);
      }

      // psi5 = D5(m_q)^dag . tmp5 
      D->unprecLinOp(psi5, tmp5, MINUS);
    }
    break;
    }
    
    if ( res.n_count == invParam.MaxCG )
      QDP_error_exit("no convergence in the inverter", res.n_count);
  
    // Project out first slice after  chi <- chi5 <- P^(-1) . psi5
    DwfFld(chi5, psi5, MINUS);
    chi = chi5[0];

    END_CODE();
  }

}
