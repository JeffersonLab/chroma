// $Id: prec_fermact_qprop.cc,v 1.2 2003-12-02 15:44:01 edwards Exp $
/*! \file
 *  \brief Propagator solver for a generic even-odd preconditioned fermion operator
 *
 *  Solve for the propagator of a even-odd non-preconditioned fermion operator
 */

#include "chromabase.h"
#include "fermact.h"
#include "actions/ferm/invert/invcg2.h"

using namespace QDP;

//! Propagator of a generic even-odd preconditioned fermion linear operator
/*! \ingroup qprop
 *
 * This routine is actually generic to all even-odd preconditioned fermions
 *
 * \param psi      initial solution ( Modify )
 * \param state    gauge field ( Read )
 * \param chi      source ( Read )
 * \param invType  inverter type ( Read (
 * \param RsdCG    CG (or MR) residual used here ( Read )
 * \param MaxCG    maximum number of CG iterations ( Read )
 * \param ncg_had  number of CG iterations ( Write )
 */

template<typename T>
static 
void qprop_t(const EvenOddPrecWilsonTypeFermAct<T>& me,
	     T& psi, 
	     const ConnectState& state, 
	     const T& chi, 
	     enum InvType invType,
	     const Real& RsdCG, 
	     int MaxCG, int& ncg_had)
{
  START_CODE("EvenOddPrecWilsonTypeFermAct::qprop");

  int n_count;
  
  /* Construct the linear operator */
  /* This allocates field for the appropriate action */
  const EvenOddPrecLinearOperatorProxy<T> A(me.linOp(state));

  /* Step (i) */
  /* chi_tmp =  chi_o - D_oe * A_ee^-1 * chi_o */
  T chi_tmp;
  {
    T tmp1, tmp2;

    A.evenEvenInvLinOp(tmp1, chi, PLUS);
    A.oddEvenLinOp(tmp2, tmp1, PLUS);
    chi_tmp[rb[1]] = chi - tmp2;
  }

  switch(invType)
  {
  case CG_INVERTER: 
  {
    /* tmp = A_dag(u) * chi_tmp */
    T  tmp;
    A(tmp, chi_tmp, MINUS);
    
    /* psi = (M^dag * M)^(-1) chi */
    InvCG2(A, tmp, psi, RsdCG, MaxCG, n_count);
  }
  break;
  
#if 0
  case MR_INVERTER:
    /* psi = M^(-1) chi_tmp */
    InvMR(A, chi_tmp, psi, MRover, RsdCG, MaxCG, n_count);
    break;

  case BICG_INVERTER:
    /* psi = M^(-1) chi_tmp */
    InvBiCG(A, chi_tmp, psi, RsdCG, MaxCG, n_count);
    break;
#endif
  
  default:
    QDP_error_exit("Unknown inverter type", invType);
  }
  
  if ( n_count == MaxCG )
    QDP_error_exit("no convergence in the inverter", n_count);
  
  ncg_had = n_count;
  
  /* Step (ii) */
  /* psi_e = A_ee^-1 * [chi_e  -  D_eo * psi_o] */
  {
    T tmp1, tmp2;

    A.evenOddLinOp(tmp1, psi, PLUS);
    tmp2[rb[0]] = chi - tmp1;
    A.evenEvenInvLinOp(psi, tmp2, PLUS);
  }
  
  END_CODE("EvenOddPrecWilsonTypeFermAct::qprop");
}


template<>
void 
EvenOddPrecWilsonTypeFermAct<LatticeFermion>::qpropT(LatticeFermion& psi, 
						     const ConnectState& state, 
						     const LatticeFermion& chi, 
						     enum InvType invType,
						     const Real& RsdCG, 
						     int MaxCG, int& ncg_had) const
{
  qprop_t(*this, psi, state, chi, invType, RsdCG, MaxCG, ncg_had);
}

template<>
void 
EvenOddPrecWilsonTypeFermAct<LatticeFermion>::qprop(LatticeFermion& psi, 
						    const ConnectState& state, 
						    const LatticeFermion& chi, 
						    enum InvType invType,
						    const Real& RsdCG, 
						    int MaxCG, int& ncg_had) const
{
  qprop_t(*this, psi, state, chi, invType, RsdCG, MaxCG, ncg_had);
}



template<>
void 
EvenOddPrecWilsonTypeFermAct<LatticeDWFermion>::qpropT(LatticeDWFermion& psi, 
						       const ConnectState& state, 
						       const LatticeDWFermion& chi, 
						       enum InvType invType,
						       const Real& RsdCG, 
						       int MaxCG, int& ncg_had) const
{
  qprop_t(*this, psi, state, chi, invType, RsdCG, MaxCG, ncg_had);
}

template<>
void 
EvenOddPrecWilsonTypeFermAct<LatticeDWFermion>::qprop(LatticeFermion& psi, 
						      const ConnectState& state, 
						      const LatticeFermion& chi, 
						      enum InvType invType,
						      const Real& RsdCG, 
						      int MaxCG, int& ncg_had) const
{
  QDPIO::cerr << "EvenOddPrecWilsonTypeFermAct<DWF>::qprop - this implementation is empty" << endl;
  QDP_abort(1);
}


