// $Id: prec_fermact_qprop_array.cc,v 1.10 2004-12-12 21:22:17 edwards Exp $
// $Log: prec_fermact_qprop_array.cc,v $
// Revision 1.10  2004-12-12 21:22:17  edwards
// Major overhaul of fermact and linearop class structure. Merged ApproxLinOp
// into LinOp. Removed dsdu stuff - now in linops.  Introduced levels of
// fermacts (base or not) that have derivs. All low level fermacts/linops
// have an additional template param for type of conjugate momenta in
// derivative. Made all fermacts and linops within chroma namespace.
//
// Revision 1.9  2004/11/17 16:52:24  edwards
// Improved some comments
//
// Revision 1.8  2004/10/08 13:20:15  bjoo
// IBM Fixes
//
// Revision 1.7  2004/09/08 02:48:26  edwards
// Big switch-over of fermact IO. New fermact startup mechanism -
// now using Singleton Factory object. Moved  quarkprop4 to be
// a virtual func with top level FermionAction. Disconnected
// zolo4d and 5d fermacts temporarily. Removing usage of old
// fermact and subsequently md,gaugeact  IO mechanisms.
//
// Revision 1.6.2.1  2004/09/02 16:00:07  bjoo
// Trolled beneath the fermact mountains and changed the invocations
// of the inverters to conform to the new inverter interface (invParam -- vs RsdCG, MaxCG and friends) - also in the qprop_files too.
//
// I HAVE REMOVED ZOLOTAREV4D and ZOLOTAREV5D from the build as these
// will need extensive reworking to keep those params inside them
// (simply too many).  I will get them to conform later. Use the
// production branch if you need access to those
//
// I have added the various LOKI based files (Alexandrescu) to the
// Makefile.am so that things install correctly.
//
// I made the propagator code go through the factory invokation to
// do a propagator calculation with prec wilson fermions.
//
// Interestingly the linkage_hack appeared to be optimised away in its
// old form. I turned it into a function that returns the "foo" within,
// so it cannot be compiled away. Interestingly, it still doesn't actually
// need to be CALLED.
//
// Main achievements:
// 	Library compiles
// 	Propagator runs with suitable fermacts (I should check DWF etc)
//
// Revision 1.6  2004/07/28 02:38:02  edwards
// Changed {START,END}_CODE("foo") to {START,END}_CODE().
//
// Revision 1.5  2004/02/05 20:01:47  kostas
// Fixed the chi_tmp bug to all inverters
//
/*! \file
 *  \brief Propagator solver for a generic even-odd preconditioned fermion operator
 *
 *  Solve for the propagator of a even-odd non-preconditioned fermion operator
 */

#include "chromabase.h"
#include "fermact.h"
#include "invtype.h"
#include "actions/ferm/invert/invcg2_array.h"

using namespace QDP;

namespace Chroma {
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
void qprop_t(const EvenOddPrecWilsonTypeFermActBase< multi1d<T> >& me,
	     multi1d<T>& psi, 
	     Handle<const ConnectState> state, 
	     const multi1d<T>& chi, 
	     const InvertParam_t& invParam,
	     int& ncg_had)
{
  START_CODE();

  int n_count;
  
  /* Construct the linear operator */
  /* This allocates field for the appropriate action */
  Handle<const EvenOddPrecLinearOperatorBase< multi1d<T> > > A(me.linOp(state));

  /* Step (i) */
  /* chi_tmp_o =  chi_o - D_oe * A_ee^-1 * chi_e */
  multi1d<T> chi_tmp(me.size());
  {
    multi1d<T> tmp1(me.size());
    multi1d<T> tmp2(me.size());

    A->evenEvenInvLinOp(tmp1, chi, PLUS);
    A->oddEvenLinOp(tmp2, tmp1, PLUS);
    for(int n=0; n < me.size(); ++n)
      chi_tmp[n][rb[1]] = chi[n] - tmp2[n];
  }

  switch(invParam.invType)
  {
  case CG_INVERTER: 
  {
    /* tmp = M_dag(u) * chi_tmp */
    multi1d<T> tmp(me.size());
    (*A)(tmp, chi_tmp, MINUS);

    /* psi = (M^dag * M)^(-1) chi_tmp */
    InvCG2 (*A, tmp, psi, invParam.RsdCG, invParam.MaxCG, n_count);
  }
  break;
  
#if 0
  case MR_INVERTER:
    /* psi = M^(-1) chi */
    InvMR (*A, chi_tmp, psi, invParam.MRover, invParam.RsdCG, invParam.MaxCG, n_count);
    break;

  case BICG_INVERTER:
    /* psi = M^(-1) chi */
    InvBiCG (*A, chi_tmp, psi, invParam.RsdCG, invParam.MaxCG, n_count);
    break;
#endif
  
  default:
    QDP_error_exit("Unknown inverter type", invParam.invType);
  }
  
  if ( n_count == invParam.MaxCG )
    QDP_error_exit("no convergence in the inverter", n_count);
  
  ncg_had = n_count;
  
  /* Step (ii) */
  /* psi_e = A_ee^-1 * [chi_e  -  D_eo * psi_o] */
  {
    multi1d<T> tmp1(me.size());
    multi1d<T> tmp2(me.size());

    A->evenOddLinOp(tmp1, psi, PLUS);
    for(int n=0; n < me.size(); ++n)
      tmp2[n][rb[0]] = chi[n] - tmp1[n];
    A->evenEvenInvLinOp(psi, tmp2, PLUS);
  }
  
  END_CODE();
}


template<>
void 
EvenOddPrecWilsonTypeFermActBase< multi1d<LatticeFermion> >::qpropT(multi1d<LatticeFermion>& psi, 
								    Handle<const ConnectState> state, 
								    const multi1d<LatticeFermion>& chi, 
								    const InvertParam_t& invParam,
								    int& ncg_had) const
{
  qprop_t<LatticeFermion>(*this, psi, state, chi, invParam, ncg_had);
}

}; // namespace Chroma
