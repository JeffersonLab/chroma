// $Id: dwf_fermact_qprop_array_w.cc,v 1.2 2004-12-09 03:58:04 edwards Exp $
/*! \file
 *  \brief Base class for unprec and even-odd preconditioned DWF qprop
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/prec_dwf_fermact_base_array_w.h"
#include "actions/ferm/fermacts/unprec_dwf_fermact_base_array_w.h"
#include "actions/ferm/linop/dwffld_w.h"


using namespace QDP;

namespace Chroma { 

//! Propagator DWF linear operator
/*!
 * \param psi      quark propagator ( Modify )
 * \param state    gauge field ( Read )
 * \param chi      source ( Read )
 * \param invParam inverter params ( Read (
 * \param ncg_had  number of CG iterations ( Write )
 */

template<typename T, template<class> class C>
void qprop_t(const C<T>& me,
	     T& psi, 
	     Handle<const ConnectState> state, 
	     const T& chi, 
	     const InvertParam_t& invParam,
	     int& ncg_had)
{
  START_CODE();

  const int  N5 = me.size();   // array size better match
  const Real m_q = me.quark_mass();
  int n_count;
  
  // Initialize the 5D fields
  multi1d<T> chi5(N5);
  {
    //  chi5 = (chi,0,0,0,..,0)^T
    chi5 = zero;
    chi5[0] = chi;

    // tmp5 = P . chi5
    multi1d<T> tmp5(N5);
    DwfFld(tmp5, chi5, PLUS);

    // chi5 = D5(1) . tmp5 =  D5(1) . P . (chi,0,0,..,0)^T 
    // Create a Pauli-Villars linop and use it for just this part
    Handle<const LinearOperator< multi1d<T> > > B(me.unprecLinOp(state,Real(1)));

    (*B)(chi5, tmp5, PLUS);
  }

  // psi5 = (psi,0,0,0,...,0)^T
  multi1d<T> psi5(N5);
  psi5 = zero;
  psi5[0] = psi;

  // Solve  D5(m_q) . psi5 = chi5
  me.qpropT(psi5, state, chi5, invParam, ncg_had);

  // Overall normalization
  Real ftmp1 = Real(1) / Real(1 - m_q);

  // Project out first slice after  chi5 <- P^(-1) . psi5
  DwfFld(chi5, psi5, MINUS);

  // Normalize and remove contact term
  psi = ftmp1*(chi5[0] - chi);

  END_CODE();
}


template<>
void 
EvenOddPrecDWFermActBaseArray<LatticeFermion>::qprop(LatticeFermion& psi, 
						     Handle<const ConnectState> state, 
						     const LatticeFermion& chi, 
						     const InvertParam_t& invParam,
						     int& ncg_had) const
{
  qprop_t(*this, psi, state, chi, invParam,ncg_had);
}


template<>
void 
UnprecDWFermActBaseArray<LatticeFermion>::qprop(LatticeFermion& psi, 
						Handle<const ConnectState> state, 
						const LatticeFermion& chi, 
						const InvertParam_t& invParam,
						int& ncg_had) const
{
  qprop_t(*this, psi, state, chi, invParam, ncg_had);
}

}; // namespace Chroma
