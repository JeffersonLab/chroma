// -*- C++ -*-
// $Id: prec_dwf_qprop_array_sse_w.h,v 1.2 2005-01-06 03:50:48 edwards Exp $
/*! \file
 *  \brief 4D style even-odd preconditioned domain-wall fermion action
 */

#ifndef __prec_dwf_qprop_array_sse_w_h__
#define __prec_dwf_qprop_array_sse_w_h__

#include "fermact.h"

using namespace QDP;

namespace Chroma
{

  //! SSE Propagator DWF qpropT
  /*! \ingroup qprop
   *
   * Propagator solver for DWF fermions
   */
  class SSEDWFQprop : public SystemSolver< multi1d<LatticeFermion> >
  {
  public:
    //! Constructor
    /*!
     * \param qpropT_    5D solver ( Read )
     * \param PV_        Pauli-Villars linear operator ( Read )
     * \param m_q_       quark mass ( Read )
     */
    SSEDWFQprop(Handle<const ConnectState> state_, 
		const Real& OverMass_,
		const Real& Mass_,
		int N5_,
		const InvertParam_t& invParam_) : 
      state(state_), OverMass(OverMass_), Mass(Mass_), N5(N5_), invParam(invParam_) {init();}

    //! Need a real destructor
    ~SSEDWFQprop();

    //! Return the subset on which the operator acts
    const OrderedSubset& subset() const {return all;}

    //! Solver the linear system
    /*!
     * \param psi      quark propagator ( Modify )
     * \param chi      source ( Read )
     * \return number of CG iterations
     */
    int operator() (multi1d<LatticeFermion>& psi, const multi1d<LatticeFermion>& chi) const;

  protected:
    //! Private internal initializer
    void init();
      
  private:
    Handle<const ConnectState> state;
    const Real OverMass;
    const Real Mass;
    const int  N5;
    const InvertParam_t invParam;
  };

}

#endif
