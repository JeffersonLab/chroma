// -*- C++ -*-
// $Id: prec_dwf_qprop_array_sse_w.h,v 1.5 2005-02-21 19:28:59 edwards Exp $
/*! \file
 *  \brief 4D style even-odd preconditioned domain-wall fermion action
 */

#ifndef __prec_dwf_qprop_array_sse_w_h__
#define __prec_dwf_qprop_array_sse_w_h__

#include "fermact.h"
#include "io/param_io.h"       // to get AnisoParam_t

namespace Chroma
{

  //! SSE Propagator DWF qpropT
  /*! \ingroup qprop
   *
   * Propagator solver for DWF fermions
   */
  class SSEDWFQpropT : public SystemSolver< multi1d<LatticeFermion> >
  {
  public:
    //! Constructor
    /*!
     * Preferred constructor
     *
     * \param m_q_       quark mass ( Read )
     */
    SSEDWFQpropT(Handle<const ConnectState> state_, 
		 const Real& OverMass_,
		 const Real& Mass_,
		 int N5_,
		 const AnisoParam_t& anisoParam_,
		 const InvertParam_t& invParam_) : 
      state(state_), OverMass(OverMass_), Mass(Mass_), 
      N5(N5_), anisoParam(anisoParam_), invParam(invParam_) 
      {init();}

    //! Alternative constructor for compatibility
    /*!
     * \param m_q_       quark mass ( Read )
     */
    SSEDWFQpropT(Handle< const EvenOddPrecLinearOperator< multi1d<LatticeFermion>, multi1d<LatticeColorMatrix> > > A,
		 Handle<const ConnectState> state_, 
		 const Real& OverMass_,
		 const Real& Mass_,
		 const AnisoParam_t& anisoParam_,
		 const InvertParam_t& invParam_) : 
      state(state_), OverMass(OverMass_), Mass(Mass_), 
      N5(A->size()), anisoParam(anisoParam_), invParam(invParam_) 
      {init();}

    //! Need a real destructor
    ~SSEDWFQpropT() {fini();}

    //! Expected length of array index
    int size() const {return N5;}

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
    void init() const;
      
    //! Private internal destructor
    void fini() const;
      
  private:
    Handle<const ConnectState> state;
    const Real OverMass;
    const Real Mass;
    const int  N5;
    const AnisoParam_t anisoParam;
    const InvertParam_t invParam;
  };

}

#endif
