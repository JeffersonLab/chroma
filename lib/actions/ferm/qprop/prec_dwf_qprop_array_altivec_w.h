// -*- C++ -*-
// $Id: prec_dwf_qprop_array_altivec_w.h,v 3.3 2006-06-11 20:28:13 edwards Exp $
/*! \file
 *  \brief 4D style even-odd preconditioned domain-wall fermion action
 */

#ifndef __prec_dwf_qprop_array_altivec_w_h__
#define __prec_dwf_qprop_array_altivec_w_h__

#include "fermact.h"
#include "io/aniso_io.h"


extern "C" 
{
  struct ALTIVEC_DWF_Gauge;  // forward decl
}

namespace Chroma
{

  //! ALTIVEC Propagator DWF qpropT
  /*! \ingroup qprop
   *
   * Propagator solver for DWF fermions
   */
  class ALTIVECDWFQpropT : public SystemSolverArray<LatticeFermion>
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

#if 0
    //! Constructor
    /*!
     * Preferred constructor
     *
     * \param m_q_       quark mass ( Read )
     */
    ALTIVECDWFQpropT(Handle< FermState<T,P,Q> > state_, 
		     const Real& OverMass_,
		     const Real& Mass_,
		     int N5_,
		     const AnisoParam_t& anisoParam_,
		     const InvertParam_t& invParam_) : 
      OverMass(OverMass_), Mass(Mass_), 
      N5(N5_), anisoParam(anisoParam_), invParam(invParam_) 
      {init(state_);}
#endif

    //! Alternative constructor for compatibility
    /*!
     * \param m_q_       quark mass ( Read )
     */
    ALTIVECDWFQpropT(Handle< EvenOddPrecConstDetLinearOperatorArray<T,P,Q> > A_,
		     Handle< FermState<T,P,Q> > state_, 
		     const Real& OverMass_,
		     const Real& Mass_,
		     const AnisoParam_t& anisoParam_,
		     const InvertParam_t& invParam_) : 
      A(A_), OverMass(OverMass_), Mass(Mass_), 
      N5(A->size()), anisoParam(anisoParam_), invParam(invParam_) 
      {init(state_);}

    //! Need a real destructor
    ~ALTIVECDWFQpropT() {fini();}

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
    SystemSolverResults_t operator() (multi1d<LatticeFermion>& psi, const multi1d<LatticeFermion>& chi) const;

  protected:
    //! Private internal initializer
    void init(Handle< FermState<T,P,Q> > state);

    //! Private internal destructor
    void fini();
      
  private:
    ALTIVEC_DWF_Gauge *g;

    Handle< EvenOddPrecConstDetLinearOperatorArray<T,P,Q> > A;
    const Real OverMass;
    const Real Mass;
    const int  N5;
    const AnisoParam_t anisoParam;
    const InvertParam_t invParam;
  };

}

#endif
