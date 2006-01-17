// -*- C++ -*-
// $Id: prec_dwf_qprop_array_altivec_w.h,v 2.2 2006-01-17 16:01:46 bjoo Exp $
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
  class ALTIVECDWFQpropT : public SystemSolver< multi1d<LatticeFermion> >
  {
  public:
    //! Constructor
    /*!
     * Preferred constructor
     *
     * \param m_q_       quark mass ( Read )
     */
    ALTIVECDWFQpropT(Handle<const ConnectState> state_, 
		 const Real& OverMass_,
		 const Real& Mass_,
		 int N5_,
		 const AnisoParam_t& anisoParam_,
		 const InvertParam_t& invParam_) : 
      OverMass(OverMass_), Mass(Mass_), 
      N5(N5_), anisoParam(anisoParam_), invParam(invParam_) 
      {init(state_);}

    //! Alternative constructor for compatibility
    /*!
     * \param m_q_       quark mass ( Read )
     */
    ALTIVECDWFQpropT(Handle< const EvenOddPrecConstDetLinearOperator< multi1d<LatticeFermion>, multi1d<LatticeColorMatrix> > > A,  // throw away
		 Handle<const ConnectState> state_, 
		 const Real& OverMass_,
		 const Real& Mass_,
		 const AnisoParam_t& anisoParam_,
		 const InvertParam_t& invParam_) : 
      OverMass(OverMass_), Mass(Mass_), 
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
    int operator() (multi1d<LatticeFermion>& psi, const multi1d<LatticeFermion>& chi) const;

  protected:
    //! Private internal initializer
    void init(Handle<const ConnectState> state);

    //! Private internal destructor
    void fini();
      
  private:
    ALTIVEC_DWF_Gauge *g;

    const Real OverMass;
    const Real Mass;
    const int  N5;
    const AnisoParam_t anisoParam;
    const InvertParam_t invParam;
  };

}

#endif
