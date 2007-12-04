// -*- C++ -*-
// $Id: unprec_s_cprec_t_wilson_linop_w.h,v 1.4 2007-12-04 16:04:42 bjoo Exp $
/*! \file
 *  \brief Unpreconditioned Wilson fermion linear operator
 */

#ifndef UNPREC_S_CPREC_T_WILSON_LINOP_H
#define UNPREC_S_CPREC_T_WILSON_LINOP_H

#include "qdp_config.h"

#if QDP_NS == 4
#if QDP_ND == 4
#if QDP_NC == 3

#include "linearop.h"
#include "central_tprec_linop.h"
#include "actions/ferm/linop/dslash_w.h"
#include "actions/ferm/linop/central_tprec_nospin_utils.h"

namespace Chroma 
{ 
  //! Wilson Dirac Operator - Unpreconditioned in Space, Centrally Preconditioned in time
  /*!
   * \ingroup linop
   *
   * This routine is specific to Wilson fermions!
   *
   *                                                      ~      ~+
   * This subroutine applies the unpreconditioned matrix  M  or  M   the vector
   * Psi,
   *
   *      	       	   {   ~
   *      	       	   {   M(U) . Psi      	       if  ISign = PLUS
   *      	   Chi  =  {
   *      	       	   {   ~   +
   *      	       	   {   M(U)  . Psi     	       if  ISign = MINUS
   *
   */
  
  class UnprecSCprecTWilsonLinOp : 
    public UnprecSpaceCentralPrecTimeLinearOperator<LatticeFermion, 
    multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;
    typedef PScalar<PColorMatrix<RComplex<REAL>, 3> > CMat;              // Useful type: ColorMat with no Outer<>
    typedef PSpinVector<PColorVector<RComplex<REAL>, 3>, 2> HVec_site;   // Useful type: Half Vec with no Outer<>
    //! Partial constructor
    UnprecSCprecTWilsonLinOp() {}


    //! Full constructor with Anisotropy
    UnprecSCprecTWilsonLinOp(Handle< FermState<T,P,Q> > fs_,
		      const Real& Mass_,
		      const AnisoParam_t& aniso_)
      {create(fs_,Mass_,aniso_);}

    //! Destructor is automatic
    ~UnprecSCprecTWilsonLinOp() {}

    //! Return the fermion BC object for this linear operator
    const FermBC<T,P,Q>& getFermBC() const {return fs->getBC();}


    //! Creation routine with Anisotropy
    void create(Handle< FermState<T,P,Q> > fs_,
		const Real& Mass_,
		const AnisoParam_t& aniso_);


    //! The time direction
    
    int tDir() const { 
      // Always 3 for now?
      return 3; 
    }

    //! Apply inv (C_L)^{-1}
    void invCLeftLinOp(T& chi, const T& psi, enum PlusMinus isign) const;

    //! Apply inv (C_R)^{-1}
    void invCRightLinOp(T& chi, const T& psi, enum PlusMinus isign) const;

    //! Apply C_L
    void cLeftLinOp(T& chi, const T& psi, enum PlusMinus isign) const;
    
    //! Apply C_R
    void cRightLinOp(T& chi, const T& psi, enum PlusMinus isign) const;


    //! Apply the the space block onto a source vector
    void spaceLinOp(T& chi, const T& psi, enum PlusMinus isign) const
    {
      Real mhalf=Real(-0.5);
      Dw3D.apply(chi, psi, isign,0);
      Dw3D.apply(chi, psi, isign,1);
      chi *= mhalf;
    }

        
    //! Apply dD_s/dU Y \outer X = Tr { X dD_s/dU Y } with X, Y fermion fields 
    void derivSpaceOp(P& ds_u, const T& X, const T& Y, 
		      enum PlusMinus isign) const { 
      Dw3D.deriv(ds_u, X, Y, isign);
    }

    //! Apply d C_R /dU Y \outer X = Tr { X d C_R /dU Y } with X, Y fermion fields 
    void derivCROp(P& ds_u, const T& X, const T& Y, enum PlusMinus isign) const {
      QDPIO::cout << "Not Yet Implemented " << endl;
      QDP_abort(1);
    }
    //! Apply d C_L /dU Y \outer X = Tr { X d C_L /dU Y } with X, Y fermion fields 
    void derivCLOp(P& ds_u, const T& X, const T& Y, enum PlusMinus isign) const{
      QDPIO::cout << "Not Yet Implemented " << endl;
      QDP_abort(1);
    }
    
    //! Apply d/dt ( log C_R ) = C_R^{-1} d/dt C_R
    void derivLogCR(P& ds_u) const {
      QDPIO::cout << "Not Yet Implemented " << endl;
      QDP_abort(1);
    }

    //! Apply d/dt ( log C_L ) = C_L^{-1} d/dt C_L
    void derivLogCL(P& ds_u) const
    { 
      QDPIO::cout << "Not Yet Implemented " << endl;
      QDP_abort(1);
    }
    
    //! Flopcounter
    unsigned long nFlops() const 
    { 
      //      QDPIO::cout << "Flopcount Not Yet Implemented " << endl;
      return 0;
    }



  private:

    AnisoParam_t aniso;
    Real fact;  // tmp holding  Nd+Mass
    Real invfact;
    Real Mass;
    Handle< FermState<T,P,Q> > fs;
    multi1d<LatticeColorMatrix> u;

    multi2d< int  > tsite; // One row for each spatial site, one column for each timeslice
    multi2d< CMat > P_mat; // Same thing but holding the matrices
    multi2d< CMat > P_mat_dag; // Same thing but holding the matrices

    multi1d< CMat > Q_mat_inv;        // Just one matrix for each spatial site ( 1+P[t=0] )^{-1}
    multi1d< CMat > Q_mat_dag_inv;    // Just one matrix for each spatial site ( 1+P_dag[t=Nt-1])^{-1}

    WilsonDslash3D Dw3D;

  };

} // End Namespace Chroma


#endif
#endif
#endif

#endif
