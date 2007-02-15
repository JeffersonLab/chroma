// -*- C++ -*-
// $Id: unprec_s_cprec_t_wilson_linop_w.h,v 1.2 2007-02-15 19:59:42 bjoo Exp $
/*! \file
 *  \brief Unpreconditioned Wilson fermion linear operator
 */

#ifndef UNPREC_S_CPREC_T_WILSON_LINOP_H
#define UNPREC_S_CPREC_T_WILSON_LINOP_H

#include "linearop.h"
#include "central_tprec_linop.h"

#include "actions/ferm/linop/dslash_w.h"


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
    void spaceLinOp(T& chi, const T& psi, enum PlusMinus isign) const;

        
    //! Apply dD_s/dU Y \outer X = Tr { X dD_s/dU Y } with X, Y fermion fields 
    void derivSpaceOp(P& ds_u, const T& X, const T& Y, 
		      enum PlusMinus isign) const { 
      QDPIO::cout << "Not Yet Implemented " << endl;
      QDP_abort(1);
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
    //! Apply the Nt*Nc matrix T
    inline
    void TOp(LatticeHalfFermion& chi, const LatticeHalfFermion& psi, enum PlusMinus isign) const;

    //! Invert T using Sherman-Morrison Woodbury
    inline
    void invTOp(LatticeHalfFermion& chi, const LatticeHalfFermion& psi, enum PlusMinus isign) const;


    //! Invert 3by3 complex matrix (not SU3)
    inline
    void invert3by3( CMat& M_inv, const CMat& M) const;  

    AnisoParam_t aniso;
    Real fact;  // tmp holding  Nd+Mass
    Real invfact;
    Real Mass;
    Handle< FermState<T,P,Q> > fs;
    multi1d<LatticeColorMatrix> u;


    int x_index;
    int y_index;
    int z_index;
    int t_index;
    int Nx;
    int Ny;
    int Nz;
    int Nt;
    
    

    

    multi3d< multi1d<int> > tsite;
    multi3d< multi1d< CMat > > P_mat;
    multi3d< multi1d< CMat > > P_mat_dag;

    multi3d< CMat > Q_mat_inv;
    multi3d< CMat > Q_mat_dag_inv;

  };

} // End Namespace Chroma


#endif
