// -*- C++ -*-
// $Id: central_tprec_linop.h,v 3.1 2007-02-13 22:25:04 bjoo Exp $
/*! @file
 * @brief Time-preconditioned Linear Operators
 */

#ifndef central_tprec_linop_h
#define central_tprec_linop_h

#include "linearop.h"

namespace Chroma
{

  //-----------------------------------------------------------------------------------
  //! Time preconditioned linear operator
  /*! @ingroup linop
   *
   * Support for time preconditioned linear operators
   * Given a matrix M written in block form:
   *
   *  M = D_t  +  D_s
   *
   * The preconditioning consists of defining matrices C_l and C_r
   * 
   * so that C_l^{-1} C_r^{-1} = D_t
   * 
   * and \gamma_5 C_l^{-1} \gamma_5 = C_r^\dagger
   * and \gamma_5 C_r^{-1} \gamma_5 = C_l^\dagger
   *
   * The preconditioned op is then 
   *
   *  M = 1 + C_l D_s C_r
   * 
   * and 
   *
   *  M^\dag = 1 + C_r^\dag D^\dag_s C_l^\dag
   *
   * The simplest way of performing the preconditioning is 
   * To have no other preconditioning than the time preconditioning
   */

  //! For now, no forces just yet - come later...
  template<typename T, typename P, typename Q>
  class UnprecSpaceCentralPrecTimeLinearOperator : public DiffLinearOperator<T,P,Q>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~UnprecSpaceCentralPrecTimeLinearOperator() {}

    //! Defined on the entire lattice
    const OrderedSubset& subset() const {return all;}

    //! Return the fermion BC object for this linear operator
    virtual const FermBC<T,P,Q>& getFermBC() const = 0;

    //! The time direction
    virtual int tDir() const = 0;

    //! Apply inv (C_L)^{-1}
    virtual void invCLeftLinOp(T& chi, const T& psi, enum PlusMinus isign) const = 0;

    //! Apply inv (C_R)^{-1}
    virtual void invCRightLinOp(T& chi, const T& psi, enum PlusMinus isign) const = 0;

    //! Apply C_L
    virtual void cLeftLinOp(T& chi, const T& psi, enum PlusMinus isign) const = 0;
    
    //! Apply C_R
    virtual void cRightLinOp(T& chi, const T& psi, enum PlusMinus isign) const = 0;


    //! Apply the the space block onto a source vector
    virtual void spaceLinOp(T& chi, const T& psi, enum PlusMinus isign) const = 0;

    //! Apply the operator onto a source vector
    virtual void operator() (T& chi, const T& psi, enum PlusMinus isign) const
    {
      T   tmp1, tmp2; moveToFastMemoryHint(tmp1); moveToFastMemoryHint(tmp2);

      switch (isign)
      {
      case PLUS:
	//  chi   = ( 1 + C_L D_s C_R ) psi
	cRightLinOp(tmp1, psi, isign);
	spaceLinOp(tmp2, tmp1, isign);
	cLeftLinOp(tmp1, tmp2, isign);
	chi = psi + tmp1;
	break;

      case MINUS:
	//  chi   =  (1 + C_R^\dag D_s C_L^\dag ) psi
	cLeftLinOp(tmp1, psi, isign);
	spaceLinOp(tmp2, tmp1, isign);
	cRightLinOp(tmp1, tmp2, isign);
	chi = psi + tmp1;

	break;

      default:
	QDPIO::cerr << "unknown sign" << endl;
	QDP_abort(1);
      }
    }

    //! Apply the UNPRECONDITIONED operator onto a source vector
    /*! Mainly intended for debugging */
    virtual void unprecLinOp(T& chi, const T& psi, 
			     enum PlusMinus isign) const
    {

      switch (isign) { 
      case PLUS:
	{
	  T   tmp1, tmp2;  moveToFastMemoryHint(tmp1); moveToFastMemoryHint(tmp2);

	  invCRightLinOp(tmp1, psi, isign);
	  (*this)(tmp2, tmp1, isign);
	  invCLeftLinOp(chi, tmp2, isign);
	}
	break;
      case MINUS:
	{
	  T   tmp1, tmp2;  moveToFastMemoryHint(tmp1); moveToFastMemoryHint(tmp2);

	  invCLeftLinOp(tmp1, psi, isign);
	  (*this)(tmp2, tmp1, isign);
	  invCRightLinOp(chi, tmp2, isign);
	}
	break;
      default:
	QDPIO::cerr << "unknown sign" << endl;
	QDP_abort(1);
      }

    }
    
    //! Apply dD_s/dU Y \outer X = Tr { X dD_s/dU Y } with X, Y fermion fields 
    virtual void derivSpaceOp(P& ds_u, const T& X, const T& Y, enum PlusMinus isign) const = 0;
    //! Apply d C_R /dU Y \outer X = Tr { X d C_R /dU Y } with X, Y fermion fields 
    virtual void derivCROp(P& ds_u, const T& X, const T& Y, enum PlusMinus isign) const =0;

    //! Apply d C_L /dU Y \outer X = Tr { X d C_L /dU Y } with X, Y fermion fields 
    virtual void derivCLOp(P& ds_u, const T& X, const T& Y, enum PlusMinus isign) const =0;

    //! Apply d/dt ( log C_R ) = C_R^{-1} d/dt C_R
    virtual void derivLogCR(P& ds_u) const =0;

    //! Apply d/dt ( log C_L ) = C_L^{-1} d/dt C_L
    virtual void derivLogCL(P& ds_u) const =0;
    
    //! Apply the d/dt of the preconditioned linop
    virtual void deriv(P& ds_u, const T& X, const T& Y, enum PlusMinus isign) const {
      // Chain Rule
      P tmp(Nd);
      ds_u.resize(Nd);
      T temp_ferm;

      switch( isign ) {
      case PLUS: 
	{
	  // M = 1 + C_L D_S C_R
	  //
	  // X ( dM/dU ) Y = X  [ d C_L/ dU ] (D_s C_R Y)
	  //              +  X C_L [ d D_s / d U ] C_R Y
	  //              +  X  C_L D_s [ d C_R / d U ] Y
	  T cR_Y;
	  T cLdag_X;

	  cRightLinOp(cR_Y, Y, PLUS);
	  cLeftLinOp(cLdag_X, X, MINUS);

          //  (D_s C_R Y)
	  spaceLinOp(temp_ferm, cR_Y, PLUS);

	  //  X  [ d C_L/ dU ] (D_s C_R Y) 
	  derivCLOp(ds_u, X, temp_ferm, PLUS);
	  // X C_L [ d D_s / d U ] C_R Y
	  derivSpaceOp(tmp, cLdag_X, cR_Y, PLUS);
	  ds_u += tmp;


	  // tmp_ferm^\dag = X^\dag C_l D_s => temp_ferm = D_s^\dag C_l^\dag X
	  spaceLinOp(temp_ferm, cLdag_X, MINUS);
	  
	  //  X  C_L D_s [ d C_R / d U ] Y
	  derivCROp(tmp, temp_ferm, Y, PLUS);
	  ds_u += tmp;
	}
	break;
      case MINUS:
	{
	  // M = 1 + C_R^\dag D_S^\dag C_L^\dag
	  //
	  // X ( dM/dU ) Y = X^\dag  [ d C_R^\dag / dU ] (D_s^\dag C_L^\dag Y)
	  //              +  X^\dag  C_R^\dag  [ d D^\dag_s / d U ] C_L^\dag Y
	  //              +  X^\dag  C_R^\dag  D^\dag _s [ d C_L^\dag / d U ] Y
	  T cR_X;
	  T cLdag_Y;

	  // C_L^\dag Y
	  cLeftLinOp(cLdag_Y, Y, MINUS);

	  // X^\dag  C_R^\dag =>  (C_R X)^\dagger
	  cRightLinOp(cR_X, X, PLUS);

          //  (D^\dag_s C_L^\dag Y )
	  spaceLinOp(temp_ferm, cLdag_Y, MINUS);

	  //  X^\dag  [ d C_R^\dag / dU ] (D_s^\dag C_L^\dag Y)
	  derivCROp(ds_u, X, temp_ferm, MINUS);

	  //  X^\dag  C_R^\dag  [ d D^\dag_s / d U ] C_L^\dag Y
	  derivSpaceOp(tmp, cR_X, cLdag_Y, MINUS );
	  ds_u += tmp;


	  // tmp_ferm^\dag = X^\dag C^\dag_R D^\dag_s => temp_ferm = D_s  C_R X
	  spaceLinOp(temp_ferm, cR_X, PLUS);
	  
	  // X^\dag  C_R^\dag  D^\dag _s [ d C_L^\dag / d U ] Y
	  derivCLOp(tmp, temp_ferm, Y, MINUS);
	  ds_u += tmp;
	  
	}
	break;
      default:
	QDPIO::cerr << "unknown sign" << endl;
	QDP_abort(1);
      };
    }

    


  };

}

#endif
