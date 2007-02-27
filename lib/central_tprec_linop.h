// -*- C++ -*-
// $Id: central_tprec_linop.h,v 3.5 2007-02-27 20:28:34 bjoo Exp $
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
    const Subset& subset() const {return all;}

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



  //-----------------------------------------------------------------------------------
  //! Time preconditioned linear operator
  /*! @ingroup linop
   *
   * Support for time preconditioned linear operators with ILU spatial preconditioning
   * Given a matrix M written in block form:
   *
   *  M = D_t  +  D_s  + A 
   *
   * Preconditioning consists of writing:
   *
   *   D_t = P_{-} T  + P_{+} T^\dagger 
   *
   *  Then define T+ = P_{+} + P_{-} T
   *              T- = P_{-} + P_{+} T^\dagger 
   *
   * T is same as before - essentailly the T+ and T- are like the 
   * C_L^{-1} and C_R^{-1} of the unpreconditioned case.
   *
   *
   * We now write   C_L^{-1}  = [ T+_{ee} Ds_{eo}  ]
   *                            [    0    T+_{oo}  ]
   *
   * and C_R^{-1} = [ T-_{ee}  0       ]
   *                [ Ds_{oe}  T-_{oo} ]
   *
   * and 
   *  M = C_L^{-1} +  C_R^{-1}  + tilda{A}
   *
   * where tilde{A} = A - 1
   *
   * So trivial wilson like case: A = 0, clover case A = -1/4 sigma_munu F_munu
   * 
   * The preconditioning is:    M = C_L^{-1} C_L M C_R C_R^{-1} = C_L^{-1} ( \tilde{M} ) C_R^{-1} 
   * with the preconditioned matrix:
   *
   *  \tilde{M} = C_R + C_L + C_L (A - 1) C_R 
   * 
   *  C_R = [ ( T-_{ee} )^{-1}                                      0           ]
   *        [-( T-_{oo} )^{-1} Ds_{oe} ( T-_{ee} )^{-1}        (T-_{oo})^{-1}   ]
   *
   * and 
   *
   *  C_L = [ ( T+_{ee} )^{-1}    -( T+_{ee} )^{-1} Ds_{eo} ( T+_{oo} )^{-1}   ]
   *        [       0                             T+_{oo}^{-1}                 ]
   *
   *
   *  why is this a preconditioning? Certainly it is an ILU form  C_R is lower, C_L is upper and there is a residual term... C_L C_R
   *
   */

  //! For now, no forces just yet - come later...
  template<typename T, typename P, typename Q>
  class ILUPrecSpaceCentralPrecTimeLinearOperator : public DiffLinearOperator<T,P,Q>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~ILUPrecSpaceCentralPrecTimeLinearOperator() {}

    //! Defined on the entire lattice
    const Subset& subset() const {return all;}

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

    //! Apply A - 1
    virtual void AMinusOneOp(T& chi, const T& psi, enum PlusMinus isign) const = 0;

    //! Apply the operator onto a source vector
    virtual void operator() (T& chi, const T& psi, enum PlusMinus isign) const
    {
      T   tmp1, tmp2; moveToFastMemoryHint(tmp1); moveToFastMemoryHint(tmp2);

      switch (isign)
      {
      case PLUS:
	//  chi   = ( C_R + C_L + C_L (A - 1 ) C_R ) psi
	//         = C_L psi + (1  + C_L (A-1))  C_R psi
	//
	//  chi = C_L psi
	cLeftLinOp(chi, psi, isign);

	// tmp1 = C_R psi
	cRightLinOp(tmp1, psi, isign);

	// chi = C_L psi + C_R psi 
	chi += tmp1;
	
	// tmp1 = C_L (A-1) C_R psi 
	AMinusOneOp(tmp2, tmp1, isign);
	cLeftLinOp(tmp1, tmp2, isign);

	chi += tmp1;
	break;

      case MINUS:
	//  chi   =  ( C_R^\dag + C_L^\dag + C_R^\dag (A-1)^\dag  C_L^\dag ) \psi
	//        =  C_R^\dag psi +(1  +  C_R^\dag (A-1)^\dag) C_L^\dag \psi
	cRightLinOp(chi, psi, isign);
	
	// tmp1 = C_L^\dag \psi
	cLeftLinOp( tmp1, psi, isign);
	
	// chi = C_R^\dag \psi + C_L^\dag \psi
	chi += tmp1;

	// tmp2 = (A-1)^\dag C_L^\dag psi
	AMinusOneOp( tmp2, tmp1, isign);
	cRightLinOp( tmp1, tmp2, isign);

	chi += tmp1;

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

    //! Apply the d/dt of the preconditioned linop
    virtual void deriv(P& ds_u, const T& X, const T& Y, enum PlusMinus isign) const = 0;

  };

  //-----------------------------------------------------------------------------------
  //! Time preconditioned linear operator
  /*! @ingroup linop
   *
   * Support for time preconditioned linear operators with ILU spatial preconditioning
   * Given a matrix M written in block form:
   *
   *  M = D_t  +  D_s 
   *
   * Preconditioning consists of writing:
   *
   *   D_t = P_{-} T  + P_{+} T^\dagger 
   *
   *  Then define C_L^{-1} = P_{+} + P_{-} T
   *              C_R^{-1} = P_{-} + P_{+} T^\dagger 
   *
   * 
   * Now we can write 
   *
   *  C_L M C_R = ( 1 +  C_L D_s C_R ) 
   *
   *
   * And now we write D_s in a 3d even odd preconditioned form:
   *
   * and M  = [ 1                   C_L D_s^{eo} C_R ]
   *          [ C_L D^{oe}_s C_R    1                ]
   *
   *
   * and this can be Schur decomposed as
   * [ 1              0   ]      [     1                                    0   ] [ 1  C_L D^{oe} C_R ]         
   * [C_L D^{oe} C_R  1   ]      [     0  1 - C_L D^{oe}_s C_R C_L D_s^{eo} C_R ] [ 0          1      ]
   *
   */


  //! For now, no forces just yet - come later...
  template<typename T, typename P, typename Q>
  class EO3DPrecSpaceCentralPrecTimeLinearOperator : public DiffLinearOperator<T,P,Q>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~EO3DPrecSpaceCentralPrecTimeLinearOperator() {}

    //! Defined on the cb 1
    const Subset& subset() const {return rb3[1];}

    //! Return the fermion BC object for this linear operator
    virtual const FermBC<T,P,Q>& getFermBC() const = 0;

    //! The time direction
    virtual int tDir() const = 0;

    //! Apply inv (C_L)^{-1}
    virtual void invCLeftLinOp(T& chi, const T& psi, enum PlusMinus isign, int cb3) const = 0;

    //! Apply inv (C_R)^{-1}
    virtual void invCRightLinOp(T& chi, const T& psi, enum PlusMinus isign, int cb3) const = 0;

    //! Apply C_L
    virtual void cLeftLinOp(T& chi, const T& psi, enum PlusMinus isign, int cb3) const = 0;
    
    //! Apply C_R
    virtual void cRightLinOp(T& chi, const T& psi, enum PlusMinus isign, int cb3) const = 0;


    //! 3D preconditioning

    //! M_{ee} where the checkerboarding is in 3d
    virtual void evenEvenLinOp(T& chi, const T& psi, enum PlusMinus isign) const = 0;

    //! M_{eo} where the checkerboarding is in 3d
    virtual void evenOddLinOp(T& chi, const T& psi, enum PlusMinus isign) const = 0;

    //! M_{oe} where the checkerboarding is in 3d
    virtual void oddEvenLinOp(T& chi, const T& psi, enum PlusMinus isign) const = 0;

    //! M_{oo} where the checkerboarding is in 3d
    virtual void oddOddLinOp(T& chi, const T& psi, enum PlusMinus isign) const = 0;

    //! M_{ee}^{-1} where the checkerboarding is in 3d - this is awkward for clover 
    // I think. But let's check it for wilson first
    virtual void evenEvenInvLinOp(T& chi, const T& psi, enum PlusMinus isign) const = 0;



    //! Apply the operator onto a source vector
    virtual void operator() (T& chi, const T& psi, enum PlusMinus isign) const
    {
      T   tmp1, tmp2; moveToFastMemoryHint(tmp1); moveToFastMemoryHint(tmp2);

      switch (isign)
      {
      case PLUS:
	// \chi = ( M_oo - M_oe M^{-1}_ee M_eo ) \psi
	
	// chi = M_oo \psi
	oddOddLinOp(chi, psi, PLUS);

	// tmp2 = M_oe M^{-1}_ee M_eo  \psi
	evenOddLinOp(tmp2, psi, PLUS);
	evenEvenInvLinOp(tmp1, tmp2, PLUS);
	oddEvenLinOp(tmp2, tmp1, PLUS);
	
	chi[rb3[1]] -= tmp2;
	break;

      case MINUS:
	// \chi = ( M_oo - M_oe M^{-1}_ee M_eo )^\dagger \psi
	//      = M^\dagger_oo \psi - M^\dagger_{oe} ( M^{-1}_ee )^\dagger M^\dagger{eo}	
	//
	// NB: Daggering acts on checkerboarding to achieve the result above.

	oddOddLinOp(chi, psi, MINUS);

	// tmp2 = M_oe M^{-1}_ee M_eo  \psi
	evenOddLinOp(tmp2, psi, MINUS);
	evenEvenInvLinOp(tmp1, tmp2, MINUS);
	oddEvenLinOp(tmp2, tmp1, MINUS);

	chi[rb3[1]] -= tmp2;
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

      T   tmp1, tmp2,tmp3;  
      moveToFastMemoryHint(tmp1); 
      moveToFastMemoryHint(tmp2);	  
      moveToFastMemoryHint(tmp3);

      switch (isign) { 
      case PLUS:
	{

	  // tmp 1 = C_R^{-1} \psi
	  invCRightLinOp(tmp1, psi, isign, 0);
	  invCRightLinOp(tmp1, psi, isign, 1);


	  // Now apply ( 1  M_ee^{-1} M_eo ) 
          //           ( 0        1        )
	  evenOddLinOp(tmp2, tmp1, PLUS);
	  evenEvenInvLinOp(tmp3,tmp2, PLUS);
	  tmp1[rb3[0]] += tmp3;

	  // rb[1] part
	  (*this)(tmp2, tmp1, PLUS);

	  // rb[0] part
	  evenEvenLinOp(tmp2, tmp1, PLUS);

	  // Now apply ( 1               0  ) 
          //           (M_oe M_ee^{-1}   1  )
	  evenEvenInvLinOp(tmp1, tmp2, PLUS);
	  oddEvenLinOp(tmp3, tmp1, PLUS);
	  tmp2[rb3[1]] += tmp3;


	  //  chi = C_L^{-1} tmp2
	  invCLeftLinOp(chi, tmp2, isign, 0);
	  invCLeftLinOp(chi, tmp2, isign, 1);
	}
	break;
      case MINUS:
	{

	  moveToFastMemoryHint(tmp1); 
	  moveToFastMemoryHint(tmp2);
	  moveToFastMemoryHint(tmp3);


	  // tmp 1 = C_R^{-\dagger} \psi
	  invCLeftLinOp(tmp1, psi, isign, 0);
	  invCLeftLinOp(tmp1, psi, isign, 1);


	  // Now apply ( 1  M_ee^{-\dagger} M_eo^\dagger ) 
          //           ( 0                      1        )
	  evenOddLinOp(tmp2, tmp1, isign);
	  evenEvenInvLinOp(tmp3,tmp2, isign);
	  tmp1[rb3[0]] += tmp3;

	  // rb[1] part
	  (*this)(tmp2, tmp1, isign);

	  // rb[0] part
	  evenEvenLinOp(tmp2, tmp1, isign);

	  // Now apply ( 1                            0  ) 
          //           (M_oe^\dagger M_ee^{-dagger}   1  )

	  evenEvenInvLinOp(tmp1, tmp2, isign);
	  oddEvenLinOp(tmp3, tmp1, isign);
	  tmp2[rb3[1]] += tmp3;


	  //  chi = C_L^{-\dagger} tmp2
	  invCRightLinOp(chi, tmp2, isign, 0);
	  invCRightLinOp(chi, tmp2, isign, 1);

	}
	break;
      default:
	QDPIO::cerr << "unknown sign" << endl;
	QDP_abort(1);
      }

    }

    //! Apply the d/dt of the preconditioned linop
    virtual void deriv(P& ds_u, const T& X, const T& Y, enum PlusMinus isign) const = 0;

  };


}

#endif
