// -*- C++ -*-
// $Id: central_tprec_linop.h,v 3.8 2007-12-18 21:06:47 bjoo Exp $
/*! @file
 * @brief Time-preconditioned Linear Operators
 */

#ifndef central_tprec_linop_h
#define central_tprec_linop_h

#include "qdp_config.h"

#if QDP_NS == 4
#if QDP_NC == 3
#if QDP_ND == 4

#include "linearop.h"

namespace Chroma
{

  //! Central Preconditioned Linear Operators

  //! For now, no forces just yet - come later...
  template<typename T, typename P, typename Q>
  class CentralTimePrecLinearOperator : public DiffLinearOperator<T,P,Q>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~CentralTimePrecLinearOperator() {}

    //! Defined on the entire lattice
    const Subset& subset() const = 0;

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


    //! Apply the operator onto a source vector
    //  This varies a lot amongst the families
    virtual void operator() (T& chi, const T& psi, enum PlusMinus isign) const = 0;


    //! Apply the UNPRECONDITIONED operator onto a source vector
    /*! Mainly intended for debugging */
    /*! This by the way defines the essence of the preconditioning ie:
     *    M_unprec = C^{-1}_L  M_prec C^{-1}_R
     *
     *    This applies whether we use
     *           no spatial preconditioning: M_prec = 1 + C_L D_s C_R
     *           ILU preconditioning         M_prec = L + U + L (A-1) U 
     *                Here L & U involve both D_s and C_L & C_R, with 
     *                A being a left over diagonal piece.
     *                However det(L) = det(C_R), det(U) = det(C_L) 
     *           EO3DPrec                    M_prec = L D U 
     *                here LDU is a Schur Decomposition of M_prec from no 
     *                spatial decomposition case.
     */
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
    virtual void deriv(P& ds_u, const T& X, const T& Y, enum PlusMinus isign) const  = 0;    

    //! Get log det ( T^\dag T )
    virtual Double logDetTDagT(void) const = 0;

    //! Get the force due to the det T^\dag T bit
    virtual void derivLogDetTDagT(P& ds_u, enum PlusMinus isign) const = 0;
  };


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
  class UnprecSpaceCentralPrecTimeLinearOperator : public CentralTimePrecLinearOperator<T,P,Q>
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

    //! Override operator() for particular preconditioning
    virtual void operator() (T& chi, const T& psi, enum PlusMinus isign) const
    {
      T   tmp1, tmp2; moveToFastMemoryHint(tmp1); moveToFastMemoryHint(tmp2);

      switch (isign)
      {
      case PLUS:
	// chi   = ( 1 + C_L D_s C_R ) psi
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


    
    //! Apply the d/dt of the preconditioned linop
    virtual void deriv(P& ds_u, const T& X, const T& Y, enum PlusMinus isign) const  = 0;    

    //! Get log det ( T^\dag T )
    virtual Double logDetTDagT(void) const = 0;

    //! Get the force due to the det T^\dag T bit
    virtual void derivLogDetTDagT(P& ds_u, enum PlusMinus isign) const = 0;


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
  class ILUPrecSpaceCentralPrecTimeLinearOperator : public CentralTimePrecLinearOperator<T,P,Q>
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

    //! Override operator()  For this preconditioning 
    virtual void operator() (T& chi, const T& psi, enum PlusMinus isign) const
    {
      T   tmp1, tmp2; moveToFastMemoryHint(tmp1); moveToFastMemoryHint(tmp2);

      switch (isign)
      {
      case PLUS:
	// Eisenstat's Trick:
	//  chi   = ( C_L + C_R + C_L (A - 1 ) C_R ) psi
	//        =   C_L ( 1 + (A - 1 ) C_R) psi +   C_R psi
	//
	//  eval as  tmp1 = C_R psi
	//           tmp2 = psi + (A-1)tmp1
	//           chi   = C_L tmp2 
	//           chi  += tmp1
	cRightLinOp(tmp1, psi, isign);
	AMinusOneOp(tmp2, tmp1,isign);
	tmp2 +=psi;
	cLeftLinOp(chi, tmp2, isign);
	chi += tmp1;

	break;

      case MINUS:
	// Eisenstat's Trick:

	//  chi   =  ( C_R^\dag + C_L^\dag + C_R^\dag (A-1)^\dag  C_L^\dag ) \psi
	//        =  C^R\dag ( 1 + (A-1)^\dag C_L^\dag ) psi + C_L^\dag \psi

	cLeftLinOp(tmp1, psi, isign);
	AMinusOneOp(tmp2, tmp1, isign);
	tmp2 += psi;
	cRightLinOp(chi, tmp2, isign);
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
    
    //! Get log det ( T^\dag T )
    virtual Double logDetTDagT(void) const = 0;

    //! Get the force due to the det T^\dag T bit
    virtual void derivLogDetTDagT(P& ds_u, enum PlusMinus isign) const = 0;

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
  class EO3DPrecSpaceCentralPrecTimeLinearOperator : public CentralTimePrecLinearOperator<T,P,Q>
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

    //! Apply inv (C_L)^{-1}
    virtual void invCLeftLinOp(T& chi, const T& psi, enum PlusMinus isign) const {
      invCLeftLinOp(chi, psi, isign, 0);
      invCLeftLinOp(chi, psi, isign, 1);
    }

    //! Apply inv (C_R)^{-1}
    virtual void invCRightLinOp(T& chi, const T& psi, enum PlusMinus isign) const {
      invCRightLinOp(chi, psi, isign, 0);
      invCRightLinOp(chi, psi, isign, 1);
    }

    //! Apply C_L
    virtual void cLeftLinOp(T& chi, const T& psi, enum PlusMinus isign) const {
      cLeftLinOp(chi, psi, isign, 0);
      cLeftLinOp(chi, psi, isign, 1);
    }
    
    //! Apply C_R
    virtual void cRightLinOp(T& chi, const T& psi, enum PlusMinus isign) const {
      cRightLinOp(chi, psi, isign, 0);
      cRightLinOp(chi, psi, isign, 1);
    }


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
	  invCRightLinOp(tmp1, psi, isign);

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
	  invCLeftLinOp(chi, tmp2, isign);
	}
	break;
      case MINUS:
	{

	  moveToFastMemoryHint(tmp1); 
	  moveToFastMemoryHint(tmp2);
	  moveToFastMemoryHint(tmp3);


	  // tmp 1 = C_R^{-\dagger} \psi
	  invCLeftLinOp(tmp1, psi, isign);

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

    //! Get log det ( T^\dag T )
    virtual Double logDetTDagT(void) const = 0;

    //! Get the force due to the det T^\dag T bit
    virtual void derivLogDetTDagT(P& ds_u, enum PlusMinus isign) const = 0;

  };


}

#endif
#endif
#endif

#endif
