// -*- C++ -*-
// $Id: eoprec_constdet_linop.h,v 3.2 2007-02-22 21:11:44 bjoo Exp $
/*! @file
 * @brief Preconditioned 4D Linop with Gauge Independent Even-Even part
 */

#ifndef __eoprec_constdet_linop_h__
#define __eoprec_constdet_linop_h__

#include "eoprec_linop.h"

using namespace QDP::Hints;

namespace Chroma 
{

  //----------------------------------------------------------------
  //! Even-odd preconditioned linear operator
  /*! @ingroup linop
   *
   * Support for even-odd preconditioned linear operators
   * Given a matrix M written in block form:
   *
   *      [      A             D        ]
   *      [       E,E           E,O     ]
   *      [                             ]
   *      [      D             A        ]
   *      [       O,E           O,O     ]
   *
   * The preconditioning consists of using the triangular matrices
   *
   *      [      1              0        ]
   *      [       E,E            E,O     ]
   *  L = [                              ]
   *      [     D     A^(-1)    1        ]
   *      [      O,E   E,E        O,O    ]
   *
   * and
   *
   *      [      1            A^{-1}  D       ]
   *      [       E,E          EE      E,O    ]
   *  U = [                                   ]
   *      [      0                    1       ]
   *      [       O,E                  O,O    ]
   *
   * The preconditioned matrix is formed from
   *
   *  M'   =  L^-1 * M * U^-1
   *
   * where
   *
   *           [      1              0        ]
   *           [       E,E            E,O     ]
   *  L^(-1) = [                              ]
   *           [   - D     A^(-1)    1        ]
   *           [      O,E   E,E        O,O    ]
   *
   * and
   *
   *           [      1            - A^(-1) D       ]
   *           [       E,E            E,E    E,O    ]
   *  U^(-1) = [                                    ]
   *           [      0                1            ]
   *           [       O,E              O,O         ]
   *
   * Resulting in a new  M
   *
   *      [      A                    0                      ]
   *      [       E,E                  E,O                   ]
   *      [                                                  ]
   *      [      0                A     -  D    A^(-1)  D    ]
   *      [       O,E              O,O      O,E   E,E    E,O ]
   *
   *
   * This class is used to implement the resulting linear operator
   *
   *      ~
   *      M  =  A(o,o) - D(o,e) . A^-1(e,e) . D(e,o)
   *
   * where A^{-1}_{ee} is independent of the gauge fields. This
   * means that the det A_{ee} is an irrelevant constant and that
   * the force term due to the A_{ee} part is zero.
   *
   * This structure suits most of the linear operators we use, and 
   * It simplifies the force term.
   * By construction, the linear operator is ONLY defined on the odd subset
   *
   */

  template<typename T, typename P, typename Q>
  class EvenOddPrecConstDetLinearOperator : public EvenOddPrecLinearOperator<T,P,Q>  
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~EvenOddPrecConstDetLinearOperator() {}

    //! Return the fermion BC object for this linear operator
    virtual const FermBC<T,P,Q>& getFermBC() const = 0;

    //! Apply the derivative of the operator onto a source vector
    /*! User should make sure deriv routines do a resize  */
    virtual void deriv(P& ds_u, const T& chi, const T& psi, 
		       enum PlusMinus isign) const
    {
      // Need deriv of  (A_oo - D_oe*Ainv_ee*D_eo*psi_e)
      enum PlusMinus msign = (isign == PLUS) ? MINUS : PLUS;

      //
      // Make sure the deriv routines do a resize !!!
      //
      T   tmp1, tmp2;  // if an array is used here, the space is not reserved
      moveToFastMemoryHint(tmp1);
      moveToFastMemoryHint(tmp2);

      P   ds_1;  // deriv routines should resize

      //
      // NOTE: even with even-odd decomposition, the ds_u will still have contributions
      // on all cb. So, no adding of ds_1 onto ds_u under a subset
      //
            
      //  ds_u  -=  chi^dag * D'_oe * Ainv_ee * D_eo * psi_o
      this->evenOddLinOp(tmp1, psi, isign);
      this->evenEvenInvLinOp(tmp2, tmp1, isign);
      this->derivOddEvenLinOp(ds_u, chi, tmp2, isign);
      

      //  ds_u  -=  chi^dag * D_oe * Ainv_ee * D'_eo * psi_o
      this->evenOddLinOp(tmp1, chi, msign);
      this->evenEvenInvLinOp(tmp2, tmp1, msign);
      this->derivEvenOddLinOp(ds_1, tmp2, psi, isign);
      ds_u += ds_1;

      // This is yucky and ws should have a *= function 
      // so that I wouldn't have to expose the fact that I expect
      // ds_u to be an Nd dimensional array
      for(int mu=0; mu < Nd; mu++) { 
	ds_u[mu] *= Real(-1);
      }

      getFermBC().zero(ds_u);
    }

    //! Apply the even-even block onto a source vector
    virtual void derivEvenEvenLinOp(P& ds_u, const T& chi, const T& psi, 
				    enum PlusMinus isign) const
    {
      QDPIO::cerr << "EvenOdd: not implemented" << endl;
      QDP_abort(1);
    }
  
    //! Apply the the even-odd block onto a source vector
    virtual void derivEvenOddLinOp(P& ds_u, const T& chi, const T& psi, 
				   enum PlusMinus isign) const
    {
      QDPIO::cerr << "EvenOdd: not implemented" << endl;
      QDP_abort(1);
    }
 
    //! Apply the the odd-even block onto a source vector
    virtual void derivOddEvenLinOp(P& ds_u, const T& chi, const T& psi, 
				   enum PlusMinus isign) const
    {
      QDPIO::cerr << "EvenOdd: not implemented" << endl;
      QDP_abort(1);
    }

    //! Apply the the odd-odd block onto a source vector
    virtual void derivOddOddLinOp(P& ds_u, const T& chi, const T& psi, 
				  enum PlusMinus isign) const
    {
      QDPIO::cerr << "EvenOdd: not implemented" << endl;
      QDP_abort(1);
    }

    // Multipole derivatives
    virtual void derivEvenEvenLinOpMP(P& ds_u, const multi1d<T>& chi, const multi1d<T>& psi, 
				       enum PlusMinus isign) const 
   {
     // implement in terms existing derivatives:
     // Trivially zero  since determinant is constant
     ds_u.resize(Nd);
     ds_u = zero;
     

     /*
       F_tmp = zero;
       for(int i=0; i < chi.size(); i++) { 
         derivEvenEvenLinOp(F_tmp, chi[i], psi[i], isign);
	 ds_u += F_tmp;
       }
     */
   }

    // Multipole derivatives
    virtual void derivEvenOddLinOpMP(P& ds_u, const multi1d<T>& chi, const multi1d<T>& psi, 
				       enum PlusMinus isign) const 
   {
     // implement in terms existing derivatives:
     ds_u.resize(Nd);
     ds_u = zero;
     
     P F_tmp; // deriv will resize
     for(int i=0; i < chi.size(); i++) { 
       derivEvenOddLinOp(F_tmp, chi[i], psi[i], isign);
       ds_u += F_tmp;
     }
   }

    virtual void derivOddEvenLinOpMP(P& ds_u, const multi1d<T>& chi, const multi1d<T>& psi, 
				       enum PlusMinus isign) const 
   {
     // implement in terms existing derivatives:
     ds_u.resize(Nd);
     ds_u = zero;
     
     P F_tmp; // deriv will resize
     for(int i=0; i < chi.size(); i++) { 
       derivOddEvenLinOp(F_tmp, chi[i], psi[i], isign);
       ds_u += F_tmp;
     }
   }

    virtual void derivOddOddLinOpMP(P& ds_u, const multi1d<T>& chi, const multi1d<T>& psi, 
				       enum PlusMinus isign) const 
   {
     // Trivially zero since we are constdet?
     ds_u.resize(Nd);
     ds_u = zero;
     
     /*
       F_tmp = zero;
       for(int i=0; i < chi.size(); i++) { 
          derivOddOddLinOp(F_tmp, chi[i], psi[i], isign);
          ds_u += F_tmp;
       }
     */
   }
				      

   //! Apply the derivative of the operator onto a source vector
    /*! User should make sure deriv routines do a resize  */
    virtual void derivMultipole(P& ds_u, const multi1d<T>& chi, const multi1d<T>& psi, 
				enum PlusMinus isign) const
    {
      // Need deriv of  (A_oo - D_oe*Ainv_ee*D_eo*psi_e)
      enum PlusMinus msign = (isign == PLUS) ? MINUS : PLUS;

      //
      // Make sure the deriv routines do a resize !!!
      //
      multi1d<T> tmp2(chi.size());  // Need this for things like  M chi_i

      P   ds_1;  // deriv routines should resize

      //
      // NOTE: even with even-odd decomposition, the ds_u will still have contributions
      // on all cb. So, no adding of ds_1 onto ds_u under a subset
      //
            
      //  ds_u  -=  chi^dag * D'_oe * Ainv_ee * D_eo * psi_o
      for(int i=0; i < chi.size(); i++) { 
	T   tmp1;
	this->evenOddLinOp(tmp1, psi[i], isign);
	this->evenEvenInvLinOp(tmp2[i], tmp1, isign);
      }

      this->derivOddEvenLinOpMP(ds_u, chi, tmp2, isign);


      //  ds_u  -=  chi^dag * D_oe * Ainv_ee * D'_eo * psi_o
      for(int i=0; i < chi.size(); i++) { 
	T tmp1;
	this->evenOddLinOp(tmp1, chi[i], msign);
	this->evenEvenInvLinOp(tmp2[i], tmp1, msign);
      }

      this->derivEvenOddLinOpMP(ds_1, tmp2, psi, isign);
      ds_u += ds_1;

      // This is yucky and ws should have a *= function 
      // so that I wouldn't have to expose the fact that I expect
      // ds_u to be an Nd dimensional array
      for(int mu=0; mu < Nd; mu++) { 
	ds_u[mu] *= Real(-1);
      }

      getFermBC().zero(ds_u);
    }

  };


  //----------------------------------------------------------------
  //! Even-odd preconditioned linear operator including derivatives for arrays
  /*! @ingroup linop
   *
   * Support for even-odd preconditioned linear operators with derivatives
   *
   * Given a matrix M written in block form:
   *
   *      [      A             D        ]
   *      [       E,E           E,O     ]
   *      [                             ]
   *      [      D             A        ]
   *      [       O,E           O,O     ]
   *
   * The preconditioning consists of using the triangular matrices
   *
   *      [      1              0        ]
   *      [       E,E            E,O     ]
   *  L = [                              ]
   *      [     D     A^(-1)    1        ]
   *      [      O,E   E,E        O,O    ]
   *
   * and
   *
   *      [      A              D       ]
   *      [       E,E            E,O    ]
   *  U = [                             ]
   *      [      0              1       ]
   *      [       O,E            O,O    ]
   *
   * The preconditioned matrix is formed from
   *
   *  M'   =  L^-1 * M * U^-1
   *
   * where
   *
   *           [      1              0        ]
   *           [       E,E            E,O     ]
   *  L^(-1) = [                              ]
   *           [   - D     A^(-1)    1        ]
   *           [      O,E   E,E        O,O    ]
   *
   * and
   *
   *           [      A^(-1)       - A^(-1) D       ]
   *           [       E,E            E,E    E,O    ]
   *  U^(-1) = [                                    ]
   *           [      0                1            ]
   *           [       O,E              O,O         ]
   *
   * Resulting in a new  M
   *
   *      [      1                    0                      ]
   *      [       E,E                  E,O                   ]
   *      [                                                  ]
   *      [      0                A     -  D    A^(-1)  D    ]
   *      [       O,E              O,O      O,E   E,E    E,O ]
   *
   *
   * This class is used to implement the resulting linear operator
   *
   *      ~
   *      M  =  A(o,o) - D(o,e) . A^-1(e,e) . D(e,o)
   *
   * By construction, the linear operator is ONLY defined on the odd subset
   *
   */

  template<typename T, typename P, typename Q>
  class EvenOddPrecConstDetLinearOperatorArray : public EvenOddPrecLinearOperatorArray<T,P,Q>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~EvenOddPrecConstDetLinearOperatorArray() {}

    //! Only defined on the odd lattice
    const Subset& subset() const {return rb[1];}

    //! Expected length of array index
    virtual int size(void) const = 0;

    //! Return the fermion BC object for this linear operator
    virtual const FermBC<T,P,Q>& getFermBC() const = 0;

    //! Apply the operator onto a source vector
    /*! User should make sure deriv routines do a resize  */
    virtual void deriv(P& ds_u, const multi1d<T>& chi, const multi1d<T>& psi, 
		       enum PlusMinus isign) const
    {
      // Need deriv of  (A_oo - D_oe*Ainv_ee*D_eo*psi_e)
      enum PlusMinus msign = (isign == PLUS) ? MINUS : PLUS;

      //
      // Make sure the deriv routines do a resize !!!
      //
      multi1d<T>   tmp1(size()), tmp2(size()), tmp3(size());  // if an array is used here, the space is not reserved
      moveToFastMemoryHint(tmp1);
      moveToFastMemoryHint(tmp2);
      //      moveToFastMemoryHint(tmp3);

      P            ds_1;  // deriv routines should resize

      //
      // NOTE: even with even-odd decomposition, the ds_u will still have contributions
      // on all cb. So, no adding of ds_1 onto ds_u under a subset
      //
      //  ds_u  =  chi^dag * A'_oo * psi
      //      derivOddOddLinOp(ds_u, chi, psi, isign);

      


      //  ds_u  -=  chi^dag * D'_oe * Ainv_ee * D_eo * psi_o
      this->evenOddLinOp(tmp1, psi, isign);
      this->evenEvenInvLinOp(tmp2, tmp1, isign);
      this->derivOddEvenLinOp(ds_u, chi, tmp2, isign);

      //  ds_u  -=  chi^dag * D_oe * Ainv_ee * D'_eo * psi_o
      this->evenOddLinOp(tmp1, chi, msign);
      this->evenEvenInvLinOp(tmp2, tmp1, msign);
      this->derivEvenOddLinOp(ds_1, tmp2, psi, isign);
      ds_u += ds_1;

      // This is yucky and ws should have a *= function 
      // so that I wouldn't have to expose the fact that I expect
      // ds_u to be an Nd dimensional array
      for(int mu=0; mu < Nd; mu++) { 
	ds_u[mu] *= Real(-1);
      }

      getFermBC().zero(ds_u);
    }

    //! Apply the even-even block onto a source vector
    virtual void derivEvenEvenLinOp(P& ds_u, const multi1d<T>& chi, const multi1d<T>& psi, 
				    enum PlusMinus isign) const
    {
      QDPIO::cerr << "EvenOdd: not implemented" << endl;
      QDP_abort(1);
    }
  
    //! Apply the the even-odd block onto a source vector
    virtual void derivEvenOddLinOp(P& ds_u, const multi1d<T>& chi, const multi1d<T>& psi, 
				   enum PlusMinus isign) const
    {
      QDPIO::cerr << "EvenOdd: not implemented" << endl;
      QDP_abort(1);
    }
 
    //! Apply the the odd-even block onto a source vector
    virtual void derivOddEvenLinOp(P& ds_u, const multi1d<T>& chi, const multi1d<T>& psi, 
				   enum PlusMinus isign) const
    {
      QDPIO::cerr << "EvenOdd: not implemented" << endl;
      QDP_abort(1);
    }

    //! Apply the the odd-odd block onto a source vector
    virtual void derivOddOddLinOp(P& ds_u, const multi1d<T>& chi, const multi1d<T>& psi, 
				  enum PlusMinus isign) const
    {
      QDPIO::cerr << "EvenOdd: not implemented" << endl;
      QDP_abort(1);
    }
  };




}
#endif
