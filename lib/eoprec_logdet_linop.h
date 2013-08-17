// -*- C++ -*-
// $Id: eoprec_logdet_linop.h,v 3.2 2007-02-22 21:11:45 bjoo Exp $
/*! @file
 * @brief Preconditioned  Linear Operators where the Even Even part depends on the gauge field.
 *
 * We assume we can evaluate Log Det E, where E is the Even Even part. 
 * Essentially this is for things like clover.
 */

#ifndef __eoprec_logdet_linop_h__
#define __eoprec_logdet_linop_h__

#include "eoprec_linop.h"

using namespace QDP::Hints;

namespace Chroma 
{

  //-------------------------------------------------------------------------
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
   * where A^{-1}_{ee} does depend on the gauge fields but we can 
   * exactly simulate the determinant without pseudofermion fields
   * since we know det A_{ee} and can write it in the action as 
   * exp( log det A ) = exp( Tr Ln A )
   * 
   * Since we can explicitly evaluate Tr Ln A for the action, we
   * can also evaluate the force contribution
   * 
   * d/dt ( Tr Ln A ) = Tr ( (d/dt)A A^{-1} )
   *
   * and  d/dt ( A^{-1} ) = A^{-1} [ (d/dt)A ] A^{-1}
   * 
   * hence we have functions  logDetEvenEvenLinOp()
   *  and                     derivEvenEvenLinOp()
   *  and                     derivLogDetEvenEvenLinOp()
   */

  template<typename T, typename P, typename Q>
  class EvenOddPrecLogDetLinearOperator : public EvenOddPrecLinearOperator<T,P,Q>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~EvenOddPrecLogDetLinearOperator() {}

    //! Only defined on the odd lattice
    const Subset& subset() const {return rb[1];}

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
      T   tmp1, tmp2, tmp3;  // if an array is used here, the space is not reserved
      moveToFastMemoryHint(tmp1);
      moveToFastMemoryHint(tmp2);
      moveToFastMemoryHint(tmp3);

      P   ds_1;  // deriv routines should resize

      //
      // NOTE: even with even-odd decomposition, the ds_u will still have contributions
      // on all cb. So, no adding of ds_1 onto ds_u under a subset
      //
      //  ds_u  =  chi^dag * A'_oo * psi
      this->derivOddOddLinOp(ds_u, chi, psi, isign);

      //  ds_u  -=  chi^dag * D'_oe * Ainv_ee * D_eo * psi_o
      this->evenOddLinOp(tmp1, psi, isign);
      this->evenEvenInvLinOp(tmp2, tmp1, isign);
      this->derivOddEvenLinOp(ds_1, chi, tmp2, isign);
      ds_u -= ds_1;

      //  ds_u  +=  chi^dag * D_oe * Ainv_ee * A'_ee * Ainv_ee * D_eo * psi_o

      // Reuse tmp2 = Ainv_ee D_eo psi_o
      this->evenOddLinOp(tmp1, chi, msign);
      this->evenEvenInvLinOp(tmp3, tmp1, msign);
      this->derivEvenEvenLinOp(ds_1, tmp3, tmp2, isign);
      ds_u += ds_1;

      //  ds_u  -=  chi^dag * D_oe * Ainv_ee * D'_eo * psi_o
      // Reuse tmp3^dag = chi^\dag D_oe Ainv_ee
      this->derivEvenOddLinOp(ds_1, tmp3, psi, isign);
      ds_u -= ds_1;

      getFermBC().zero(ds_u);
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
      T   tmp1;
      multi1d<T> tmp2(chi.size());
      multi1d<T> tmp3(chi.size());  // if an array is used here, the space is not reserved

      P   ds_1;  // deriv routines should resize

      //
      // NOTE: even with even-odd decomposition, the ds_u will still have contributions
      // on all cb. So, no adding of ds_1 onto ds_u under a subset
      //
      //  ds_u  =  chi^dag * A'_oo * psi
      this->derivOddOddLinOpMP(ds_u, chi, psi, isign);

      //  ds_u  -=  chi^dag * D'_oe * Ainv_ee * D_eo * psi_o
      for(int i=0; i < chi.size(); i++) { 
	this->evenOddLinOp(tmp1, psi[i], isign);
	this->evenEvenInvLinOp(tmp2[i], tmp1, isign);
      }
      this->derivOddEvenLinOpMP(ds_1, chi, tmp2, isign);
      ds_u -= ds_1;


      //  ds_u  +=  chi^dag * D_oe * Ainv_ee * A'_ee * Ainv_ee * D_eo * psi_o
      for(int i=0; i < chi.size(); i++) { 
	// Reuse tmp2 = Ainv_ee D_eo psi_o
	this->evenOddLinOp(tmp1, chi[i], msign);
	this->evenEvenInvLinOp(tmp3[i], tmp1, msign);
      }
      this->derivEvenEvenLinOpMP(ds_1, tmp3, tmp2, isign);
      ds_u += ds_1;

      //  ds_u  -=  chi^dag * D_oe * Ainv_ee * D'_eo * psi_o
      // Reuse tmp3^dag = chi^\dag D_oe Ainv_ee
      this->derivEvenOddLinOpMP(ds_1, tmp3, psi, isign);
      ds_u -= ds_1;

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
     ds_u.resize(Nd);
     ds_u = zero;
     
     P F_tmp; // deriv will resize
     for(int i=0; i < chi.size(); i++) { 
       this->derivEvenEvenLinOp(F_tmp, chi[i], psi[i], isign);
       ds_u += F_tmp;
     }
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
       this->derivOddEvenLinOp(F_tmp, chi[i], psi[i], isign);
       ds_u += F_tmp;
     }
   }

    virtual void derivOddOddLinOpMP(P& ds_u, const multi1d<T>& chi, const multi1d<T>& psi, 
				       enum PlusMinus isign) const 
   {
     // Trivially zero since we are constdet?
     ds_u.resize(Nd);
     ds_u = zero;
     
     P F_tmp;
     for(int i=0; i < chi.size(); i++) { 
       this->derivOddOddLinOp(F_tmp, chi[i], psi[i], isign);
       ds_u += F_tmp;
     }
   }

    //! Get the force from the EvenEven Trace Log
    virtual void derivLogDetEvenEvenLinOp(P& ds_u, enum PlusMinus isign) const
    {
      QDPIO::cerr << "EvenOdd: not implemented" << endl;
      QDP_abort(1);
    }

    //! Get the log det of the even even part
    virtual Double logDetEvenEvenLinOp(void) const = 0;
  };


  //-------------------------------------------------------------------------
  //! Even-odd preconditioned 5D linear operator
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
   * where A^{-1}_{ee} does depend on the gauge fields but we can 
   * exactly simulate the determinant without pseudofermion fields
   * since we know det A_{ee} and can write it in the action as 
   * exp( log det A ) = exp( Tr Ln A )
   * 
   * Since we can explicitly evaluate Tr Ln A for the action, we
   * can also evaluate the force contribution
   * 
   * d/dt ( Tr Ln A ) = Tr ( (d/dt)A A^{-1} )
   *
   * and  d/dt ( A^{-1} ) = A^{-1} [ (d/dt)A ] A^{-1}
   * 
   * hence we have functions  lnDetEvenEven()
   *  and                     derivEvenEven()
   *  and                     derivEvenEvenInv()
   *  --- this needs work, may not even need more than lnDetEvenEven
   *
   *
   */

  template<typename T, typename P, typename Q>
  class EvenOddPrecLogDetLinearOperatorArray : public EvenOddPrecLinearOperatorArray<T,P,Q>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~EvenOddPrecLogDetLinearOperatorArray() {}

    //! Only defined on the odd lattice
    const Subset& subset() const {return rb[1];}

    //! Get the szie expected of arrays
    virtual int size() const = 0;

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
      moveToFastMemoryHint(tmp3);

      P            ds_1;  // deriv routines should resize

      //
      // NOTE: even with even-odd decomposition, the ds_u will still have contributions
      // on all cb. So, no adding of ds_1 onto ds_u under a subset
      //
      //  ds_u  =  chi^dag * A'_oo * psi
      derivOddOddLinOp(ds_u, chi, psi, isign);

      //  ds_u  -=  chi^dag * D'_oe * Ainv_ee * D_eo * psi_o
      evenOddLinOp(tmp1, psi, isign);
      evenEvenInvLinOp(tmp2, tmp1, isign);
      derivOddEvenLinOp(ds_1, chi, tmp2, isign);
      ds_u -= ds_1;

      //  ds_u  +=  chi^dag * D_oe * Ainv_ee * A'_ee * Ainv_ee * D_eo * psi_o
      evenOddLinOp(tmp1, psi, isign);
      evenEvenInvLinOp(tmp2, tmp1, isign);
      evenOddLinOp(tmp1, chi, msign);
      evenEvenInvLinOp(tmp3, tmp1, msign);
      derivEvenEvenLinOp(ds_1, tmp3, tmp2, isign);
      ds_u += ds_1;

      //  ds_u  -=  chi^dag * D_oe * Ainv_ee * D'_eo * psi_o
      evenOddLinOp(tmp1, chi, msign);
      evenEvenInvLinOp(tmp3, tmp1, msign);
      derivEvenOddLinOp(ds_1, tmp3, psi, isign);
      ds_u -= ds_1;

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

    //! Get the force from the EvenEven Trace Log
    virtual void derivLogDetEvenEvenLinOp(P& ds_u, enum PlusMinus isign) const
    {
      QDPIO::cerr << "EvenOdd: not implemented" << endl;
      QDP_abort(1);
    }

    //! Get the log det of the even even part
    virtual Double logDetEvenEvenLinOp(void) const = 0;
  };



}
#endif
