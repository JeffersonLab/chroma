// -*- C++ -*-
/*! @file
 * @brief Preconditioned 4D Linop with Gauge Independent Even-Even part
 */

#ifndef __seoprec_constdet_linop_h__
#define __seoprec_constdet_linop_h__

#include "seoprec_linop.h"

using namespace QDP::Hints;

namespace Chroma 
{

  //----------------------------------------------------------------
  //! Symmetric even-odd preconditioned linear operator
  /*! @ingroup linop
   *
   * Support for symmetric even-odd preconditioned linear operators
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
   *      [         1                  0       ]
   *      [          E,E                E,O    ]
   *  L = [                                    ]
   *      [  A^(-1) D     A^(-1)       A       ]
   *      [   O,O    O,E   E,E          O,O    ]
   *
   * and
   *
   *      [      1             A^{-1}  D       ]
   *      [       E,E           EE      E,O    ]
   *  U = [                                    ]
   *      [      0                     1       ]
   *      [       O,E                   O,O    ]
   *
   * The preconditioned matrix is formed from
   *
   *  M'   =  L^-1 * M * U^-1
   *
   * where
   *
   *           [         1                  0       ]
   *           [          E,E                E,O    ]
   *  L^(-1) = [                                    ]
   *           [ -A^(-1) D     A^(-1)       A^(-1)  ]
   *           [   O,O    O,E   E,E          O,O    ]
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
   *      [      A                0                                 ]
   *      [       E,E              O,O                              ]
   *      [                                                         ]
   *      [      0                1    -  A^(-1)  D    A^(-1)  D    ]
   *      [       O,E              O,O     O,O     O,E   E,E    E,O ]
   *
   *
   * This class is used to implement the resulting linear operator
   *
   *      ~
   *      M  =  I(o,o)  -  A^-1(o,o) . D(o,e) . A^-1(e,e) . D(e,o)
   *
   * where A^{-1}_{ee} and A^{-1}_(o,o) are independent of the gauge fields. 
   * This means that their force term contributions are zero.
   *
   * This structure suits most of the linear operators we use, and 
   * It simplifies the force term.
   * By construction, the linear operator is ONLY defined on the odd subset
   *
   */

  template<typename T, typename P, typename Q>
  class SymEvenOddPrecConstDetLinearOperator : public EvenOddPrecLinearOperator<T,P,Q>  
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~SymEvenOddPrecConstDetLinearOperator() {}

    //! Return the fermion BC object for this linear operator
    virtual const FermBC<T,P,Q>& getFermBC() const = 0;

    //! Apply the derivative of the operator onto a source std::vector
    /*! User should make sure deriv routines do a resize  */
    virtual void deriv(P& ds_u, const T& chi, const T& psi, 
		       enum PlusMinus isign) const
    {
      // Need deriv of  (I_oo  -  Ainv_oo*D_oe*Ainv_ee*D_eo*psi_e)
      enum PlusMinus msign = (isign == PLUS) ? MINUS : PLUS;

      //
      // Make sure the deriv routines do a resize !!!
      //
      ds_u.resize(Nd);
      ds_u = zero;

      T   tmp1, Dpsi, ADpsi, Achi, DAchi, ADAchi;
      P   ds_1;  // deriv routines should resize

      //
      // NOTE: even with even-odd decomposition, the ds_u will still have contributions
      // on all cb. So, no adding of ds_1 onto ds_u under a subset
      //
      // Dpsi  = D_eo psi_o
      this->evenOddLinOp(Dpsi, psi, isign);

      // ADpsi = Ainv_ee D_eo psi_o
      this->evenEvenInvLinOp(ADpsi, Dpsi, isign);

      // Achi = Ainv_oo chi_o
      this->oddOddInvLinOp(Achi, chi, msign);

      // DAchi = D_eo Ainv_oo chi_o
      this->evenOddLinOp(DAchi, Achi, msign);

      // ADAchi  = Ainv_ee D_eo Ainv_oo chi_o
      this->evenEvenInvLinOp(ADAchi, DAchi, msign);
      
      //  ds_u  -=  chi^dag * Ainv_oo * D'_oe * Ainv_ee * D_eo * psi_o
      this->derivOddEvenLinOp(ds_1, Achi, ADpsi, isign);
      ds_u -= ds_1;

      //  ds_u  -=  chi^dag * Ainv_oo * D_oe * Ainv_ee * D'_eo * psi_o
      this->derivEvenOddLinOp(ds_1, ADAchi, psi, isign);
      ds_u -= ds_1;

      getFermBC().zero(ds_u);
    }

    //! Apply the even-even block onto a source std::vector
    virtual void derivEvenEvenLinOp(P& ds_u, const T& chi, const T& psi, 
				    enum PlusMinus isign) const
    {
      QDPIO::cerr << "EvenOdd: not implemented" << std::endl;
      QDP_abort(1);
    }
  
    //! Apply the the even-odd block onto a source std::vector
    virtual void derivEvenOddLinOp(P& ds_u, const T& chi, const T& psi, 
				   enum PlusMinus isign) const
    {
      QDPIO::cerr << "EvenOdd: not implemented" << std::endl;
      QDP_abort(1);
    }
 
    //! Apply the the odd-even block onto a source std::vector
    virtual void derivOddEvenLinOp(P& ds_u, const T& chi, const T& psi, 
				   enum PlusMinus isign) const
    {
      QDPIO::cerr << "EvenOdd: not implemented" << std::endl;
      QDP_abort(1);
    }

    //! Apply the the odd-odd block onto a source std::vector
    virtual void derivOddOddLinOp(P& ds_u, const T& chi, const T& psi, 
				  enum PlusMinus isign) const
    {
      QDPIO::cerr << "EvenOdd: not implemented" << std::endl;
      QDP_abort(1);
    }

    //! Apply the derivative of the operator onto a source std::vector
    /*! User should make sure deriv routines do a resize  */
    virtual void derivMultipole(P& ds_u, const multi1d<T>& chi, const multi1d<T>& psi, 
				enum PlusMinus isign) const
    {
      // Need deriv of  (I_oo  -  Ainv_oo*D_oe*Ainv_ee*D_eo*psi_e)
      enum PlusMinus msign = (isign == PLUS) ? MINUS : PLUS;

      //
      // Make sure the deriv routines do a resize !!!
      //
      ds_u.resize(Nd);
      ds_u = zero;

      T   tmp1, Dpsi, ADpsi, Achi, DAchi, ADAchi;
      P   ds_1;  // deriv routines should resize

      //
      // NOTE: even with even-odd decomposition, the ds_u will still have contributions
      // on all cb. So, no adding of ds_1 onto ds_u under a subset
      //
      for(int i=0; i < chi.size(); i++) { 
	// Dpsi  = D_eo psi_o
	this->evenOddLinOp(Dpsi, psi[i], isign);

	// ADpsi = Ainv_ee D_eo psi_o
	this->evenEvenInvLinOp(ADpsi, Dpsi, isign);

	// Achi = Ainv_oo chi_o
	this->oddOddInvLinOp(Achi, chi[i], msign);

	// DAchi = D_eo Ainv_oo chi_o
	this->evenOddLinOp(DAchi, Achi, msign);

	// ADAchi  = Ainv_ee D_eo Ainv_oo chi_o
	this->evenEvenInvLinOp(ADAchi, DAchi, msign);
      
	//  ds_u  -=  chi^dag * Ainv_oo * D'_oe * Ainv_ee * D_eo * psi_o
	this->derivOddEvenLinOp(ds_1, Achi, ADpsi, isign);
	ds_u -= ds_1;

	//  ds_u  -=  chi^dag * Ainv_oo * D_oe * Ainv_ee * D'_eo * psi_o
	this->derivEvenOddLinOp(ds_1, ADAchi, psi[i], isign);
	ds_u -= ds_1;
      }

      getFermBC().zero(ds_u);
    }

  };



  //----------------------------------------------------------------
  //! Symmetric even-odd preconditioned linear operator
  /*! @ingroup linop
   *
   * Support for symmetric even-odd preconditioned linear operators
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
   *      [         1                  0       ]
   *      [          E,E                E,O    ]
   *  L = [                                    ]
   *      [  A^(-1) D     A^(-1)       A       ]
   *      [   O,O    O,E   E,E          O,O    ]
   *
   * and
   *
   *      [      1             A^{-1}  D       ]
   *      [       E,E           EE      E,O    ]
   *  U = [                                    ]
   *      [      0                     1       ]
   *      [       O,E                   O,O    ]
   *
   * The preconditioned matrix is formed from
   *
   *  M'   =  L^-1 * M * U^-1
   *
   * where
   *
   *           [         1                  0       ]
   *           [          E,E                E,O    ]
   *  L^(-1) = [                                    ]
   *           [ -A^(-1) D     A^(-1)       A^(-1)  ]
   *           [   O,O    O,E   E,E          O,O    ]
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
   *      [      A                0                                 ]
   *      [       E,E              O,O                              ]
   *      [                                                         ]
   *      [      0                1    -  A^(-1)  D    A^(-1)  D    ]
   *      [       O,E              O,O     O,O     O,E   E,E    E,O ]
   *
   *
   * This class is used to implement the resulting linear operator
   *
   *      ~
   *      M  =  I(o,o)  -  A^-1(o,o) . D(o,e) . A^-1(e,e) . D(e,o)
   *
   * where A^{-1}_{ee} and A^{-1}_(o,o) are independent of the gauge fields. 
   * This means that their force term contributions are zero.
   *
   * This structure suits most of the linear operators we use, and 
   * It simplifies the force term.
   * By construction, the linear operator is ONLY defined on the odd subset
   *
   */

  template<typename T, typename P, typename Q>
  class SymEvenOddPrecConstDetLinearOperatorArray : public SymEvenOddPrecLinearOperatorArray<T,P,Q>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~SymEvenOddPrecConstDetLinearOperatorArray() {}

    //! Only defined on the odd lattice
    const Subset& subset() const {return rb[1];}

    //! Expected length of array index
    virtual int size(void) const = 0;

    //! Return the fermion BC object for this linear operator
    virtual const FermBC<T,P,Q>& getFermBC() const = 0;

    //! Apply the operator onto a source std::vector
    /*! User should make sure deriv routines do a resize  */
    virtual void deriv(P& ds_u, const multi1d<T>& chi, const multi1d<T>& psi, enum PlusMinus isign) const
    {
      QDPIO::cerr << "deriv: not implemented" << std::endl;
      QDP_abort(1);
    }

    //! Apply the even-even block onto a source std::vector
    virtual void derivEvenEvenLinOp(P& ds_u, const multi1d<T>& chi, const multi1d<T>& psi, 
				    enum PlusMinus isign) const
    {
      QDPIO::cerr << "EvenOdd: not implemented" << std::endl;
      QDP_abort(1);
    }
  
    //! Apply the the even-odd block onto a source std::vector
    virtual void derivEvenOddLinOp(P& ds_u, const multi1d<T>& chi, const multi1d<T>& psi, 
				   enum PlusMinus isign) const
    {
      QDPIO::cerr << "EvenOdd: not implemented" << std::endl;
      QDP_abort(1);
    }
 
    //! Apply the the odd-even block onto a source std::vector
    virtual void derivOddEvenLinOp(P& ds_u, const multi1d<T>& chi, const multi1d<T>& psi, 
				   enum PlusMinus isign) const
    {
      QDPIO::cerr << "EvenOdd: not implemented" << std::endl;
      QDP_abort(1);
    }

    //! Apply the the odd-odd block onto a source std::vector
    virtual void derivOddOddLinOp(P& ds_u, const multi1d<T>& chi, const multi1d<T>& psi, 
				  enum PlusMinus isign) const
    {
      QDPIO::cerr << "EvenOdd: not implemented" << std::endl;
      QDP_abort(1);
    }
  };

}
#endif
