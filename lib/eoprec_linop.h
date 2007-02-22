// -*- C++ -*-
/*! @file
 * @brief Base class for even-odd preconditioned 4D and 5D Linop
 */

#ifndef __eoprec_linop_h__
#define __eoprec_linop_h__

#include "chromabase.h"
#include "linearop.h"

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
  class EvenOddPrecLinearOperator : public DiffLinearOperator<T,P,Q>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~EvenOddPrecLinearOperator() {}

    //! Only defined on the odd lattice
    const Subset& subset() const {return rb[1];}

    //! Return the fermion BC object for this linear operator
    virtual const FermBC<T,P,Q>& getFermBC() const = 0;

    //! Apply the even-even block onto a source vector
    /*! This does not need to be optimized */
    virtual void evenEvenLinOp(T& chi, const T& psi, 
			       enum PlusMinus isign) const = 0;
  
    //! Apply the inverse of the even-even block onto a source vector
    virtual void evenEvenInvLinOp(T& chi, const T& psi, 
				  enum PlusMinus isign) const = 0;
  
    //! Apply the the even-odd block onto a source vector
    virtual void evenOddLinOp(T& chi, const T& psi, 
			      enum PlusMinus isign) const = 0;

    //! Apply the the odd-even block onto a source vector
    virtual void oddEvenLinOp(T& chi, const T& psi, 
			      enum PlusMinus isign) const = 0;

    //! Apply the the odd-odd block onto a source vector
    virtual void oddOddLinOp(T& chi, const T& psi, 
			     enum PlusMinus isign) const = 0;

    //! Apply the operator onto a source vector
    virtual void operator() (T& chi, const T& psi, 
			     enum PlusMinus isign) const
    {
      T   tmp1, tmp2; moveToFastMemoryHint(tmp1); moveToFastMemoryHint(tmp2);

      /*  Tmp1   =  D     A^(-1)     D    Psi  */
      /*      O      O,E        E,E   E,O    O */
      evenOddLinOp(tmp1, psi, isign);
      evenEvenInvLinOp(tmp2, tmp1, isign);
      oddEvenLinOp(tmp1, tmp2, isign);

      /*  Chi   =  A    Psi  -  Tmp1  */
      /*     O      O,O    O        O */
      oddOddLinOp(chi, psi, isign);
      chi[rb[1]] -= tmp1;

      getFermBC().modifyF(chi, rb[1]);
    }


    //! Apply the UNPRECONDITIONED operator onto a source vector
    /*! Mainly intended for debugging */
    virtual void unprecLinOp(T& chi, const T& psi, 
			     enum PlusMinus isign) const
    {
      T   tmp1, tmp2;  moveToFastMemoryHint(tmp1); moveToFastMemoryHint(tmp2);

      /*  Chi   =   A    Psi   +    D    Psi   */
      /*     E       E,E    O        E,O    O  */
      evenEvenLinOp(tmp1, psi, isign);
      evenOddLinOp(tmp2, psi, isign);
      chi[rb[0]] = tmp1 + tmp2;

      /*  Chi   =  A    Psi  -  Tmp1  */
      /*     O      O,O    O        O */
      oddEvenLinOp(tmp1, psi, isign);
      oddOddLinOp(tmp2, psi, isign);
      chi[rb[1]] = tmp1 + tmp2;

      getFermBC().modifyF(chi);
    }


    //! Apply the even-even block onto a source vector
    virtual void derivEvenEvenLinOp(P& ds_u, const T& chi, const T& psi, 
				    enum PlusMinus isign) const = 0;
   
    //! Apply the the even-odd block onto a source vector
    virtual void derivEvenOddLinOp(P& ds_u, const T& chi, const T& psi, 
				   enum PlusMinus isign) const = 0;
 
    //! Apply the the odd-even block onto a source vector
    virtual void derivOddEvenLinOp(P& ds_u, const T& chi, const T& psi, 
				   enum PlusMinus isign) const = 0;

    //! Apply the the odd-odd block onto a source vector
    virtual void derivOddOddLinOp(P& ds_u, const T& chi, const T& psi, 
				  enum PlusMinus isign) const = 0;

    /*! Mainly intended for debugging */
    virtual void derivUnprecLinOp(P& ds_u, const T& chi, const T& psi, 
				  enum PlusMinus isign) const
    {
      T   tmp1, tmp2;  moveToFastMemoryHint(tmp1); moveToFastMemoryHint(tmp2);
      P   ds_tmp;  // deriv routines should resize

      //  ds_u  =  chi^dag * A'_ee * psi
      derivEvenEvenLinOp(ds_u, chi, psi, isign);

      //  ds_u +=  chi^dag * D'_eo * psi
      derivEvenOddLinOp(ds_tmp, chi, psi, isign);
      ds_u += ds_tmp;

      //  ds_u +=  chi^dag * D'_oe * psi
      derivOddEvenLinOp(ds_tmp, chi, psi, isign);
      ds_u += ds_tmp;

      //  ds_u +=  chi^dag * A'_oo * psi
      derivOddOddLinOp(ds_tmp, chi, psi, isign);
      ds_u += ds_tmp;

      getFermBC().zero(ds_u);
    }


    //! Apply the derivative of the operator onto a source vector
    /*! User should make sure deriv routines do a resize.
     *  This function is left pure virtual - as derived 
     *  functions need to override it 
     */
    virtual void deriv(P& ds_u, const T& chi, const T& psi, 
		       enum PlusMinus isign) const = 0;


    //! Return flops performed by the evenEvenLinOp
    virtual unsigned long evenEvenNFlops() const { return 0; }
    
    //! Return flops performed by the evenOddLinOp
    virtual unsigned long evenOddNFlops() const { return 0; }

    //! Return flops performed by the oddEvenLinOp
    virtual unsigned long oddEvenNFlops() const { return 0; }

    //! Return flops performed by the oddOddLinOp
    virtual unsigned long oddOddNFlops() const { return 0; }

    //! Return flops performed by the evenEvenInvLinOp
    virtual unsigned long evenEvenInvNFlops() const { return 0; }

    //! Return flops performed by the operator()
    virtual unsigned long nFlops() const { 
      return 0;
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
  class EvenOddPrecLinearOperatorArray : public DiffLinearOperatorArray<T,P,Q>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~EvenOddPrecLinearOperatorArray() {}

    //! Only defined on the odd lattice
    const Subset& subset() const {return rb[1];}

    //! Expected length of array index
    virtual int size(void) const = 0;

    //! Return the fermion BC object for this linear operator
    virtual const FermBC<T,P,Q>& getFermBC() const = 0;

    //! Apply the even-even block onto a source vector
    /*! This does not need to be optimized */
    virtual void evenEvenLinOp(multi1d<T>& chi, const multi1d<T>& psi, 
			       enum PlusMinus isign) const = 0;
  
    //! Apply the inverse of the even-even block onto a source vector
    virtual void evenEvenInvLinOp(multi1d<T>& chi, const multi1d<T>& psi, 
				  enum PlusMinus isign) const = 0;
  
    //! Apply the the even-odd block onto a source vector
    virtual void evenOddLinOp(multi1d<T>& chi, const multi1d<T>& psi, 
			      enum PlusMinus isign) const = 0;

    //! Apply the the odd-even block onto a source vector
    virtual void oddEvenLinOp(multi1d<T>& chi, const multi1d<T>& psi, 
			      enum PlusMinus isign) const = 0;

    //! Apply the the odd-odd block onto a source vector
    virtual void oddOddLinOp(multi1d<T>& chi, const multi1d<T>& psi, 
			     enum PlusMinus isign) const = 0;

    //! Apply the operator onto a source vector
    virtual void operator() (multi1d<T>& chi, const multi1d<T>& psi, 
			     enum PlusMinus isign) const
    {
      multi1d<T>  tmp1(size());  moveToFastMemoryHint(tmp1);
      multi1d<T>  tmp2(size());  moveToFastMemoryHint(tmp2);

      /*  Tmp1   =  D     A^(-1)     D    Psi  */
      /*      O      O,E        E,E   E,O    O */
      evenOddLinOp(tmp1, psi, isign);
      evenEvenInvLinOp(tmp2, tmp1, isign);
      oddEvenLinOp(tmp1, tmp2, isign);

      /*  Chi   =  A    Psi  -  Tmp1  */
      /*     O      O,O    O        O */
      oddOddLinOp(chi, psi, isign);
      for(int n=0; n < size(); ++n)
	chi[n][rb[1]] -= tmp1[n];

      getFermBC().modifyF(chi, rb[1]);
    }


    //! Apply the UNPRECONDITIONED operator onto a source vector
    /*! Mainly intended for debugging */
    virtual void unprecLinOp(multi1d<T>& chi, const multi1d<T>& psi, 
			     enum PlusMinus isign) const
    {
      multi1d<T>  tmp1(size()); moveToFastMemoryHint(tmp1);
      multi1d<T>  tmp2(size()); moveToFastMemoryHint(tmp2);

      /*  Chi   =   A    Psi   +    D    Psi   */
      /*     E       E,E    O        E,O    O  */
      evenEvenLinOp(tmp1, psi, isign);
      evenOddLinOp(tmp2, psi, isign);
      for(int n=0; n < size(); ++n)
	chi[n][rb[0]] = tmp1[n] + tmp2[n];

      /*  Chi   =   D    Psi    +    A    Psi   */
      /*     O       O,E    E         O,O    O  */
      oddEvenLinOp(tmp1, psi, isign);
      oddOddLinOp(tmp2, psi, isign);
      for(int n=0; n < size(); ++n)
	chi[n][rb[1]] = tmp1[n] + tmp2[n];

      getFermBC().modifyF(chi);
    }


    //! Apply the even-even block onto a source vector
    virtual void derivEvenEvenLinOp(P& ds_u, const multi1d<T>& chi, const multi1d<T>& psi, 
				    enum PlusMinus isign) const = 0;
  
    //! Apply the the even-odd block onto a source vector
    virtual void derivEvenOddLinOp(P& ds_u, const multi1d<T>& chi, const multi1d<T>& psi, 
				   enum PlusMinus isign) const = 0;
 
    //! Apply the the odd-even block onto a source vector
    virtual void derivOddEvenLinOp(P& ds_u, const multi1d<T>& chi, const multi1d<T>& psi, 
				   enum PlusMinus isign) const = 0;

    //! Apply the the odd-odd block onto a source vector
    virtual void derivOddOddLinOp(P& ds_u, const multi1d<T>& chi, const multi1d<T>& psi, 
				  enum PlusMinus isign) const = 0;

    /*! Mainly intended for debugging */
    virtual void derivUnprecLinOp(P& ds_u, const multi1d<T>& chi, const multi1d<T>& psi, 
				  enum PlusMinus isign) const
    {
      T   tmp1, tmp2;  // if an array is used here, the space is not reserved
      moveToFastMemoryHint(tmp1);
      moveToFastMemoryHint(tmp2);

      P   ds_tmp;  // deriv routines should resize

      //  ds_u  =  chi^dag * A'_ee * psi
      derivEvenEvenLinOp(ds_u, chi, psi, isign);

      //  ds_u +=  chi^dag * D'_eo * psi
      derivEvenOddLinOp(ds_tmp, chi, psi, isign);
      ds_u += ds_tmp;

      //  ds_u +=  chi^dag * D'_oe * psi
      derivOddEvenLinOp(ds_tmp, chi, psi, isign);
      ds_u += ds_tmp;

      //  ds_u +=  chi^dag * A'_oo * psi
      derivOddOddLinOp(ds_tmp, chi, psi, isign);
      ds_u += ds_tmp;

      getFermBC().zero(ds_u);
    }

    //! Apply the operator onto a source vector
    /*! User should make sure deriv routines do a resize.
     *  This function is left pure virtual - as derived 
     *  functions need to override it 
     */
    virtual void deriv(P& ds_u, const multi1d<T>& chi, const multi1d<T>& psi, 
		       enum PlusMinus isign) const = 0;

    //! Apply the derivative of the UNPRECONDITIONED operator onto a source vector

    //! Return flops performed by the evenEvenLinOp
    virtual unsigned long evenEvenNFlops() const { return 0; }
    
    //! Return flops performed by the evenOddLinOp
    virtual unsigned long evenOddNFlops() const { return 0; }

    //! Return flops performed by the oddEvenLinOp
    virtual unsigned long oddEvenNFlops() const { return 0; }

    //! Return flops performed by the oddOddLinOp
    virtual unsigned long oddOddNFlops() const { return 0; }

    //! Return flops performed by the evenEvenInvLinOp
    virtual unsigned long evenEvenInvNFlops() const { return 0; }

    //! Return flops performed by the operator()
    virtual unsigned long nFlops() const { 
      return (this->oddOddNFlops()
	      +this->oddEvenNFlops()
	      +this->evenEvenInvNFlops()
	      +this->evenOddNFlops());
    }


  };


} // End namespace Chroma

#endif
