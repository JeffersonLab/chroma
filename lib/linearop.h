// -*- C++ -*-
// $Id: linearop.h,v 3.5 2009-04-17 02:05:30 bjoo Exp $

/*! @file
 * @brief Linear Operators
 */

#ifndef __linearop_h__
#define __linearop_h__

#include "chromabase.h"
#include "fermbc.h"

using namespace QDP::Hints;

namespace Chroma
{

  //-----------------------------------------------------------------------------------
  //! Linear Operator
  /*! @ingroup linop
   *
   * Supports creation and application for linear operators that
   * hold things like Dirac operators, etc.
   */
  template<typename T>
  class LinearOperator
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~LinearOperator() {}

    //! Apply the operator onto a source vector
    virtual void operator() (T& chi, const T& psi, enum PlusMinus isign) const = 0;

    //! Apply the operator onto a source vector to some precision
    virtual void operator() (T& chi, const T& psi, enum PlusMinus isign, 
			     Real epsilon) const
    {
      (*this)(chi,psi,isign);
    }

    //! Return the subset on which the operator acts
    virtual const Subset& subset() const = 0;


    //! Return the number of flops performed by operator()
    //! Base linop returns 0 and this can be overridden
    virtual unsigned long nFlops() const { return 0; }
  };


  //-----------------------------------------------------------------------------------
  //! Linear Operator to arrays
  /*! @ingroup linop
   *
   * Supports creation and application for linear operators that
   * hold things like Dirac operators, etc.
   */
  template<typename T>
  class LinearOperatorArray
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~LinearOperatorArray() {}

    //! Expected length of array index
    virtual int size() const = 0;

    //! Apply the operator onto a source vector
    virtual void operator() (multi1d<T>& chi, const multi1d<T>& psi, 
			     enum PlusMinus isign) const = 0;

    //! Apply the operator onto a source vector to some precision
    virtual void operator() (multi1d<T>& chi, const multi1d<T>& psi, 
			     enum PlusMinus isign, Real epsilon) const
    {
      (*this)(chi,psi,isign);
    }

    //! Return the subset on which the operator acts
    virtual const Subset& subset() const = 0;

    //! Return the number of flops performed by operator()
    //! Base linop returns 0 and this can be overridden
    virtual unsigned long nFlops() const { return 0; };
  };


  //-----------------------------------------------------------------------------------
  //! Differentiable Linear Operator
  /*! @ingroup linop
   *
   * Supports creation and application for linear operators that
   * hold things like Dirac operators, etc. that are differentiable
   */
  template<typename T, typename P, typename Q>
  class DiffLinearOperator : public LinearOperator<T>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~DiffLinearOperator() {}

    //! Return the fermion BC object for this linear operator
    virtual const FermBC<T,P,Q>& getFermBC() const = 0;

    //! Apply the derivative of the operator onto a source vector
    /*! Default implementation */
    virtual void deriv(P& ds_u, const T& chi, const T& psi, 
		       enum PlusMinus isign) const
    {
      QDPIO::cerr << "deriv: not implemented" << endl;
      QDP_abort(1);
    }



    //! Apply the derivative of the operator onto a source vector to some precision
    /*! Default implementation */
    virtual void deriv(P& ds_u, const T& chi, const T& psi, 
		       enum PlusMinus isign, const Real& epsilon) const
    {
      deriv(ds_u,chi,psi,isign);
    }

    //! Return the force for multiple poles
    virtual void derivMultipole(P& ds_u, const multi1d<T>& chi, const multi1d<T>&psi, enum PlusMinus isign) const
    {
      // Code up in terms of derivs for category default
      ds_u.resize(Nd);     
      ds_u = zero;

      P F_tmp;
      for(int i=0; i < chi.size(); i++) { 
	deriv(F_tmp, chi[i], psi[i], isign);
	ds_u += F_tmp;
      }
    }

  };


  //-----------------------------------------------------------------------------------
  //! Differentiable Linear Operator
  /*! @ingroup linop
   *
   * Supports creation and application for linear operators that
   * hold things like Dirac operators, etc. that are differentiable
   */
  template<typename T, typename P, typename Q>
  class DiffLinearOperatorArray : public LinearOperatorArray<T>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~DiffLinearOperatorArray() {}

    //! Return the fermion BC object for this linear operator
    virtual const FermBC<T,P,Q>& getFermBC() const = 0;

    //! Apply the derivative of the operator onto a source vector
    /*! Default implementation */
    virtual void deriv(P& ds_u, const multi1d<T>& chi, const multi1d<T>& psi, 
		       enum PlusMinus isign) const
    {
      QDPIO::cerr << "deriv: not implemented" << endl;
      QDP_abort(1);
    }

    //! Apply the derivative of the operator onto a source vector to some precision
    /*! Default implementation */
    virtual void deriv(P& ds_u, const multi1d<T>& chi, const multi1d<T>& psi, 
		       enum PlusMinus isign, const Real& epsilon) const
    {
      deriv(ds_u,chi,psi,isign);
    }
  };


  //----------------------------------------------------------------
  //! Unpreconditioned linear operator including derivatives
  /*! @ingroup linop
   *
   * Support for unpreconditioned linear operators with derivative
   */
  template<typename T, typename P, typename Q>
  class UnprecLinearOperator : public DiffLinearOperator<T,P,Q>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~UnprecLinearOperator() {}

    //! Only defined on the entire lattice
    const Subset& subset() const {return all;}
  };


  //----------------------------------------------------------------
  //! Unpreconditioned linear operator including derivatives
  /*! @ingroup linop
   *
   * Support for unpreconditioned linear operators with derivative
   */
  template<typename T, typename P, typename Q>
  class UnprecLinearOperatorArray : public DiffLinearOperatorArray<T,P,Q>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~UnprecLinearOperatorArray() {}

    //! Only defined on the entire lattice
    const Subset& subset() const {return all;}
  };


  //-----------------------------------------------------------------------------------
  //! Dslash-like Linear Operator
  /*! @ingroup linop
   *
   * These are concessions/optimizations for red-black checkboarding 
   */
  template<typename T, typename P, typename Q>
  class DslashLinearOperator : public DiffLinearOperator<T,P,Q>
  {
  public:
    //! Virtual destructor to help in cleanup
    virtual ~DslashLinearOperator() {}

    //! Return the fermion BC object for this linear operator
    virtual const FermBC<T,P,Q>& getFermBC() const = 0;

    //! Apply operator on both checkerboards (entire lattice)
    virtual void operator() (T& d, const T& psi, enum PlusMinus isign) const
    {
      apply(d, psi, isign, 0);
      apply(d, psi, isign, 1);
    }

    //! Apply checkerboarded linear operator
    /*! 
     * To avoid confusion (especially of the compilers!), call the checkerboarded
     * apply instead of operator()
     */
    virtual void apply (T& chi, const T& psi, enum PlusMinus isign, int cb) const = 0;


    //! Take deriv of D
    /*!
     * \param chi     left vector                                 (Read)
     * \param psi     right vector                                (Read)
     * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
     *
     * \return Computes   \f$chi^\dag * \dot(D} * psi\f$
     */
    virtual void deriv(P& ds_u, const T& chi, const T& psi, 
		       enum PlusMinus isign) const
    {
      QDPIO::cerr << "deriv: not implemented" << endl;
      QDP_abort(1);
    }

    //! Take deriv of D
    /*!
     * \param chi     left vector on cb                           (Read)
     * \param psi     right vector on 1-cb                        (Read)
     * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
     * \param cb      Checkerboard of chi vector                  (Read)
     *
     * \return Computes   \f$chi^\dag * \dot(D} * psi\f$
     */
    virtual void deriv(P& ds_u, const T& chi, const T& psi, 
		       enum PlusMinus isign, int cb) const
    {
      QDPIO::cerr << "deriv: not implemented" << endl;
      QDP_abort(1);
    }

    virtual void derivMultipole(P& ds_u, const multi1d<T>& chi, const multi1d<T>&psi, enum PlusMinus isign) const
    {
      // Code up in terms of derivs for category default
      ds_u.resize(Nd);
      ds_u = zero;

      P F_tmp;
      for(int i=0; i < chi.size(); i++) { 
	deriv(F_tmp, chi[i], psi[i], isign);
	ds_u += F_tmp;
      }

    }


    virtual void derivMultipole(P& ds_u, const multi1d<T>& chi, const multi1d<T>&psi, enum PlusMinus isign,int cb) const
    {
      // Code up in terms of derivs for category default
      ds_u.resize(Nd);
      ds_u = zero;

      P F_tmp;
      for(int i=0; i < chi.size(); i++) { 
	deriv(F_tmp, chi[i], psi[i], isign,cb);
	ds_u += F_tmp;
      }

    }
  };



  //-----------------------------------------------------------------------------------
  //! Dslash-like Linear Operator for arrays
  /*! @ingroup linop
   *
   * These are concessions/optimizations for red-black checkboarding 
   */
  template<typename T, typename P, typename Q>
  class DslashLinearOperatorArray : public DiffLinearOperatorArray<T,P,Q>
  {
  public:
    //! Virtual destructor to help in cleanup
    virtual ~DslashLinearOperatorArray() {}

    //! Return the fermion BC object for this linear operator
    virtual const FermBC<T,P,Q>& getFermBC() const = 0;

    //! Expected length of array index
    virtual int size() const = 0;

    //! Apply operator on both checkerboards (entire lattice)
    virtual void operator() (multi1d<T>& d, const multi1d<T>& psi, enum PlusMinus isign) const
    {
      d.resize(size());
      apply(d, psi, isign, 0);
      apply(d, psi, isign, 1);
    }

    //! Apply checkerboarded linear operator
    /*! 
     * To avoid confusion (especially of the compilers!), call the checkerboarded
     * apply instead of operator()
     */
    virtual void apply (multi1d<T>& chi, const multi1d<T>& psi, enum PlusMinus isign, int cb) const = 0;

    //! Take deriv of D
    /*!
     * \param chi     left vector                                 (Read)
     * \param psi     right vector                                (Read)
     * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
     *
     * \return Computes   \f$chi^\dag * \dot(D} * psi\f$
     */
    virtual void deriv(P& ds_u, const multi1d<T>& chi, const multi1d<T>& psi, 
		       enum PlusMinus isign) const
    {
      QDPIO::cerr << "deriv: not implemented" << endl;
      QDP_abort(1);
    }

    //! Take deriv of D
    /*!
     * \param chi     left vector on cb                           (Read)
     * \param psi     right vector on 1-cb                        (Read)
     * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
     * \param cb      Checkerboard of chi vector                  (Read)
     *
     * \return Computes   \f$chi^\dag * \dot(D} * psi\f$
     */
    virtual void deriv(P& ds_u, const multi1d<T>& chi, const multi1d<T>& psi, 
		       enum PlusMinus isign, int cb) const
    {
      QDPIO::cerr << "deriv: not implemented" << endl;
      QDP_abort(1);
    }
  };

}



#endif
