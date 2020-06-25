// -*- C++ -*-
/*! @file
 * @brief Linear system solvers
 */

#ifndef __syssolver_h__
#define __syssolver_h__

#include "chromabase.h"
#include <vector>
#include <memory>

namespace Chroma
{
  //-----------------------------------------------------------------------------------
  //! Holds return info from SystemSolver call
  /*! @ingroup solvers */
  struct SystemSolverResults_t
  {
    SystemSolverResults_t() {n_count=0; resid=zero;}

    int  n_count;      /*!< Number of iterations */
    Real resid;        /*!< (True) Residual of unpreconditioned problem, 
			*    resid = sqrt(norm2(rhs - A.soln)) */
  };


  //-----------------------------------------------------------------------------------
  //! Linear system solvers
  /*! @ingroup solvers
   *
   * Solves linear systems of equations. The solver may only live on a subset.
   */
  template<typename T>
  class SystemSolver
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~SystemSolver() {}

    //! Apply the operator onto a source std::vector
    /*! 
     * Solves   A*psi = chi  or  psi = A^(-1)*chi up to some accuracy.
     * There is the interesting possibility of generalizing to support PLUS/MINUS
     *
     * Should the accuracy be specified here ???
     */
    virtual SystemSolverResults_t operator() (T& psi, const T& chi) const = 0;

    virtual std::vector<SystemSolverResults_t> operator() (const std::vector<std::shared_ptr<T>>& psi, const std::vector<std::shared_ptr<const T>>& chi) const
    {
       QDPIO::cout << "Solving the linear systems one by one" << std::endl;
       assert(psi.size() == chi.size());
       std::vector<SystemSolverResults_t> res(psi.size());
       for(int ncols = psi.size(), col=0; col < ncols; ++col)
         res[col] = this->operator()(*psi[col], *chi[col]);
       return res;
    }

    //! Return the subset on which the operator acts
    virtual const Subset& subset() const = 0;
  };


  //-----------------------------------------------------------------------------------
  //! Projector, that is, P*P*x = P*x
  /*! @ingroup solvers
   *
   * Apply a projector. The projector may only live on a subset.
   */
  template<typename T>
  class Projector
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~Projector() {}

    //! Apply the oblique projector A*V*inv(U^H*A*V)*U^H
    /*! 
     * Returns A*V*inv(U^H*A*V)*U^H*chi = psi
     */
    virtual void AVUObliqueProjector(T& psi, const T& chi) const = 0;

    //! Apply the oblique projector V*inv(U^H*A*V)*U^H*A
    /*! 
     * Returns V*inv(V^H*A*V)*U^H*A*chi = psi
     */
    virtual void VUAObliqueProjector(T& psi, const T& chi) const = 0;

    //! Rank of the projector, which is the rank of U and V also
    virtual unsigned int rank() const = 0;

    //! Return U[i]
    virtual void U(unsigned int i, T& psi) const = 0;

     //! Return V[i]
    virtual void V(unsigned int i, T& psi) const = 0;

    //! Return U[i]^H*A*V[i]
    virtual void lambda(unsigned int i, DComplex& lambda) const = 0;

    //! Return the subset on which the operator acts
    virtual const Subset& subset() const = 0;
  };



  //-----------------------------------------------------------------------------------
  //! Linear system solvers of arrays
  /*! @ingroup solvers
   *
   * Solves linear systems of equations. The solver may only live on a subset.
   */
  template<typename T>
  class SystemSolverArray
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~SystemSolverArray() {}

    //! Expected length of array index
    virtual int size() const = 0;

    //! Apply the operator onto a source std::vector
    /*! 
     * Solves   A*psi = chi  or  psi = A^(-1)*chi up to some accuracy.
     * There is the interesting possibility of generalizing to support PLUS/MINUS
     *
     * Should the accuracy be specified here ???
     */
    virtual SystemSolverResults_t operator() (multi1d<T>& psi, const multi1d<T>& chi) const = 0;

    //! Return the subset on which the operator acts
    virtual const Subset& subset() const = 0;
  };



  //-----------------------------------------------------------------------------------
  //! Linear multi-system solvers
  /*! @ingroup solvers
   *
   * Solves multi-shift linear systems of equations. The solver may only live on a subset.
   */
  template<typename T>
  class MultiSystemSolver
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~MultiSystemSolver() {}

    //! Apply the operator onto a source std::vector
    /*! 
     * Solves   A*psi = chi  or  psi = A^(-1)*chi up to some accuracy.
     * There is the interesting possibility of generalizing to support PLUS/MINUS
     *
     * Should the accuracy be specified here ???
     *
     * Should the shifts be here or in the constructor???
     */
    virtual SystemSolverResults_t operator() (multi1d<T>& psi, 
					      const multi1d<Real>& shifts, 
					      const T& chi) const = 0;

    //! Return the subset on which the operator acts
    virtual const Subset& subset() const = 0;
  };




  //-----------------------------------------------------------------------------------
  //! Linear multi-system solvers of arrays
  /*! @ingroup solvers
   *
   * Solves multi-shift linear systems of equations. The solver may only live on a subset.
   */
  template<typename T>
  class MultiSystemSolverArray
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~MultiSystemSolverArray() {}

    //! Expected length of array index
    virtual int size() const = 0;

    //! Apply the operator onto a source std::vector
    /*! 
     * Solves   A*psi = chi  or  psi = A^(-1)*chi up to some accuracy.
     * There is the interesting possibility of generalizing to support PLUS/MINUS
     *
     * Should the accuracy be specified here ???
     *
     * Should the shifts be here or in the constructor???
     */
    virtual SystemSolverResults_t operator() (multi1d< multi1d<T> >& psi, 
					      const multi1d<Real>& shifts, 
					      const multi1d<T>& chi) const = 0;

    //! Return the subset on which the operator acts
    virtual const Subset& subset() const = 0;
  };


  //-----------------------------------------------------------------------------------
  //! Linear multi-system solvers with accumulation
  /*! @ingroup solvers
   *
   * A tweak of solvers for multi-shift linear systems of equations. 
   * so that the results are accumulated a la:
   * 
   *   s = A sum_sigma t_i X_i 
   * where x_i are the solutions.
   * The solver may only live on a subset.
   */
  template<typename T>
  class MultiSystemSolverAccumulate
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~MultiSystemSolverAccumulate() {}

    //! Apply the operator onto a source std::vector
    virtual SystemSolverResults_t operator() (T& psi, 
					      const Real& norm,
					      const multi1d<Real>& residues,
					      const multi1d<Real>& poles, 
					      const T& chi) const = 0;

    //! Return the subset on which the operator acts
    virtual const Subset& subset() const = 0;
  };

  //-----------------------------------------------------------------------------------
  //! Linear multi-system solvers with accumulation
  /*! @ingroup solvers
   *
   * A tweak of solvers for multi-shift linear systems of equations. 
   * so that the results are accumulated a la:
   * 
   *   s = A sum_sigma t_i X_i 
   * where x_i are the solutions.
   * The solver may only live on a subset.
   */
  template<typename T>
  class MultiSystemSolverAccumulateArray
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~MultiSystemSolverAccumulateArray() {}

    //! Apply the operator onto a source std::vector
    virtual SystemSolverResults_t operator() (multi1d<T>& psi, 
					      const Real& norm,
					      const multi1d<Real>& residues,
					      const multi1d<Real>& poles, 
					      const multi1d<T>& chi) const = 0;

    //! Return the subset on which the operator acts
    virtual const Subset& subset() const = 0;
  };

}



#endif
