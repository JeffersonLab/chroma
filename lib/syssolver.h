// -*- C++ -*-
// $Id: syssolver.h,v 3.3 2007-02-22 21:11:45 bjoo Exp $
/*! @file
 * @brief Linear system solvers
 */

#ifndef __syssolver_h__
#define __syssolver_h__

#include "chromabase.h"

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

    //! Apply the operator onto a source vector
    /*! 
     * Solves   A*psi = chi  or  psi = A^(-1)*chi up to some accuracy.
     * There is the interesting possibility of generalizing to support PLUS/MINUS
     *
     * Should the accuracy be specified here ???
     */
    virtual SystemSolverResults_t operator() (T& psi, const T& chi) const = 0;

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

    //! Apply the operator onto a source vector
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

    //! Apply the operator onto a source vector
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

    //! Apply the operator onto a source vector
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

}



#endif
