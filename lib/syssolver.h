// -*- C++ -*-
// $Id: syssolver.h,v 1.2 2005-01-14 20:13:04 edwards Exp $
/*! @file
 * @brief Linear system solvers
 */

#ifndef __syssolver_h__
#define __syssolver_h__


namespace Chroma
{
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
    virtual int operator() (T& psi, const T& chi) const = 0;

    //! Return the subset on which the operator acts
    virtual const OrderedSubset& subset() const = 0;
  };

}



#endif
