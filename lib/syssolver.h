// -*- C++ -*-
// $Id: syssolver.h,v 1.1 2005-01-03 15:36:50 edwards Exp $
/*! @file
 * @brief Linear system solvers
 */

#ifndef __syssolver_h__
#define __syssolver_h__

using namespace QDP;

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

using namespace Chroma;


#endif
