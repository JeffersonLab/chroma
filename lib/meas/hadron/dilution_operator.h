// -*- C++ -*-
// $Id: dilution_operator.h,v 1.1 2007-12-14 06:53:42 edwards Exp $
/*! \file
 *  \brief Dilution operators
 */

#ifndef __dilution_operator_h__
#define __dilution_operator_h__

#include "chromabase.h"
#include <iterator>

namespace Chroma
{
  //! Abstract dilution operator
  /*! @ingroup hadron
   *
   * Supports creation of (abstract) dilution operators used in 
   * stochastic sources and solutions
   */
  template<typename T>
  class DilutionOperator
  {
  private:
    //! Private type for iterators over this kind of object
    /*! NOTE: THIS NEEDS WORK */
    typedef std::iterator<std::input_iterator_tag, int> T0;
    
    //! Private type for iterators over this kind of object
    /*! NOTE: THIS NEEDS WORK */
    typedef std::iterator<std::forward_iterator_tag, int> T1;
    
  public:
    //! Public type for iterators over this kind of object
    typedef T0  const_iterator;

    //! Public type for iterators over this kind of object
    typedef T1  iterator;

    //! Virtual destructor to help with cleanup;
    virtual ~DilutionOperator() {}

    //! Get an iterator for this dilution
    /*! This uses the t0 for the time slice where the creation operator is created */
    virtual const_iterator begin(int t0) const = 0;

    //! Get an iterator for this dilution
    virtual const_iterator end() const = 0;
 
    //! Get an iterator for this dilution
    /*! Loop to modify the source and solution vectors independent of t0 */
    virtual iterator begin() = 0;

    //! Get an iterator for this dilution
    /*! Loop to modify the source and solution vectors independent of t0 */
    virtual iterator end() = 0;
 
    //! The decay direction
    virtual int getDecayDir() const = 0;

    //! The seed identifies this quark
    virtual const Seed& getSeed() const = 0;

    //! Is this element of the dilution zero?
    virtual bool zero(const_iterator iter) const = 0;

    //! Return the original source vector
    /*! NOTE: this might be slow */
    virtual T source() const = 0;
    
    //! Return the diluted source vector
    /*! MAYBE THIS SHOULD BE A CONST REFERENCE?? PROBABLY NO */
    virtual const T& dilutedSource(const_iterator iter) const = 0;
    
    //! Return the solution vector corresponding to the diluted source
    /*! MAYBE THIS SHOULD BE A CONST REFERENCE?? POSSIBLY YES */
    virtual const T& dilutedSolution(const_iterator iter) const = 0;

    //! Return the diluted source vector
    virtual T& dilutedSource(iterator iter) const = 0;
    
    //! Return the solution vector corresponding to the diluted source
    virtual T& dilutedSolution(iterator iter) = 0;
  };

} // namespace Chroma


#endif
