// -*- C++ -*-
// $Id: ev_state.h,v 1.1 2003-12-02 21:30:53 edwards Exp $
/*! @file
 * @brief Connection state holding eigenvectors
 *
 * Holds gauge fields and eigenvectors for overlap-ish thingies
 */

#ifndef __overlap_state_h__
#define __overlap_state_h__

#include "handle.h"

using namespace QDP;

//! Connection-State also containing eigenvectors
/*! @ingroup fermact
 *
 * Holds a gauge field and eigenvectors for some auxilliary fermion action
 *
 * The template type is that of the eigenvectors
 */
template<typename T>
class EVConnectState : public ConnectState
{
public:
  /*! 
   * Used here is a slick way of getting the underlying word type and 
   * corresponding QDP type for some lattice object
   */
//  typedef SimpleScalar< WordType<T> >  WordBase_t;
  typedef Real  WordBase_t;    // for now, use the simpler version

  //! Partial constructor
  EVConnectState(const multi1d<LatticeColorMatrix>& u_)
    {
      u = u_;   // hhmm, for now make a copy
    }

  //! Full constructor
  EVConnectState(const multi1d<LatticeColorMatrix>& u_,  
		 const multi1d<WordBase_t>& val_, 
		 const multi1d<T>& vec_)
    {
      // hhmm, for now make a copy
      u = u_;
      eigVal = val_;
      eigVec = vec_;
    }

  //! Copy constructor
  /*! Written as a deep copy, but not expected to be used */
  EVConnectState(const EVConnectState& a) : u(a.u), eigVal(a.eigVal_), eigVec(a.eigVec_) 
    {}

  //! Destructor
  ~EVConnectState() {}

  //! Return the link fields needed in constructing linear operators
  const multi1d<LatticeColorMatrix>& getLinks() const {return u;}

  //! Return the eigenvalues
  const multi1d<WordBase_t>& getEigVal() const {return eigVal;}

  //! Return the eigenvectors
  const multi1d<T>& getEigVec() const {return eigVec;}

private:
  EVConnectState() {}  // hide default constructur
  void operator=(const EVConnectState&) {} // hide =

private:
  multi1d<LatticeColorMatrix> u;
  multi1d<WordBase_t>  eigVal;
  multi1d<T>           eigVec;
};


#endif
