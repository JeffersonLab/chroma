// -*- C++ -*-
// $Id: overlap_state.h,v 1.2 2004-01-08 11:53:08 bjoo Exp $
/*! @file
 * @brief Connection state holding eigenvectors
 *
 * Holds gauge fields and eigenvectors for overlap-ish thingies
 */

#ifndef __overlap_state_h__
#define __overlap_state_h__

#include "state.h"

using namespace QDP;

template<typename T>
class OverlapConnectStateBase : public ConnectState 
{
public: 
  //! Return the eigenvalues
  virtual const multi1d<Real>& getEigVal() const = 0;

  //! Return the eigenvectors
  virtual const multi1d<T>& getEigVec() const = 0;

  //! Return the eigenvectors
  virtual const Real& getEigValMax() const = 0;

  //! Return approximation lower bound
  virtual const Real& getApproxMin() const = 0;
  virtual const Real& getApproxMax() const = 0;

};


template<typename T>
class OverlapConnectState : public OverlapConnectStateBase<T>
{
public:

  typedef Real WordBase_t;
  
  //! Constructor with no eigenvalues
  OverlapConnectState(const multi1d<LatticeColorMatrix>& u_,  // gauge field
			const WordBase_t& approxMin_,          // epsilon
			const WordBase_t& approxMax_           // approx max
    )  {
    u = u_;
    eigValMax = 0;
    approxMin = approxMin_ ;
    approxMax = approxMax_ ;
  }

  //! Constructor with e-values and e-vectors
  OverlapConnectState(const multi1d<LatticeColorMatrix>& u_,
			const multi1d<WordBase_t>& val_, 
			const multi1d<T>& vec_,
			const WordBase_t& val_max_,
			const WordBase_t& approxMin_,
			const WordBase_t& approxMax_)  {

    // hhmm, for now make a copy
    u = u_;
    eigVal = val_;
    eigVec = vec_;
    eigValMax = val_max_;
    approxMin = approxMin_;
    approxMax = approxMax_;
  }

  OverlapConnectState(const OverlapConnectState& a) : u(a.u), eigVal(a.eigVal), eigVec(a.eigVec), eigValMax(a.eigValMax), approxMin(a.approxMin), approxMax(a.approxMax) {}

  ~OverlapConnectState() {};

  //! Return the link fields needed in constructing linear operators
  const multi1d<LatticeColorMatrix>& getLinks() const {return u;}
  

  //! Return the eigenvalues
  const multi1d<WordBase_t>& getEigVal() const {return eigVal;}

  //! Return the eigenvectors
  const multi1d<T>& getEigVec() const {return eigVec;}
  
  //! Return the max eigenvalues
  const WordBase_t& getEigValMax() const {return eigValMax;}

  const WordBase_t& getApproxMin() const { return approxMin; }
  const WordBase_t& getApproxMax() const { return approxMax; }

private:
  OverlapConnectState() {}  // hide default constructur
  void operator=(const OverlapConnectState&) {} // hide =

private:
  multi1d<LatticeColorMatrix> u;
  multi1d<WordBase_t>  eigVal;
  multi1d<T>           eigVec;
  WordBase_t           eigValMax;
  WordBase_t           approxMin;
  WordBase_t           approxMax;

};
	

#endif
