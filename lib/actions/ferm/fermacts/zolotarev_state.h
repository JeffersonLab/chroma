#ifndef __zolotarev_state_h__
#define __zolotarev_state_h__

#include "handle.h"



using namespace QDP;

template<typename T>
class ZolotarevConnectStateBase : public ConnectState {
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
class ZolotarevConnectStateProxy : public ZolotarevConnectStateBase<T>
{
public:
  //! Initialize pointer with existing pointer

  /*! Requires that the pointer p is a return value of new */
  explicit ZolotarevConnectStateProxy(const ZolotarevConnectStateBase<T>* p=0) : state(p) {}

  //! Copy pointer (one more owner)
  ZolotarevConnectStateProxy(const ZolotarevConnectStateProxy& p) : state(p.state) {}

  //! Access the value to which the pointer refers
  const ZolotarevConnectStateBase<T>& operator*() const {return state.operator*();}
  const ZolotarevConnectStateBase<T>* operator->() const {return state.operator->();}

  //! Return the link fields needed in constructing linear operators
  const multi1d<LatticeColorMatrix>& getLinks() const 
    {state->getLinks();}

  //! Return the eigenvalues
  const multi1d<Real>& getEigVal() const
    {state->getEigVal();}

  //! Return the eigenvectors
  const multi1d<T>& getEigVec() const
    {state->getEigVec();}

  //! Return the eigenvectors
  const Real& getEigValMax() const
    {state->getEigValMax();}

  //! Return the epsilon
  const Real& getApproxMin() const
    { state->getApproxMin(); }

  //! Return the upper end of the approximation
  const Real& getApproxMax() const
    { state->getApproxMax(); }

protected:
  //! Assignment
  /*! Could easily be supported, but not sure why to do so... */
  ZolotarevConnectStateProxy& operator=(const ZolotarevConnectStateProxy& p) {}

private:
  Handle<const ZolotarevConnectStateBase<T> >  state;
};

template<typename T>
class ZolotarevConnectState : public ZolotarevConnectStateBase<T>
{
 public:

  typedef Real WordBase_t;
  
  //! Constructor with "Default" approx max
  ZolotarevConnectState(const multi1d<LatticeColorMatrix>& u_,
			const WordBase_t& approxMin_)
    {
      u = u_;
      eigValMax = 0;
      approxMin = approxMin_;
      approxMax = Real(2)*Real(Nd);
    }


  //! Constructor with no eigenvalues
  ZolotarevConnectState(const multi1d<LatticeColorMatrix>& u_,  // gauge field
			const WordBase_t& approxMin_,          // epsilon
			const WordBase_t& approxMax_           // approx max
			)
    {
      u = u_;
      eigValMax = 0;
      approxMin = approxMin_ ;
      approxMax = approxMax_ ;
    }

  //! Constructor with e-values and e-vectors
  ZolotarevConnectState(const multi1d<LatticeColorMatrix>& u_,
			const multi1d<WordBase_t>& val_, 
			const multi1d<T>& vec_,
			const WordBase_t& val_max_)
    {

      // hhmm, for now make a copy
      u = u_;
      eigVal = val_;
      eigVec = vec_;
      eigValMax = val_max_;

      int n_val = val_.size();
      if( n_val > 0 ) { 
	approxMin = eigVal[ n_val - 1];
	approxMax = eigValMax;
      }
      else {
	//
	QDP_error_exit("ZolotarevConnectState: E-vector based constructor called with no e-vectors");
      }
    }

  ZolotarevConnectState(const ZolotarevConnectState& a) : u(a.u), eigVal(a.eigVal), eigVec(a.eigVec), eigValMax(a.eigValMax), approxMin(a.approxMin), approxMax(a.approxMax) {}

  ~ZolotarevConnectState() {};

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
  ZolotarevConnectState() {}  // hide default constructur
  void operator=(const ZolotarevConnectState&) {} // hide =

private:
  multi1d<LatticeColorMatrix> u;
  multi1d<WordBase_t>  eigVal;
  multi1d<T>           eigVec;
  WordBase_t           eigValMax;
  WordBase_t           approxMin;
  WordBase_t           approxMax;

};
	

#endif
