#ifndef __zolotarev_state_h__
#define __zolotarev_state_h__

#include "handle.h"
#include "actions/ferm/linop/improvement_terms_s.h"

using namespace QDP;

// Basic "Connect State" for ASQTAD
template<typename T>
class AsqtadConnectStateBase : public ConnectState {
 public: 

 //! Return the eigenvalues
  virtual const multi1d<LatticeColorMatrix>& getFatLinks() const = 0;
  virtual const multi1d<LatticeColorMatrix>& getTripleLinks() const = 0;

};

template<typename T>
class AsqtadConnectStateProxy : public AsqtadConnectStateBase<T>
{
public:
  //! Initialize pointer with existing pointer

  /*! Requires that the pointer p is a return value of new */
  explicit AsqtadConnectStateProxy(const AsqtadConnectStateBase<T>* p=0) : state(p) {}

  //! Copy pointer (one more owner)
  AsqtadConnectStateProxy(const AsqtadConnectStateProxy& p) : state(p.state) {}

  //! Access the value to which the pointer refers
  const AsqtadConnectStateBase<T>& operator*() const {return state.operator*();}
  const AsqtadConnectStateBase<T>* operator->() const {return state.operator->();}

  //! Return the link fields needed in constructing linear operators
  const multi1d<LatticeColorMatrix>& getLinks() const
    {state->getLinks(); }

  const multi1d<LatticeColorMatrix>& getFatLinks() const 
    {state->getFatLinks(); }

  const multi1d<LatticeColorMatrix>& getTripleLinks() const 
    { state->getTripleLinks(); }
  
 protected:
  //! Assignment
  /*! Could easily be supported, but not sure why to do so... */
  AsqtadConnectStateProxy& operator=(const AsqtadConnectStateProxy& p) {}

 private:
  Handle<const AsqtadConnectStateBase<T> >  state;
};

//! The actual thing
template<typename T>
class AsqtadConnectState : public AsqtadConnectStateBase<T>
{
 public:

  typedef Real WordBase_t;
  
  //! Full Constructor
  AsqtadConnectState(const multi1d<LatticeColorMatrix>& u_,
		     const WordBase_t& u0,
		     const multi1d<LatticeInteger>& phases) : u(u_)
    {

      // Compute fat and triple links. Multiply in K-S Phases.
      u_fat.resize(Nd);
      u_triple.resize(Nd);
      Fat7_Links(u, u_fat, u0);
      Triple_Links(u, u_triple, u0);

      int mu;
      for(mu=0; mu < Nd; mu++) { 
	u_fat(mu) *= phases(mu);
	u_triple(mu) *= phases(mu);
      }

    }

  AsqtadConnectState(const AsqtadConnectState& a) : u(a.u), u_fat(a.u_fat), u_triple(a.u_triple) {}

  ~AsqtadConnectState() {};

  //! Return the link fields needed in constructing linear operators
  const multi1d<LatticeColorMatrix>& getLinks() const { return u; }
  const multi1d<LatticeColorMatrix>& getFatLinks() const { return u_fat; }

  const multi1d<LatticeColorMatrix>& getTripleLinks() const { return u_triple; }
 private:
  AsqtadConnectState() {}  // hide default constructur
  void operator=(const AsqtadConnectState&) {} // hide =

private:
  multi1d<LatticeColorMatrix> u;
  multi1d<LatticeColorMatrix> u_fat;
  multi1d<LatticeColorMatrix> u_triple;
};
	

#endif
