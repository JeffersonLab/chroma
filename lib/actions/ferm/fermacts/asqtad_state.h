#ifndef __asqtad_state_h__
#define __asqtad_state_h__


#include "handle.h"
#include "state.h"
#include "actions/ferm/linop/improvement_terms_s.h"


namespace Chroma 
{ 
//! Basic "Connect State" for ASQTAD
/*! 
 * \ingroup fermact
 */
template<typename T>
class AsqtadConnectStateBase : public ConnectState {
 public: 

  virtual const multi1d<LatticeColorMatrix>& getFatLinks() const = 0;
  virtual const multi1d<LatticeColorMatrix>& getTripleLinks() const = 0;

};


//! The actual Asqtad thing
/*! 
 * \ingroup fermact
 */
template<typename T>
class AsqtadConnectState : public AsqtadConnectStateBase<T>
{
 public:

  typedef Real WordBase_t;
  
  //! Full Constructor
  // Make deep copies here
  AsqtadConnectState(const multi1d<LatticeColorMatrix>& u_,
		     const multi1d<LatticeColorMatrix>& u_fat_,
		     const multi1d<LatticeColorMatrix>& u_triple_)
    : u(u_), u_fat(u_fat_), u_triple(u_triple_)  { };


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
	
}; // End Namespace Chroma


#endif
