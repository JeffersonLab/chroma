#ifndef __asqtad_fermact_s_h__
#define __asqtad_fermact_s_h__

#include <fermact.h>
#include "actions/ferm/fermacts/asqtad_state.h"

using namespace QDP;

class AsqtadFermAct : public EvenOddStaggeredTypeFermAct<LatticeFermion>
{
public:
  //! Partial constructor
  AsqtadFermAct() {}

  //! Full constructor
  AsqtadFermAct(const Real Mass_, const Real u0_) {

    // Copy mass an u0
    Mass = Mass_;
    u0 = u0_;

  }
  
  const AsqtadConnectStateBase<LatticeFermion>* createState(const multi1d<LatticeColorMatrix>& u_) const;
    
    

  //! Produce a linear operator for this action
  const EvenOddLinearOperator<LatticeFermion>* linOp(const ConnectState& state_) const;

  //! Produce a linear operator M^dag.M for this action
  
  const LinearOperator<LatticeFermion>* lMdagM(const ConnectState& state_) 
    const;

  //! accessors 
  const Real getQuarkMass() const { 
    return Mass;
  }

  
  Real getU0() { 
    return u0;
  }


  //! Compute dS_f/dU      DO I NEED THIS??
  //  multi1d<LatticeColorMatrix> dsdu(const multi1d<LatticeColorMatrix>& u,
  //				   const LatticeFermion& psi) const;

  //! Destructor is automatic
  ~AsqtadFermAct() {}

private:
  Real Mass;
  Real u0;
};

#endif
