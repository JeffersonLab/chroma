#ifndef __asqtad_fermact_s_h__
#define __asqtad_fermact_s_h__

#include <fermact.h>
#include "actions/ferm/fermacts/asqtad_state.h"

using namespace QDP;

class EvenOddPrecAsqtadFermAct : public EvenOddPrecAsqtadFermTypeAction<LatticeFermion>
{
public:
  //! Partial constructor
  EvenOddPrecAsqtadFermAct() {}

  //! Full constructor
  EvenOddPrecAsqtadFermAct(const Real Mass_, const Real u0_) {

    // Copy mass an u0
    Mass = Mass_;
    u0 = u0_;


    // Compute KS phases
    phases.resize(Nd);

    // Auxiliary: Coordinates to use in "where" clauses
    multi1d<LatticeInteger> x(Nd);
      
    int mu;
    
    // Fill x with lattice coordinates
    for( mu = 0; mu < Nd; mu++) { 
      x[ mu ] = Layout::latticeCoordinate(mu);
    }
  
    // Compute the actual phases. 
    phases[0] = LatticeInteger(1);
    phases[1] = where((x[0] % 2) == 0, LatticeInteger(1), LatticeInteger(-1));
    phases[2] = where( ((x[0]+x[1])%2) == 0, LatticeInteger(1), LatticeInteger(-1));
    phases[3] = where( ((x[0]+x[1]+x[2])%2) == 0, LatticeInteger(1), LatticeInteger(-1));
    
  }
  

  //! Produce a linear operator for this action
  const EvenOddPrecLinearOperator<LatticeFermion>* linOp(const ConnectState& state_) const;

  //! Produce a linear operator M^dag.M for this action
  
  const EvenOddPrecLinearOperator<LatticeFermion>* lMdagM(const ConnectState& state_) const {} ;
  //! accessors 
  Real getQuarkMass() { 
    return Mass;
  }

  
  Real getU0() { 
    return u0;
  }

  //! Create Connect State from u
  const AsqtadConnectStateBase<LatticeFermion>* createState(const multi1d<LatticeColorMatrix>& u) {
    return new AsqtadConnectState<LatticeFermion>(u, u0, phases);
  }

#if 0
// THIS SHOULD NOT BE NEEDED
  //! Compute quark propagator
  void qprop(LatticeFermion& psi, 
	     const multi1d<LatticeColorMatrix>& u, 
	     const LatticeFermion& chi, 
	     enum InvType invType,
	     const Real& RsdCG, 
	     int MaxCG, int& ncg_had) const;
#endif

  //! Compute dS_f/dU      DO I NEED THIS??
//  multi1d<LatticeColorMatrix> dsdu(const multi1d<LatticeColorMatrix>& u,
//				   const LatticeFermion& psi) const;

  //! Destructor is automatic
  ~EvenOddPrecAsqtadFermAct() {}

private:
  Real Mass;
  Real u0;
  multi1d<LatticeInteger> phases;
};

#endif
