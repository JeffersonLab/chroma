#ifndef __asqtad_fermact_s_h__
#define __asqtad_fermact_s_h__

#include "fermact.h"
#include "asqtad_linop_s.h"

using namespace QDP;

class AsqtadFermAct : public AsqtadFermTypeAction<LatticeFermion>
{
public:
  //! Partial constructor
  AsqtadFermAct() {}

  //! Full constructor
  AsqtadFermAct(const Real& _Mass)
    {create(_Mass);}

  //! Creation routine
  void create(const Real& _Mass);

  //! Produce a linear operator for this action
  const LinearOperator<LatticeFermion>* linOp(const multi1d<LatticeColorMatrix>& u_fat, const multi1d<LatticeColorMatrix>& u_triple) const;

  //! Produce a linear operator M^dag.M for this action
  const LinearOperator<LatticeFermion>* lMdagM(const multi1d<LatticeColorMatrix>& u_fat, const multi1d<LatticeColorMatrix>& u_triple) const;

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
  ~AsqtadFermAct() {}

private:
  Real Mass;
};

#endif
