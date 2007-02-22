// -*- C++ -*-
// $Id: wils_hamsys_w.h,v 3.1 2007-02-22 21:11:49 bjoo Exp $
/*! \file
 *  \brief HMC
 */

#ifndef __hybrid_h__
#define __hybrid_h__

#include "io/param_io.h"
#include "gaugeact.h"
#include "fermact.h"
#include "update/molecdyn/hamilton.h"
#include "actions/ferm/fermacts/unprec_wilson_fermact_w.h"

//! A simple HamSys
/*! \ingroup molecdyn
 *
 *  A barebones Hamiltonian System
 *
 *  New better names...
 *
 *  Implements 2 flavors of unprec. Wilson fermions on a Wilson gauge actiob
 */
class WilsonGaugeTwoFlavorWilsonFermHS : public SingleExactHamSys<WilsonGaugeAct, LatticeFermion, UnprecWilsonFermAct >
{
public:

  // Constructor
  WilsonGaugeTwoFlavorWilsonFermHS(const WilsonGaugeAct& gaugeAct_,
				   const UnprecWilsonFermAct& fermAct_, 
				   const InvertParam_t& invParam_) : gaugeAct(gaugeAct_), fermAct(fermAct_), invParam(invParam_) {
    p.resize(Nd);
    u.resize(Nd);
    pf.resize(1);
    psi.resize(1);
  }

  // Destructor
  ~WilsonGaugeTwoFlavorWilsonFermHS() {}

  // Get Mom
  // =======

  // Non mutable
  const  multi1d<LatticeColorMatrix>& getMom() const { return p; }

  // Mutable
  multi1d<LatticeColorMatrix>& getMom() { return p; }

  // Get U
  // =====

  // Non Mutable
  const multi1d<LatticeColorMatrix>& getU() const { return u; }

  // Mutable
  multi1d<LatticeColorMatrix>& getU() {   return u; }

  // Get PF
  // =======
  
  // Non mutable
  const multi1d<LatticeFermion>& getPF() const { return pf; }
  // Mutable
  multi1d<LatticeFermion>& getPF() { return pf; }

  
  // Get PSI
  // =======
  
  // Non mutable
  const multi1d<LatticeFermion>& getPsi() const { return psi; }
  // Mutable
  multi1d<LatticeFermion>& getPsi() { return psi; }

  // Get Gauge Action
  const WilsonGaugeAct& getGaugeAct() const { return gaugeAct; }
  
  // Get Fermion Action
  const UnprecWilsonFermAct& getFermAct() const { return fermAct; }

  // Get number of pseudofermions
  int numPF() { return 1; }

  Real getRsdCG() const { return invParam.RsdCG; }
  int getMaxCG() const { return invParam.MaxCG; }


  // Measure Kinetic Energy on a definite state (ie not mid MD)
  // per d.o.f 

  // KE = || P^2 || / (V * Nd * (Nc*Nc-1));
  Double mesKE(void) {
 
    Double KE = 0;
    for(int mu=0; mu < Nd; mu++) { 
      KE += norm2(p[mu]);
    }

    KE /= Double(Layout::vol()*Nd); // Per DOF
    // NB this will have to be scaled up for total energy contributions

    return KE;
  }

  // Measure Fermionic Energy on a definite state (ie not mid-MD)
  // per d.o.f 

  Double mesFE(void) {

    // Get a handle on a state with fermionic BC's
    Handle<const ConnectState> g_state(fermAct.createState(u));

    const Subset& s = (fermAct.linOp(g_state))->subset();
    Double FE=0;

    // FE = phi^{dag} ( M^dag M )^{-1} phi
    // where phi = M^[dagger} eta
    //       
    // this is equal to <pf | psi>

    FE =innerProductReal(pf[0], psi[0], s);
    FE /= Double(Layout::vol()*Nc*Ns);
    return FE;
  }

  Double mesPE(void) { 
    Double w_plaq, s_plaq, t_plaq, link;
    
    // Get Gauge field U with gaugeBC's
    Handle<const ConnectState> state(getGaugeAct().createState(getU()));

    MesPlq(state->getLinks(), w_plaq, s_plaq, t_plaq, link);
    Double PE=Double(Nd-1)/Double(2*(Nc*Nc-1));
    PE *= gaugeAct.getBeta();
    PE *= w_plaq;
    return PE;

  }

private:
  WilsonGaugeAct gaugeAct;
  UnprecWilsonFermAct fermAct;
  multi1d<LatticeColorMatrix> p;
  multi1d<LatticeColorMatrix> u;
  multi1d<LatticeFermion> pf;
  multi1d<LatticeFermion> psi;
  InvertParam_t invParam;
};


#endif
