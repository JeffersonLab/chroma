// -*- C++ -*-
// $Id: hybrid.h,v 3.0 2006-04-03 04:59:10 edwards Exp $
/*! \file
 *  \brief HMC
 */

#ifndef __hybrid_h__
#define __hybrid_h__

#include "update/molecdyn/hamilton.h"
#include "gaugeact.h"
#include "fermact.h"


//! A simple HamSys
/*! \ingroup molecdyn
 *
 *  A barebones Hamiltonian System
 *
 *  New better names...
 *
 *  Implements 2 flavors of unprec. Wilson fermions on a Wilson gauge actiob
 */
class TwoFlavorWilsonGaugeFermHS : public DirectExactHamSys<WilsonGaugeAct,UnprecWilsonFermAct>
{
public:
  //! General FermBC
  WilsonGaugeFermHS(const HybTrj& hybint_) : hybint(hybint_) {}

  //! Do an integration
  void operator()(multi1d<LatticeColorMatrix>& u) const;

  //! Destructor is automatic
  ~WilsonGaugeFermHS() {}

private:
  HybTrj  hybint;
};






//! Parameters for Leap-Frog integration
struct HybTrjParam_t
{
  Real         Nf;         // Number of Flavours
  Real         dt;         // Step size
  bool         MonitordH;  // Measure dH for every step on the trajectory
  Real         tau0;       // Average trajectory length
  AlgETrj_t    AlgETrj;    // Specify trajectory length algorithm
  int          IntrplOrd;  // Order of interpolation
  int          Algorithm;  // No Metropolis step (i.e., use Hybrid
  bool         AlgLPStp;   // Algorithm requires P and psi at end of traj.
  bool         RefNextTrj; // Refresh P or chi for next trajectory?
  bool         RalgP;      // R-algorithm mode?

  int          MaxCG;
  Real         RsdCGMD;
  Real         RsdCGMDMax;
  Real         RsdCGMDFctr;
};


//! A single trajectory leap-frog integrator
/*!
 * \ingroup molecdyn
 *
 */
template<class GA, class FA, template<class,int> class HS>
class HybTrj : public HybInt<GA,FA,HS>
{
public:
  //! Constructor
  HybTrj(const HybTrjParam_t& param_) : param(param_) {}

  //! Do an integration
  void operator()(HS<GA,FA>& ham);

  //! Leaps P forward step eps with optional monitor of DelH
  void leapPMX(HS<GA,FA>& ham);

  //! Leap a step in the momenta
  void leapP(HS<GA,FA>& ham);

  //! Leap a step in the gauge fields
  void leapU(HS<GA,FA>& ham);

  //! destructor
  ~HybInt() {}

private:
  HybTrjParam_t  param;
};


//! Complete a campostrini-style MD trajectory
template<class GA, class FA, template<class,int> class HS>
class CampTrj : public HybInt<GA,FA,HS>;

//! Complete a R-algorithm MD trajectory
template<class GA, class FA, template<class,int> class HS>
class RHybTrj : public HybInt<GA,FA,HS>;




//! HMC
/*! \ingroup molecdyn
 *
 *  A barebones HMC integrator
 */
class SimpleHMC : public HMCTrj
{
public:
  //! General FermBC
  SimpleHMC(const HybTrj& hybint_) : hybint(hybint_) {}

  //! Do an integration
  void operator()(multi1d<LatticeColorMatrix>& u) const;

  //! Destructor is automatic
  ~SimpleHMC() {}

private:
  HybTrj  hybint;
};

#endif
