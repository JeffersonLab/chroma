// -*- C++ -*-
// $Id: hamilton.h,v 1.3 2003-12-30 19:50:10 edwards Exp $
/*! \file
 *  \brief Hamiltonian systems
 */

#ifndef __hamilton_h__
#define __hamilton_h__

#include "fermact.h"

using namespace QDP;

//! Abstract Hamiltonian system
/*!
 * \ingroup molecdyn
 *
 * This is a really crude first attempt at a Hamiltonian system.
 * Lots more thought needed here.
 *
 * Maybe these functions should be global? No, I think there is a case
 * to be made for having an "MD Integrator" class with funcs.
 * An HMDIntegrator should take a HamSys template param'd by action types.
 * I suspect I should have a typelist of actions, (e.g., can have
 * more than 1 ferm action). The HamSys constructor would then take
 * instances of those actions. 
 *
 * There should be some general notion of "Coordinates" and "Conjugate Momenta".
 * How pseudo-ferm cleanly fit in that picture in this code is not clear
 * (to me, at least).
 *
 * NOTE: dsdu is found in the fermion (or gauge) action
 */

template<class GA, class FA> 
class HamSys
{
public:
  //! Compute dS/dU
  virtual void dsdu(multi1d<LatticeColorMatrix>& ds_u) const = 0;

  //! Virtual destructor
  virtual ~HamSys() {}
};


//! Abstract molecular dynamics integrator
/*!
 * \ingroup molecdyn
 *
 * This is a really crude first attempt at a molecular dynamics integrator
 * Lots more thought needed here.
 *
 * An HMDIntegrator should take a HamSys template param'd by action types.
 * I suspect I should have a typelist of actions, (e.g., can have
 * more than 1 ferm action). The HamSys constructor would then take
 * instances of those actions. 
 *
 * There should be some general notion of "Coordinates" and "Conjugate Momenta".
 * How pseudo-ferm cleanly fit in that picture in this code is not clear
 * (to me, at least).
 */

template<class GA, class FA, template<class,int> class HS>
class HybInt
{
public:
  //! Do an integration
  virtual void operator()(HS<GA,FA>& ham) = 0;

  //! Leaps P forward step eps with optional monitor of DelH
  virtual void leapPMX(HS<GA,FA>& ham) = 0;

  //! Leap a step in the momenta
  virtual void leapP(HS<GA,FA>& ham) = 0;

  //! Leap a step in the gauge fields
  virtual void leapU(HS<GA,FA>& ham) = 0;

  //! Virtual destructor
  virtual ~HybInt() {}
};



//! A leap-frog integrator
/*!
 * \ingroup molecdyn
 *
 */
template<class GA, class FA, template<class,int> class HS>
class HybTrj : public HybInt<GA,FA,HS>
{
public:
  //! Constructor
  HybTrj(const Real& dt) {}

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
  Real dt;  // probably want some struct here of params
};


//! Complete a campostrini-style trajectory
template<class GA, class FA, template<class,int> class HS>
class CampTrj : public HybInt<GA,FA,HS>;


#endif
