// -*- C++ -*-
// $Id: hamilton.h,v 1.2 2003-12-29 20:06:54 edwards Exp $

/*! \file
 *  \brief Hamiltonian systems
 */

#ifndef __hamilton_h__
#define __hamilton_h__

#include "fermact.h"

using namespace QDP;

//! Hamiltonian system
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

class HamiltonianSystem
{
public:
  // Need some constructor for a GaugeAction and a FermionAction
  //! Complete a leap-frog trajectory
  virtual void hybTrj() = 0;

  //! Complete a campostrini-style trajectory
  virtual void campTrj() = 0;

  //! Leaps P forward step eps with optional monitor of DelH
  /*! 
   * At some point, this will need to create a linear operator 
   * and call an inverter.
   * Will hand a FermionAction in here which serves as a foundry.
   */
  virtual void leapPMX() = 0;

  //! Leap a step in the momenta
  virtual void leapP() = 0;

  //! Leap a step in the gauge fields
  virtual void leapU() = 0;

  //! Virtual destructor
  virtual ~HamiltonianSystem() {}
};

#endif
