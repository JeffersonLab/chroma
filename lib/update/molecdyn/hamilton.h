// -*- C++ -*-
// $Id: hamilton.h,v 1.1 2003-04-09 01:39:12 edwards Exp $

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
 * Maybe these functions should be global and not members except
 * for dsdu??
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

  //! Compute dS/dU
  virtual multi1d<LatticeColorMatrix> dsdu(const multi1d<LatticeComplex>& u) const = 0;

  //! Virtual destructor
  virtual ~HamiltonianSystem() {}
};

#endif
