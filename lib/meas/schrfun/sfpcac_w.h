// -*- C++ -*-
// $Id: sfpcac_w.h,v 3.3 2007-08-24 19:23:04 edwards Exp $
/*! \file
 *  \brief Schroedinger functional application of PCAC
 *
 */

#ifndef __sfpcac_w_h__
#define __sfpcac_w_h__

#include "fermact.h"
#include "util/ft/sftmom.h"

namespace Chroma 
{
  
  //! Schroedinger functional stuff
  /*!
   * @ingroup schrfun
   *
   * Compute correlation functions between axial current or pseudescalar
   * density and boundary fields using Schroedinger BC.
   *
   * Also computed, on demand, are correlation functions between both
   * boundaries with zero, one (vector current) and two (axial current or
   * pseudoscalar density) insertions in the bulk. These currents are
   * controlled by the ZVfactP and ZAfactP boolean flags.
   *
   * Compute quark propagators by using the qprop SystemSolver.
   * The initial guess for the inverter is zero.
   *
   * The results are written to the xml file.
   *
   * For further details see the comments in the dependent subroutines.
   *
   * \param state         gauge field state ( Read )
   * \param qprop         propagator solver ( Read )
   * \param phases        object holds list of momenta and Fourier phases ( Read )
   * \param ZVfactP       flag for doing Z_V measurements ( Read )
   * \param ZAfactP       flag for doing Z_A measurements ( Read )
   * \param x0            time slices with axial current insertions ( Read ) 
   * \param y0            time slices with axial current insertions ( Read ) 
   * \param xml           xml file object ( Write )
   * \param xml_group     string used for writing xml data ( Read )
   */
  void SFpcac(Handle< SystemSolver<LatticeFermion> > qprop,
	      Handle< FermState<LatticeFermion, multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> > > state,
	      const SftMom& phases,
	      bool ZVfactP, bool ZAfactP, 
	      int x0, int y0,
	      XMLWriter& xml,
	      const string& xml_group);

}

#endif
