// $Id: phfctr.cc,v 1.1 2003-12-16 22:18:48 edwards Exp $
/*! \file
 *  \brief  This routine is specific to Wilson fermions!
 */

#include "chromabase.h"
#include "geom/setph.h"

using namespace QDP;

//! Private copy of the phases
static bool initP;
static multi1d<LatticeInteger> phases;


//!  Initialize phase factors
/*!
 * \ingroup gauge
 *
 * Setup the "phases" used by routine PHFCTR that multiplies the gauge
 * fields so as to handle boundary conditions.
 *
 * Arguments:
 *
 *  \param boundary       boundary conditions of the lattice       (Read)
 */

void setph(const multi1d<int>& boundary)
{
  START_CODE("setph");
  
  if (boundary.size() != Nd)
    QDP_error_exit("setph: error in size of boundary");

  /* Setup the "phases" used by routine PHFCTR that multiplies the gauge */
  /* fields so as to handle boundary conditions. */
  
  phases.resize(Nd);
  initP = true;

  // phases for all the directions
  for(int m = 0; m < Nd; ++m)
  {
    /* if (coord(m) == nrow(m)-1 ) then boundary[m] else 1 */
    phases[m] = where(Layout::latticeCoordinate(m) == Layout::lattSize()[m]-1,
		      LatticeInteger(boundary[m]), LatticeInteger(1));
  }
  
  // Time to skeedattle
      
  END_CODE("setph");
}



//! Multiply the gauge fields by the phase factors (-1)^X
/*!
 * \ingroup gauge
 *
 * Arguments:
 *
 *  \param u          Gauge field               (Modify)
 */

void phfctr(multi1d<LatticeColorMatrix>& u)
{
  START_CODE("phfctr");

  if (! initP)
  {
    QDPIO::cerr << "phfctr: phases not initialized" << endl;
    QDP_abort(1);
  }

  if (u.size() != Nd)
  {
    QDPIO::cerr << "phfctr: error in size of u" << endl;
    QDP_abort(1);
  }


#if 1
  for(int mu = 0; mu < Nd; ++mu)
    u[mu] *= phases[mu];

#else
// THIS CODE BLOCK IS NOT IMPLEMENTED - IT HANDLES THE COMPLEX CASE
// IT REQUIRES SOME MORE THOUGHT

  /* Only the Schroedinger functional has complex phases */
  if (SchrFun == 0)
  {
    for(int mu = 0; mu < Nd; ++mu)
      u[mu] *= real(phases[mu]);
  }
  else
  {
    switch (opt)
    {
    case FORWARD:
      for(int mu = 0; mu < Nd; ++mu)
	u[mu] *= phases[mu];

      break;

    case BACKWARD:
      for(int mu = 0; mu < Nd; ++mu)
	u[mu] *= adj(phases[mu]);

      break;

    default:
      QDP_error_exit("invalid option", opt);
    }
  }
#endif

  END_CODE("phfctr");
}
