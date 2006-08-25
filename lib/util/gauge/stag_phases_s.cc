// $Id: stag_phases_s.cc,v 3.1 2006-08-25 23:46:37 edwards Exp $
/*! \file
 *  \brief Create phases used by staggered fermions
 */

#include "chromabase.h"
#include "util/gauge/stag_phases_s.h"

namespace Chroma 
{

  namespace StagPhases 
  { 
  
    // Define these (as opposed to previous declaration)
    // This should give these members "linkage"
    multi1d<LatticeInteger> alphaClass::phases;
    multi1d<LatticeInteger> betaClass::phases;
    bool alphaClass::initP = false;
    bool betaClass::initP = false;

    // Init function for staggered K-S phases
    void alphaClass::init()
    {
      START_CODE();

      phases.resize(Nd);

      // Auxiliary: Coordinates to use in "where" clauses
      multi1d<LatticeInteger> x(Nd);
      int mu;
    
      // Fill x with lattice coordinates
      for( mu = 0; mu < Nd; mu++) {
	x[ mu ] = Layout::latticeCoordinate(mu);
      }
                                                                                
      switch(Nd) {
      case 4:

	// MILC Conventions: eta_t(x) = 1
	//                   eta_x    = -1^{t}
//      phases[3] = LatticeInteger(1);
//      phases[0] = where( x[3] % 2 == 0, LatticeInteger(1), LatticeInteger(-1));
//
//      phases[1] = where( (x[3]+x[0]) % 2 == 0, LatticeInteger(1), LatticeInteger(-1));
//      phases[2] = where( (x[3]+x[0]+x[1] ) % 2  == 0, LatticeInteger(1), LatticeInteger(-1));
//      break;

	// CPS conventions: eta_x = 1
	//			eta_y = (-1)^x

	phases[0] = LatticeInteger(1);
	phases[1] = where( x[0] % 2 == 0, LatticeInteger(1), LatticeInteger(-1));
	phases[2] = where( (x[0]+x[1] ) % 2  == 0, LatticeInteger(1), LatticeInteger(-1));
	phases[3] = where( (x[0]+x[1]+x[2] ) % 2  == 0, LatticeInteger(1), LatticeInteger(-1));
	break;

      default:
	QDP_error_exit("Staggered phases  only supported for Nd=4 just now: Nd = %d\n", Nd);
	break;
      }

      END_CODE();
    }
  
    // Init functions for Mesonic Phases
    void betaClass::init(void) 
    {
      START_CODE();

      phases.resize(Nd);
      multi1d<LatticeInteger> x(Nd);
      int mu;
    
      // Fill x with lattice coordinates
      for( mu = 0; mu < Nd; mu++) {
	x[ mu ] = Layout::latticeCoordinate(mu);
      }
                                                                                
      switch(Nd) {
      case 4:
	phases[0] = where( ((x[1]+x[2]+x[3])%2) == 0, LatticeInteger(1), LatticeInteger(-1));
	phases[1] = where( ((x[2] + x[3])%2) == 0, LatticeInteger(1), LatticeInteger(-1));
	phases[2] = where( (x[3] % 2) == 0, LatticeInteger(1), LatticeInteger(-1) );
      
	phases[3] = LatticeInteger(1);
	break;

      default:
	QDP_error_exit("Staggered phases  only supported for Nd=4 just now: Nd = %d\n", Nd);
	break;
      }
    
      END_CODE();
    }

    const LatticeInteger& alphaClass::alpha(const int i)
    {
      if ( ! initP ) { 
	init();
	initP = true;
      }
    
      return phases[i];
    }

    const LatticeInteger& betaClass::beta(const int i)
    {
      if ( ! initP ) { 
	init();
	initP = true;
      }
    
      return phases[i];
    }
    

  }  // end namespace stagphases

}  // end namespace Chroma
