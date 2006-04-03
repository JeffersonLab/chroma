// $Id: mesphas_s.cc,v 3.0 2006-04-03 04:59:00 edwards Exp $


/* This routine is specific to staggered fermions! */

/* Compute the phase factors for the Nd staggered mesons. */

/* meson_phases -- meson phase factors ( Write ) */
/* wall_meson_phases -- meson phase factors for wall mesons ( Write ) */
/* j_decay -- direction of meson propagators ( Read ) */

#include "chromabase.h"
#include "mesphas_s.h"

namespace Chroma {

  void MesPhas(multi1d<LatticeReal>& meson_phases, 
	       int j_decay)
  {
    LatticeReal one;
    LatticeReal negone;
    multi1d<LatticeInteger> x(Nd);
    multi1d<Real> sgnn(2);
    int i;
    int m;
    
   
    meson_phases.resize(Nd);


            
    /* Start off by getting the coordinates of x(0), x(1), ..., x(Nd-3) */
    i = 0;
    for(m = 0;m  < ( Nd); ++m ) {
      if ( m != j_decay ) {
	x[i] = Layout::latticeCoordinate(m);
	i = i + 1;
      }
    }
  
    /* For convenience, define -1 and 1. */
    one = 1;
    negone = -one;
    
    sgnn[0] = 1;
    sgnn[1] = -sgnn[0];
  
  
    /* phase(0) = +1 */
    meson_phases[0] = one;
    
    /* Deal with all other directions separately. Yuk. */
  
    /* NOTE: in the comments below I mean x0 is the true *non*-checkerboard */
    /*   x coordinate. Also, I assume the direction j_decay has been  */
    /*   permuted such that j_decay == Nd. */
    switch(Nd)
    {
      /*# phase(1) = (-1)**(x0) */
    case 2:

      meson_phases[2] = where((x[0] & 1) == 0, LatticeReal(1), LatticeReal(-1));

      break;


      /*# phase(1) = (-1)**(x0) + (-1)**(x1) */
      /*# phase(2) = (-1)**(x0+x1) */
    case 3:
      /*## phase(1) = (-1)**(x0) + (-1)**(x1) */
      /* even context => cb=0, x=even; cb=1, x=odd */

      meson_phases[1] = where((x[0] & 1) == 0, LatticeReal(1), LatticeReal(-1))
	+ where((x[1] & 1) == 0, LatticeReal(1), LatticeReal(-1));

      /*## phase(2) = (-1)**(x0+x1) */
      /* even context => cb=0, x=even; cb=1, x=odd */

      meson_phases[2] = where(((x[0]+x[1]) & 1) == 0, LatticeReal(1), LatticeReal(-1));

      break;


      /*# phase(1) = (-1)**(x0) + (-1)**(x1) + (-1)**(x2) */
      /*# phase(2) = (-1)**(x0+x1) + (-1)**(x0+x2) + (-1)**(x1+x2) */
      /*# phase(3) = (-1)**(x0+x1+x2) */
    case 4:
            
      /*## phase(1) = (-1)**(x0) + (-1)**(x1) + (-1)**(x2) */
      meson_phases[1] = where((x[0] & 1) == 0, LatticeReal(1), LatticeReal(-1))
	+ where((x[1] & 1) == 0, LatticeReal(1), LatticeReal(-1))
	+ where((x[2] & 1) == 0, LatticeReal(1), LatticeReal(-1));

      /*## phase(2) = (-1)**(x0+x1) + (-1)**(x0+x2) + (-1)**(x1+x2) */
      meson_phases[2] = where(((x[0]+x[1]) & 1) == 0, LatticeReal(1), LatticeReal(-1))
	+ where(((x[0]+x[2]) & 1) == 0, LatticeReal(1), LatticeReal(-1))
	+ where(((x[1]+x[2]) & 1) == 0, LatticeReal(1), LatticeReal(-1));

      /*## phase(3) = (-1)**(x0+x1+x2) */
      meson_phases[3] = where(((x[0]+x[1]+x[2]) & 1) == 0, LatticeReal(1), LatticeReal(-1));
    
      break;

    default:
      QDP_error_exit("Can only handle d=2 3 4 dimensions", Nd);
    }
  
            
  }

}  // end namespace Chroma
