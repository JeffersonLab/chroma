// $Id: mesphas_s.cc,v 1.1 2003-08-13 03:01:28 edwards Exp $


THIS CODE IS NOT FUNCTIONAL


/* This routine is specific to staggered fermions! */

/* Compute the phase factors for the Nd staggered mesons. */

/* meson_phases -- meson phase factors ( Write ) */
/* wall_meson_phases -- meson phase factors for wall mesons ( Write ) */
/* j_decay -- direction of meson propagators ( Read ) */

#include "qdp.h"



#define NEW_HOTNESS

void MesPhas(multi1d<LatticeReal>& meson_phases, 
	     multi2d<Real>& wall_meson_phases, 
	     int j_decay)
{
  LatticeInteger x1px2;
  LatticeReal negone_x1px2;
  LatticeReal temp;
  LatticeReal one;
  LatticeReal negone;
  LatticeReal sum;
  multi1d<LatticeInteger> x(Nd);
  LatticeInteger xmod;
  LatticeBoolean lbtmp;
  multi1d<Real> sgnn(2);
  int itmp;
  int i;
  int m;
  int i0;
  int i1;
  int i2;
  int ii;
  int two;
  
  START_CODE("mesphas_s");
  
  int n_cubes = 1 << (Nd-1);
  
  meson_phases.resize(Nd);
  wall_meson_phases.resize(Nd,n_cubes);

            
  /* Start off by getting the coordinates of x(0), x(1), ..., x(Nd-3) */
  i = 0;
  for(m = 0;m  < ( Nd); ++m )
  {
    if ( m != j_decay )
    {
      x[i] = latticeCoords(m);
      i = i + 1;
    }
  }
  
  /* For convenience, define -1 and 1. */
  one = 1;
  negone = -one;
    
  sgnn[0] = 1;
  sgnn[1] = -sgnn[0];
  two = 2;
  
  /* phase(0) = +1 */
  meson_phases[0] = one;
  
  for(i = 0; i < n_cubes; ++i)
    wall_meson_phases[i] = 1;
  
  /* Deal with all other directions separately. Yuk. */
  
  /* NOTE: in the comments below I mean x0 is the true *non*-checkerboard */
  /*   x coordinate. Also, I assume the direction j_decay has been  */
  /*   permuted such that j_decay == Nd. */
  switch(Nd)
  {
    /*# phase(1) = (-1)**(x0) */
  case 2:
    /* even context => cb=0, x=even; cb=1, x=odd */
#if defined(NEW_HOTNESS)
    meson_phases[1] = where((x[0] & 1) == 0, LatticeReal(1), LatticeReal(-1));
#else
    wall_meson_phases[1][0] = sgnn[0];
    wall_meson_phases[1][1] = sgnn[1];
#endif
    break;


    /*# phase(1) = (-1)**(x0) + (-1)**(x1) */
    /*# phase(2) = (-1)**(x0+x1) */
  case 3:
#if ! defined(NEW_HOTNESS)
    /* sum = (-1)**(x1) */
    itmp = 2;
    xmod = mod(x[1], itmp);
    itmp = 0;
    lbtmp = xmod == itmp;
    copymask(sum, lbtmp, one, REPLACE);
    lbtmp = !lbtmp;
    copymask(sum, lbtmp, negone, REPLACE);
#endif
        
    /*## phase(1) = (-1)**(x0) + (-1)**(x1) */
    /* even context => cb=0, x=even; cb=1, x=odd */
#if defined(NEW_HOTNESS)
    meson_phases[1] = where((x[0] & 1) == 0, LatticeReal(1), LatticeReal(-1))
                    + where((x[1] & 1) == 0, LatticeReal(1), LatticeReal(-1));
#else
    copymask(meson_phases[1][0], lattice_even_context, one, REPLACE);
    copymask(meson_phases[1][1], lattice_even_context, negone, REPLACE);
    copymask(meson_phases[1][0], lattice_even_context, sum, ADD);
    copymask(meson_phases[1][1], lattice_even_context, sum, ADD);

    /* odd context => cb=0, x=odd; cb=1, x=even */
    copymask(meson_phases[1][0], lattice_odd_context, negone, REPLACE);
    copymask(meson_phases[1][1], lattice_odd_context, one, REPLACE);
    copymask(meson_phases[1][0], lattice_odd_context, sum, ADD);
    copymask(meson_phases[1][1], lattice_odd_context, sum, ADD);
#endif

    /*## phase(2) = (-1)**(x0+x1) */
    /* even context => cb=0, x=even; cb=1, x=odd */
#if defined(NEW_HOTNESS)
    meson_phases[2] = where(((x[0]+x[1]) & 1) == 0, LatticeReal(1), LatticeReal(-1));
#else
    copymask(meson_phases[2][0], lattice_even_context, sum, REPLACE);
    copymask(meson_phases[2][1], lattice_even_context, sum, NEGATE);

    /* odd context => cb=0, x=odd; cb=1, x=even */
    copymask(meson_phases[2][0], lattice_odd_context, sum, NEGATE);
    copymask(meson_phases[2][1], lattice_odd_context, sum, REPLACE);
#endif

    for(i = 0;i  < ( n_cubes); ++i )
    {
      i0 = mod(i, two);
      ii = i / two;
      i1 = mod(ii, two);

      wall_meson_phases[1][i] = sgnn[i0];
      wall_meson_phases[1][i] += sgnn[i1];

      ii = i0 + i1;
      ii = mod(ii, two)
	wall_meson_phases[2][i] = sgnn[ii];
    }
    break;


    /*# phase(1) = (-1)**(x0) + (-1)**(x1) + (-1)**(x2) */
    /*# phase(2) = (-1)**(x0+x1) + (-1)**(x0+x2) + (-1)**(x1+x2) */
    /*# phase(3) = (-1)**(x0+x1+x2) */
  case 4:
#if ! defined(NEW_HOTNESS)
    /* sum = (-1)**(x1) + (-1)**(x2) */
            
    /* sum <- (-1)**x(1) */
    itmp = 2;
    xmod = mod(x[1], itmp);
    itmp = 0;
    lbtmp = xmod == itmp;
    copymask(sum, lbtmp, one, REPLACE);
    lbtmp = !lbtmp;
    copymask(sum, lbtmp, negone, REPLACE);

    /* temp <- (-1)**x(1) */
    itmp = 2;
    xmod = mod(x[2], itmp);
    itmp = 0;
    lbtmp = xmod == itmp;
    copymask(temp, lbtmp, one, REPLACE);
    lbtmp = !lbtmp;
    copymask(temp, lbtmp, negone, REPLACE);
    sum += temp;

            
    /* negone_x1px2 = (-1)**(x1+x2) */
        
    /* WARNING!! watch out that the addition goes unsigned */
    x1px2 = x[1];
    x1px2 += x[2];

    itmp = 2;
    xmod = mod(x1px2, itmp);
    itmp = 0;
    lbtmp = xmod == itmp;
    copymask(negone_x1px2, lbtmp, one, REPLACE);
    lbtmp = !lbtmp;
    copymask(negone_x1px2, lbtmp, negone, REPLACE);
#endif

            
    /*## phase(1) = (-1)**(x0) + (-1)**(x1) + (-1)**(x2) */
#if defined(NEW_HOTNESS)
    meson_phases[1] = where((x[0] & 1) == 0, LatticeReal(1), LatticeReal(-1))
                    + where((x[1] & 1) == 0, LatticeReal(1), LatticeReal(-1))
                    + where((x[2] & 1) == 0, LatticeReal(1), LatticeReal(-1));
#else
    /* even context => cb=0, x=even; cb=1, x=odd */
    copymask(meson_phases[1][0], lattice_even_context, one, REPLACE);
    copymask(meson_phases[1][1], lattice_even_context, negone, REPLACE);
    copymask(meson_phases[1][0], lattice_even_context, sum, ADD);
    copymask(meson_phases[1][1], lattice_even_context, sum, ADD);

    /* odd context => cb=0, x=odd; cb=1, x=even */
    copymask(meson_phases[1][0], lattice_odd_context, negone, REPLACE);
    copymask(meson_phases[1][1], lattice_odd_context, one, REPLACE);
    copymask(meson_phases[1][0], lattice_odd_context, sum, ADD);
    copymask(meson_phases[1][1], lattice_odd_context, sum, ADD);
#endif

    /*## phase(2) = (-1)**(x0+x1) + (-1)**(x0+x2) + (-1)**(x1+x2) */
#if defined(NEW_HOTNESS)
    meson_phases[2] = where(((x[0]+x[1]) & 1) == 0, LatticeReal(1), LatticeReal(-1))
                    + where(((x[1]+x[2]) & 1) == 0, LatticeReal(1), LatticeReal(-1))
                    + where(((x[1]+x[2]) & 1) == 0, LatticeReal(1), LatticeReal(-1));
#else
    /* even context => cb=0, x=even; cb=1, x=odd */
    copymask(meson_phases[2][0], lattice_even_context, negone_x1px2, REPLACE);
    copymask(meson_phases[2][1], lattice_even_context, negone_x1px2, REPLACE);
    copymask(meson_phases[2][0], lattice_even_context, sum, ADD);
    copymask(meson_phases[2][1], lattice_even_context, sum, SUBTRACT);

    /* odd context => cb=0, x=odd; cb=1, x=even */
    copymask(meson_phases[2][0], lattice_odd_context, negone_x1px2, REPLACE);
    copymask(meson_phases[2][1], lattice_odd_context, negone_x1px2, REPLACE);
    copymask(meson_phases[2][0], lattice_odd_context, sum, SUBTRACT);
    copymask(meson_phases[2][1], lattice_odd_context, sum, ADD);
#endif

    /*## phase(3) = (-1)**(x0+x1+x2) */
#if defined(NEW_HOTNESS)
    meson_phases[3] = where(((x[0]+x[1]+x[2]) & 1) == 0, LatticeReal(1), LatticeReal(-1));
#else
    /* even context => cb=0, x=even; cb=1, x=odd */
    copymask(meson_phases[3][0], lattice_even_context, negone_x1px2, REPLACE);
    copymask(meson_phases[3][1], lattice_even_context, negone_x1px2, NEGATE);

    /* odd context => cb=0, x=odd; cb=1, x=even */
    copymask(meson_phases[3][0], lattice_odd_context, negone_x1px2, NEGATE);
    copymask(meson_phases[3][1], lattice_odd_context, negone_x1px2, REPLACE);
#endif
    
    for(i = 0;i  < ( n_cubes); ++i )
    {
      i0 = mod(i, two);
      ii = i / two;
      i1 = mod(ii, two);
      ii = ii / two;
      i2 = mod(ii, two);

      wall_meson_phases[1][i] = sgnn[i0];
      wall_meson_phases[1][i] += sgnn[i1];
      wall_meson_phases[1][i] += sgnn[i2];

      ii = i0 + i1;
      ii = mod(ii, two);
      wall_meson_phases[2][i] = sgnn[ii];

      ii = i0 + i2;
      ii = mod(ii, two);
      wall_meson_phases[2][i] += sgnn[ii];

      ii = i1 + i2;
      ii = mod(ii, two);
      wall_meson_phases[2][i] += sgnn[ii];

      ii = i0 + i1 + i2;
      ii = mod(ii, two);
      wall_meson_phases[3][i] = sgnn[ii];
    }
    break;

  default:
    QDP_error_exit("Can only handle d=2 3 4 dimensions", Nd);
  }
  
            
  END_CODE("subroutine");;
}
