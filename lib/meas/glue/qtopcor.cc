// $Id: qtopcor.cc,v 3.0 2006-04-03 04:58:58 edwards Exp $

#error "NOT FULLY CONVERTED"

/* QTOPCOR: Calculate topological charge correlation function. The routine */
/*          only calculates the squared piece (assumes the disconnected piece */
/*          averages to 0. */

/* qtop_corr -- topological charge correlation function ( Write ) */
/* qtop_den  -- topological charge density ( Read ) */

include(types.mh)

SUBROUTINE(qtopcor, qtop_corr, qtop_den, num_points)

multi1d<Real> qtop_corr(num_points);
LatticeReal qtop_den;
int num_points;
{
  include(COMMON_DECLARATIONS)
  
  /* Local Variables */
  LatticeComplex cqtop;
  LatticeComplex tmp_1;
  LatticeReal tmp_2;
  LatticeReal qtop_corr_den;
  LatticeReal lzero;
  LatticeBoolean lbit;
  LatticeBoolean lmask;
  LatticeInteger lcoord;
  LatticeInteger lradius_sq;
  LatticeInteger ltmp_1;
  LatticeInteger ltmp_2;
  
  int itmp;
  int norm;
  Real dummy;
  Double sum;
  Real rdummy;
  int cb;
  int mu;
  int r0;
  int two;
  
  START_CODE();
  
  two = 2;
  
  /* Initialize FFT and allocate memory */
  fftinit ();
      
  /* Initialize correlation function to 0 */
  qtop_corr = 0;
  
  
  /* FFT the topological charge density */
  
    lzero = 0;
  
  cqtop[0] = cmplx(qtop_den[0],lzero);
  cqtop[1] = cmplx(qtop_den[1],lzero);
  
    
  /* NOTE: I use some decay direction > Nd-1, so I fft in all directions */
  ccfft (cqtop[0], cqtop[1], Nd, FORWARD);
  
  
  /* Construct squared piece (the disconnected piece should average to zero) */
  
  /* tmp_1 = adj(FFT[qtop_den]) * FFT(qtop_den) */
    tmp_1 = adj(cqtop) * cqtop;
  
  /* FFT back to coordinate space.  */
  /* NOTE: I use some decay direction > Nd-1, so I fft in all directions */
  ccfft (tmp_1[0], tmp_1[1], Nd, BACKWARD);
  
  qtop_corr_den = real(tmp_1);
  
        
  /* Gather the coordinates of all sites and compute the radius_sq */
      lradius_sq[0] = 0;
  
  if ( Nd > 1 )
    for(mu=1; mu<=Nd-1; ++mu)
    {
      lcoord = Layout::latticeCoordinate(mu);
      itmp = nrow[mu];
      FILL(ltmp_1, itmp);
      ltmp_1 -= lcoord; /* ltmp_1 = L(mu) - lcoord */
      itmp = itmp / two;
      lbit = lcoord > itmp;	/* lbit = lcoord > L(mu)/2 */
      copymask(lcoord, lbit, ltmp_1, REPLACE);
      lradius_sq[0] += lcoord * lcoord;
    }
  
  lradius_sq[1] = lradius_sq[0];
      
  /* x-direction. x coord is set according to the checkerboard */
        
  mu = 0;
  
  ltmp_1 = 1;
  lcoord = Layout::latticeCoordinate(mu);
  lcoord = lcoord * two;
  ltmp_2 = lcoord;
  copymask(ltmp_2, lattice_odd_context, ltmp_1, ADD);    /* ltmp_2 = true x coord */
  itmp = nrow[mu];
  FILL(ltmp_1, itmp);                          
  ltmp_1 -= ltmp_2;                  /* ltmp_1 = L(mu) - lcoord */
  itmp = nrow[mu] / two;
  lbit = ltmp_2 > itmp;                          /* lbit = lcoord > L(mu)/2 */
  copymask(ltmp_2, lbit, ltmp_1, REPLACE);
  lradius_sq[0] += ltmp_2 * ltmp_2;
  
  ltmp_1 = 1;
  ltmp_2 = lcoord;
  copymask(ltmp_2, lattice_even_context, ltmp_1, ADD);
  itmp = nrow[mu];
  FILL(ltmp_1, itmp);                          
  ltmp_1 -= ltmp_2;                  /* ltmp_1 = L(mu) - lcoord */
  itmp = nrow[mu] / two;
  lbit = ltmp_2 > itmp;                          /* lbit = lcoord > L(mu)/2 */
  copymask(ltmp_2, lbit, ltmp_1, REPLACE);
  lradius_sq[1] += ltmp_2 * ltmp_2;
  
            
  
  /* For each squared radial length, sum all the qualifying density values */
  
  for(r0=0; r0 < num_points; ++r0)
  {
    /* Find the coordinates matching the current radial length */
    lmask[0] = lradius_sq[0] == r0;
    lmask[1] = lradius_sq[1] == r0;

    /* Find the total number of sites matching this radial condition */
    norm = sum(lmask[0]);
    norm += sum(lmask[1]);

    /* If some coordinates were found, sum over them */
    if ( norm > 0 )
    {
      tmp_2 = 0;
      copymask(tmp_2, lmask, qtop_corr_den, REPLACE);

      /* hypercubic sum */
      sum = sum(tmp_2[0]);
      sum += sum(tmp_2[1]);

      rdummy = FLOAT(sum);
      dummy = FLOAT(norm);
      qtop_corr[r0] = rdummy / dummy;
    }
  }
  
  /* Clean up */
                fftclr ();
  
  END_CODE();
}
