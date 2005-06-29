/* + */
/* $Id: baryon_s.cc,v 1.5 2005-06-29 12:43:38 mcneile Exp $ ($Date: 2005-06-29 12:43:38 $) */

/* This routine is specific to staggered fermions! */

/* Construct baryon propagator and Write in (pseudo) NAMELIST format */

/* quark_propagator -- quark propagator ( Read ) */
/* t_source -- cartesian coordinates of the source ( Read ) */
/* j_decay -- direction of the exponential decay ( Read ) */
/* bc_spec  -- boundary condition for spectroscopy ( Read ) */

/*        ____ */
/*        \ */
/* b(t) =  >  < b(t_source, 0) b(t + t_source, x) > */
/*        /                     */
/*        ---- */
/*          x */

/* The colour components are contracted with the totally antisymmetric */
/* 'tensor' eps(a,b,c) = antisym_tensor(a,b,c). */

/*
  Run through reorder by mcneile
  Also move output to outside this file.

*/

#include "chromabase.h"

namespace Chroma {


//! Function object used for constructing the time-slice set
  class TimeSliceFunc : public SetFunc
  {
  public:
    TimeSliceFunc(int dir): dir_decay(dir) {}

    int operator() (const multi1d<int>& coordinate) const {return coordinate[dir_decay];}
    int numSubsets() const {return Layout::lattSize()[dir_decay];}

    int dir_decay;

  private:
    TimeSliceFunc() {}  // hide default constructor
  };



  void baryon_s(LatticeStaggeredPropagator & quark_propagator_in, 
		multi1d<Complex> & barprop,
		multi1d<int> & t_source,
		int j_decay, int bc_spec)
  { 
    LatticeColorMatrix  quark_propagator;

    quark_propagator =
      peekSpin(quark_propagator_in,0,0);


    LatticeComplex b_prop;
    LatticeComplex uu_quark;

    LatticeComplex quark_cpt_A ; 
    LatticeComplex quark_cpt_B ;

    Complex cdummy;
    int t0;
    int t;
    int t_eff;
    int c;
    int cc;
    int tmp;

    /* Loop counters */
    int ci_1;
    int ci_2;
    int ci_3;
    int cf_1;
    int cf_2;
    int cf_3;



    /*# 3 dimensional totally antisymmetric tensor */

    int antisym_tensor[3][3][3] ; 
    for(int ii=0 ; ii < 3 ; ++ii)
      for(int jj=0 ; jj < 3 ; ++jj)
	for(int kk=0 ; kk < 3 ; ++kk)
	  antisym_tensor[ii][jj][kk] = 0 ; 

    antisym_tensor[2][1][0] =  1;
    antisym_tensor[1][2][0] = -1;
    antisym_tensor[2][0][1] = -1;
    antisym_tensor[0][2][1] =  1;
    antisym_tensor[1][0][2] =  1;
    antisym_tensor[0][1][2] = -1;




    // --- pasted from mesons_w ---
    // Create the time-slice set
    UnorderedSet timeslice;
    timeslice.make(TimeSliceFunc(j_decay));

    // Length of lattice in j_decay direction
    int length = timeslice.numSubsets();


    // not sure why this had to be Dcomplex
    multi1d<DComplex> hsum(length);

    t0 = t_source[j_decay]; 	     /* Note j_decay = 0 is not permitted! */


    /* Set baryon propagator to zero */
    b_prop = 0;

    ci_3 = 2;
    ci_1 = 0;
    ci_2 = 1;
    c = 1;


    /* Sum over second sink u-quark colour */
    for(cf_3 = 0;cf_3  < ( Nc); ++cf_3 )
    {

      /* Sum over sink d-u diquark colour */
      for(cf_1 = 0;cf_1  < ( Nc); ++cf_1 )
	for(cf_2 = 0;cf_2  < ( Nc); ++cf_2 )
	{
	  tmp = antisym_tensor[cf_3][cf_2][cf_1];
	  if(tmp != 0)            /* otherwise no contribution! */
	  {
	    cc = c * tmp;


	    /* Build the u-u "diquark" [only direct term (23-->23) is needed!] */

	    quark_cpt_A = peekColor(quark_propagator,cf_2,ci_2 );
	    quark_cpt_B = peekColor(quark_propagator,cf_3,ci_3);

	    uu_quark =  quark_cpt_A * quark_cpt_B ;

	    /* Finally tie the u-u "diquark" to the d-quark (1) to form the baryon */

	    quark_cpt_A = peekColor(quark_propagator,cf_1,ci_1);
	    switch(cc)
	    {
	    case +1:
	      b_prop +=  quark_cpt_A * uu_quark;
	      break;
	    case -1:
	      b_prop -=  quark_cpt_A * uu_quark;
	      break;
	    }


	  }    /* end if cc (sink diquark colour) */
	}      /* end sum sink diquark colour (cf_1, cf_2) */
    }        /* end sum second sink u-quark colour (cf_3) */

    /* Project on zero momentum: Do a slice-wise sum. */
    hsum = sumMulti(b_prop, timeslice);

    for(t = 0;t  < ( length); ++t )
    {
      t_eff = (t - t0 + length) % length;

      if ( bc_spec < 0 && (t_eff+t0) >= length)
      {
	barprop[t_eff] = -hsum[t] ; 
      }
      else
	barprop[t_eff] = hsum[t];
    }



  }

  //
  //  Version needed for smearing
  //
  //
  //


  void baryon_s(
		LatticeStaggeredPropagator & quark_propagator_in_a, 
		LatticeStaggeredPropagator & quark_propagator_in_b, 
		LatticeStaggeredPropagator & quark_propagator_in_c, 
		multi1d<Complex> & barprop,
		multi1d<int> & t_source,
		int j_decay, int bc_spec)
  { 
    LatticeColorMatrix  quark_propagator_a;
    LatticeColorMatrix  quark_propagator_b;
    LatticeColorMatrix  quark_propagator_c;

    quark_propagator_a =
      peekSpin(quark_propagator_in_a,0,0);
    quark_propagator_b =
      peekSpin(quark_propagator_in_b,0,0);
    quark_propagator_c =
      peekSpin(quark_propagator_in_c,0,0);


    LatticeComplex b_prop;
    LatticeComplex uu_quark;

    LatticeComplex quark_cpt_A ; 
    LatticeComplex quark_cpt_B ;

    Complex cdummy;
    int t0;
    int t;
    int t_eff;
    int c;
    int cc;
    int tmp;

    /* Loop counters */
    int ci_1;
    int ci_2;
    int ci_3;
    int cf_1;
    int cf_2;
    int cf_3;



    /*# 3 dimensional totally antisymmetric tensor */

    int antisym_tensor[3][3][3] ; 
    for(int ii=0 ; ii < 3 ; ++ii)
      for(int jj=0 ; jj < 3 ; ++jj)
	for(int kk=0 ; kk < 3 ; ++kk)
	  antisym_tensor[ii][jj][kk] = 0 ; 

    antisym_tensor[2][1][0] =  1;
    antisym_tensor[1][2][0] = -1;
    antisym_tensor[2][0][1] = -1;
    antisym_tensor[0][2][1] =  1;
    antisym_tensor[1][0][2] =  1;
    antisym_tensor[0][1][2] = -1;




    // --- pasted from mesons_w ---
    // Create the time-slice set
    UnorderedSet timeslice;
    timeslice.make(TimeSliceFunc(j_decay));

    // Length of lattice in j_decay direction
    int length = timeslice.numSubsets();


    // not sure why this had to be Dcomplex
    multi1d<DComplex> hsum(length);

    t0 = t_source[j_decay]; 	     /* Note j_decay = 0 is not permitted! */


    /* Set baryon propagator to zero */
    b_prop = 0;

    ci_3 = 2;
    ci_1 = 0;
    ci_2 = 1;
    c = 1;


    /* Sum over second sink u-quark colour */
    for(cf_3 = 0;cf_3  < ( Nc); ++cf_3 )
    {

      /* Sum over sink d-u diquark colour */
      for(cf_1 = 0;cf_1  < ( Nc); ++cf_1 )
	for(cf_2 = 0;cf_2  < ( Nc); ++cf_2 )
	{
	  tmp = antisym_tensor[cf_3][cf_2][cf_1];
	  if(tmp != 0)            /* otherwise no contribution! */
	  {
	    cc = c * tmp;


	    /* Build the u-u "diquark" [only direct term (23-->23) is needed!] */

	    quark_cpt_A = peekColor(quark_propagator_a,cf_2,ci_2 );
	    quark_cpt_B = peekColor(quark_propagator_b,cf_3,ci_3);

	    uu_quark =  quark_cpt_A * quark_cpt_B ;

	    /* Finally tie the u-u "diquark" to the d-quark (1) to form the baryon */

	    quark_cpt_A = peekColor(quark_propagator_c,cf_1,ci_1);
	    switch(cc)
	    {
	    case +1:
	      b_prop +=  quark_cpt_A * uu_quark;
	      break;
	    case -1:
	      b_prop -=  quark_cpt_A * uu_quark;
	      break;
	    }


	  }    /* end if cc (sink diquark colour) */
	}      /* end sum sink diquark colour (cf_1, cf_2) */
    }        /* end sum second sink u-quark colour (cf_3) */

    /* Project on zero momentum: Do a slice-wise sum. */
    hsum = sumMulti(b_prop, timeslice);

    for(t = 0;t  < ( length); ++t )
    {
      t_eff = (t - t0 + length) % length;

      if ( bc_spec < 0 && (t_eff+t0) >= length)
      {
	barprop[t_eff] = -hsum[t] ; 
      }
      else
	barprop[t_eff] = hsum[t];
    }



  }



  //
  //  Some utility routines for the baryon code
  //  These may need to be moved someplace else soon.
  //







}  // end namespace Chroma
