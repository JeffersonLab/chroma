// $Id: fuzwilp.cc,v 3.0 2006-04-03 04:58:58 edwards Exp $ fuzwilp.h,v 1.1 2004/04/26 16:12:49 mcneile Exp $

// // version with added tmax (ACI)
/*! \file
 *  \brief Calculate ape-fuzzed Wilson loops
 */

#include "chromabase.h"
#include "meas/smear/ape_smear.h"
#include "meas/glue/fuzwilp.h"

namespace Chroma { 
//! Calculate ape-fuzzed Wilson loops
/*!
 * \ingroup glue
 *
 * Computes time-like APE_fuzzed Wilson loops, including non-planar loops,
 *
 * This version makes APE-smeared links with no blocking as required
 *          for potential clculations
 *
 * Warning: This version is VERY Slow as it has non-recursive shifting
 *          of some link products
 * Warning: this works only for Nc = 2 and 3 ! (Projection of
 *                                              smeared/blocked links)
 *
 * Warning: this version assumes the space-like directions (perpendicular
 *          to j_decay) to have equal length.
 *
 * \param u         gauge field ( Read )
 * \param j_decay   'time' direction for 'fuzzed' Wilson loops ( Read )
  * \param tmax      maximum time-extent loops ( Read )
 * \param n_smear   number of applying smearing to the gauge links ( Read )
 * \param sm_fact   "smearing" factor = weight of old link w. r. to staples ( Read )
 * \param BlkAccu   accuracy in fuzzy link projection ( Read )
 * \param BlkMax    maximum number of iterations in fuzzy link projection ( Read ) 
 */

void fuzwilp( const multi1d<LatticeColorMatrix>& u, 
        int j_decay, int tmax, int n_smear,
	const Real& sm_fact, const Real& BlkAccu, int BlkMax, 
	XMLWriter& xml, const string& xml_group)

{ 
  START_CODE();

    
  int lengthr; 
  int lengtht; 
  int lengthrs;

  multi1d<int> nrow(Nd);
  nrow = Layout::lattSize();


  lengtht  = nrow[j_decay] / 2;
  lengthr  = nrow[0]/2;
  lengthrs = (lengthr) * (lengthr+1) / 2;

  if (tmax < lengtht) 
    {
      lengtht = tmax;
      QDPIO::cout << " using lengtht = " << lengtht << endl;
      
    }
  multi1d<LatticeColorMatrix> u_smear(Nd);
  multi1d<LatticeColorMatrix> u_tmp(Nd);
  multi2d<LatticeColorMatrix> u_prod( (Nd-1), lengthr);
  LatticeColorMatrix   tmp_tog;
  LatticeColorMatrix   up_t;
  LatticeColorMatrix   u_corn;
  LatticeColorMatrix   tmp_1;
  LatticeColorMatrix   tmp_2;
  LatticeReal wl_trace;

  multi2d<Double> fuz_wlp1(lengthr, lengtht);
  multi2d<Double> fuz_wlp2(lengthrs, lengtht);
  Double ddummy;
  Real rdummy;
  Real vol;
  int bl_level;
  int mu;
  int nu;
  int mum;
  int nun;
  int t;
  int r;
  int s;
  int i;
  int n;

/* Check that all direction perpendicular to j_decay are equal. */
  for(mu = 1;mu  < ( Nd); ++mu )
    if( mu != j_decay && nrow[mu] != nrow[0] )
    {
      QDPIO::cout << " j_decay = " << j_decay << "and lengthr =" << nrow[0] << endl;
      
      QDP_error_exit("Wrong lattice size for Wilson loops: ", mu, nrow[mu]);
    }

  /* First construct the smeared links for mu != j_decay */
  
  /* Copy u's to u_smear */
  u_smear = u;

  bl_level = 0;
  
  /* Smear the space-like links n_smear times */
  for(i = 1;i  <= ( n_smear); ++i )
  {
    for(mu = 0;mu  < ( Nd); ++mu )
    {
      if( mu != j_decay )
	  APE_Smear (u_smear, u_tmp[mu], mu, bl_level, sm_fact, BlkAccu, BlkMax, j_decay);
    }
    u_smear = u_tmp;
  }

  // Compute the average link
  Double link = 0.0;
  
  for(int mu=0; mu < Nd; ++mu)
    link += sum(real(trace(u_smear[mu])));

  link /= double(Layout::vol()*Nd*Nc);
  QDPIO::cout << "Average link after smearing: " << link << endl;

  
/* Construct products of smeared links (for mu != j_decay) */

//     u_prod[mum][r]    r=0,1,2,....
//    mum =mu direction with mu.ne.j_decay, i.e. 0,1,2
//       ----> mum
//     Us  Us  Us
//     --- --- --- ....
//
  mum = -1;
  for(mu = 0;mu  < ( Nd); ++mu )
    if( mu != j_decay )
    {
      mum = mum + 1;
      u_prod[mum][0] = u_smear[mu];
      for(r=1;r < lengthr; ++r)
      {
       tmp_tog = shift(u_prod[mum][r-1],FORWARD, mu);
       tmp_1   = u_smear[mu] * tmp_tog ; 
//
//    Up[r]   = Us.[Us.Us.Us....Us]
//                  [Up[r-1]shifted +1 in mu
//
       u_prod[mum][r] = tmp_1;
      }
    }
    
  fuz_wlp1 = 0;			/* initialize the planar  Wilson loops */
  fuz_wlp2 = 0;			/* initialize the nplanar Wilson loops */

  /* Compute 'fuzzied' Wilson loops (with normal links in 'time' direction!) */

//
//   For t x r planar case will use:
//		U3	
//	--------->------
//	|		|
//  	|		|		  	^
//  U4	^		^  U2			t |
//	|		| 	        (j_decay) |---->
//	|-------->------|		     	    r (mu)
//		U1
// so that
//	  W = Tr[ U1 * U2 * U3+ * U4+ ]
//

  for(t = 0;t  < ( lengtht); ++t )
  {

    /* Compute product of un-fuzzed links in t (j_decay) direction */
   QDPIO::cout << "t= " << t << endl;
    if( t == 0 )
    {
	up_t = u[j_decay];
    }
    else
    {
	tmp_tog = shift(up_t, FORWARD, j_decay);
	up_t    = u[j_decay] * tmp_tog; 
//
//	| U
//	| U    last Up shifted up one in t direction
//	| U
//	| U
//      ----
//      | U
//		use new up_t for U4 
//
    }

    mum = -1;
    for(mu = 0;mu  < ( Nd); ++mu )
    {
      if( mu != j_decay )
      {
	mum = mum + 1;
//   QDPIO::cout << "mu= " << mu << endl;
	for(r = 0;r  < ( lengthr); ++r )
	{
	    /* Gather 'time-like' link segment from r-direction */

	    tmp_tog = shift(up_t, FORWARD, mu);
	    n = r;
	    while( n > 0 )
	    {
	      tmp_1 = shift(tmp_tog, FORWARD, mu);
	      tmp_tog = tmp_1;
	      n = n - 1;
	    }
	    tmp_1 = tmp_tog;   
//
//	.		|
//	.		^   tmp_1  (Up shifted r+1 in mu direction)
//	.		|
//	.		|
//				use tmp_1 for U2

	    /* Gather 'space-like' link segment from t-direction */

	    tmp_tog = shift(u_prod[mum][r], FORWARD, j_decay);
	    n = t;
	    while( n > 0 )
	    {
	      tmp_2 = shift( tmp_tog, FORWARD, j_decay);
	      tmp_tog = tmp_2;
	      n = n - 1;
	    }
//
//
//	  tmp_tog (u_prod shifted t+1 in t direction)
//	-------->--------
//
//
//
//	.................
//				use tmp_tog for U3
//
//
	    /* Now complete the planar Wilson loop */

	    tmp_2 = tmp_1 * adj(tmp_tog);
//
//		U2 * U3+
//
	    tmp_1 = tmp_2 * adj(up_t);	
//
//		U2 * U3+ * U4+
//
 	    tmp_2 = u_prod[mum][r] * tmp_1;
//
//		U1 * U2 * U3+ * U4+
//
	    wl_trace = real(trace(tmp_2));
	    fuz_wlp1[r][t] += sum(wl_trace);

//  QDPIO::cout << "t= " << t << "   r= " << r << endl;
//  QDPIO::cout << "fuz_wlp1= "  << 2*fuz_wlp1[r][t]/ 
//		double((Layout::vol())*(Nd-1)*(Nd-2)*Nc) << endl;

	  /* Now do non-planar loops */

//
//	We are currently inside t, mu, and r loops
//	   _____________
//	  /		|
//	 /		|
//	/		|
//	|		|
//	|		|
//	^ t (j_decay)	|
//	|   ____________|
//	|  /   s (nu)
//	| /r (mu)	
//	|/ 
//
	  nun = -1;
	  for(nu = 0;nu  < ( Nd); ++nu )
	  {
	    if ( nu != j_decay && nu == mu )
	      /* advance nun, since next "if" is not satisfied! */
	      nun = nun + 1;

	    if ( nu != j_decay && nu != mu )
	    {
	      nun = nun + 1;
//   QDPIO::cout << "nu= " << nu << endl;
	      for(s = 0;s  <= ( r); ++s )
	      {
		/* Construct the 'forward' space-like segments */

		  /* Gather 'nu link segment' from r-direction */

		  tmp_tog = shift(u_prod[nun][s], FORWARD, mu);

		  n = r;
		  while( n > 0 )
		  {
		    tmp_1 = shift( tmp_tog, FORWARD, mu);
		    tmp_tog = tmp_1;
		    n = n - 1;
		  }
		  u_corn = u_prod[mum][r] * tmp_tog;
//
//	   ___________
//	  /     s (nu)
//	 / r (mu)  		1st u_corn contribution  
//	/
//

		  /* Gather 'mu link segment' from s-direction */

		  tmp_tog = shift(u_prod[mum][r], FORWARD, nu);
		  n = s;
		  while( n > 0 )
		  {
		    tmp_1 = shift(tmp_tog, FORWARD, nu);
		    tmp_tog = tmp_1;
		    n = n -1;
		  }
		  u_corn += u_prod[nun][s] * tmp_tog;
//
//	              /  
//	             / r (mu)  	2nd u_corn contribution
//	____________/
//	   s (nu)

		/* Now collect to construct the Wilson loop */

		  /* Gather 'time-like' link segment from r-direction first */

		  tmp_tog = shift(up_t, FORWARD, mu);
		  n = r;
		  while( n > 0 )
		  {
		    tmp_1 = shift(tmp_tog, FORWARD, mu);
		    tmp_tog = tmp_1;
		    n = n - 1;
		  }
		  tmp_1 = tmp_tog;

		  /* Gather 'time-like' link segment from s-direction next */

		  tmp_tog = shift(tmp_1, FORWARD, nu);
		  n = s;
		  while( n > 0 )
		  {
		    tmp_2 = shift(tmp_tog, FORWARD, nu);
		    tmp_tog = tmp_2;
		    n = n - 1;
		  }
		  tmp_1 = tmp_tog;
//
//  		.	|
//		.	|
//	.	.	|
//    	.	.	^   tmp_1 is now Up shifted by s,t in mu,nu 
//    	.	.	|
//    	.	.	|
//    	.	. . . . |
//    	.   .
//	.
//

		  /* Gather 'space-like' link segment from t-direction */

		  tmp_tog = shift(u_corn, FORWARD, j_decay);
		  n = t;
		  while( n > 0 )
		  {
		    tmp_2 = shift(tmp_tog, FORWARD, j_decay);
		    tmp_tog = tmp_2;
		    n = n - 1;
         	  }
//
//	   _____________
//	  /		.
//	 /		.
//	/		.
//			.	tmp_tog is u_corn shifted by t+1
//	.		.
//	.       	.
//	.   .............
//	.  .   s (nu)
//	. . r (mu)	
//	.. 
//
		  /* Now complete the 'forward' non-planar Wilson loop */

		  tmp_2 = tmp_1 * adj(tmp_tog);
		  tmp_1 = tmp_2 * adj(up_t);
		  tmp_2 = u_corn * tmp_1;


//			U3 (tmp_tog)
//		   _____________
//		  /		|
//		 /		|
//		/		|
//		|		| U2 (tmp_1)
//		|		|
// U4 (Up)	|		|
//		|   ____________|
//		|  /   
//		| /   U1 (u_corn)	
//		|/ 
//
//    so W = Tr [ U1 * U2 * U3+ * U4+ ]
//
		  wl_trace = real(trace(tmp_2));
		  n = r * (r+1) / 2 + s;
		  fuz_wlp2[n][t] += sum(wl_trace);

		/* Construct the 'backward' space-like segments */

		  /* Gather 'nu link segment' from r-direction */

		  tmp_tog = shift(u_prod[nun][s], FORWARD, mu);

		  n = r;
		  while( n > 0 )
 		  {
		    tmp_1 = shift(tmp_tog, FORWARD, mu);
		    tmp_tog = tmp_1;
		    n = n - 1;
		  }
		  tmp_1 = tmp_tog;

		  /* Now fetch this from backward s-direction */

		  tmp_tog = shift(tmp_1, BACKWARD, nu);
		  n = s;
		  while( n > 0 )
		  {
		    tmp_2 = shift(tmp_tog, BACKWARD, nu);
                    tmp_tog = tmp_2;
		    n = n - 1;
		  }
		  u_corn = u_prod[mum][r] * adj(tmp_tog);

		  /* Gather 'mu link segment' from s-direction */

		  tmp_tog = shift(u_prod[mum][r], BACKWARD, nu);
		  n = s;
		  while( n > 0 )
		  {
		   tmp_1 = shift(tmp_tog, BACKWARD, nu);
		   tmp_tog = tmp_1;
		    n = n - 1;
		  }
		  tmp_1 = tmp_tog;

		  /* Gather 'nu link segment' from backward s-direction */

		  tmp_tog = shift(u_prod[nun][s], BACKWARD, nu);
		  n = s;
		  while( n > 0 )
		  {
		    tmp_2 = shift(tmp_tog, BACKWARD, nu);
		    tmp_tog = tmp_2;
		    n = n - 1;
		  }
		  u_corn += adj(tmp_tog) * tmp_1;


		/* Now collect to construct the Wilson loop */

		  /* Gather 'time-like' link segment from r-direction first */

		  tmp_tog = shift(up_t, FORWARD, mu);
		  n = r;
		  while( n > 0 )
		  {
		    tmp_1 = shift(tmp_tog, FORWARD, mu);
		    tmp_tog = tmp_1;
		    n = n - 1;
		  }
		  tmp_1 = tmp_tog;

		  /* Gather 'time-like' link segment from backward s-direction next */

		  tmp_tog = shift(tmp_1, BACKWARD, nu);
		  n = s;
		  while( n > 0 )
		  {
		    tmp_2 = shift(tmp_tog, BACKWARD, nu);
		    tmp_tog = tmp_2;
		    n = n - 1;
		  }
		  tmp_1 = tmp_tog;

		  /* Gather 'space-like' link segment from t-direction */

		  tmp_tog = shift(u_corn, FORWARD, j_decay);
		  n = t;
		  while( n > 0 )
		  {
		    tmp_2 = shift(tmp_tog, FORWARD, j_decay);
		    tmp_tog = tmp_2;
		    n = n - 1 ;
		  }

		  /* Now complete the 'backward' non-planar Wilson loop */

		  tmp_2 = tmp_1 * adj(tmp_tog);
		  tmp_1 = tmp_2 * adj(up_t);
		  tmp_2 = u_corn * tmp_1;
		  wl_trace = real(trace(tmp_2));
		  n = r * (r+1) / 2 + s;
		  fuz_wlp2[n][t] += sum(wl_trace);

	      }    /* end s loop */
	    }      /* end nu != j_decay & nu != mu */
	  }        /* end nu loop */
	}          /* end r loop */
      }            /* end mu != j_decay */
    }              /* end mu loop */
  }                /* end t loop */

  ddummy = 1.0 / double (Layout::vol()*Nc*(Nd-1)) ;

  push(xml,"fuz_wlp1");			// XML tag for fuz_wlp1
  write(xml, "lengthr", lengthr);

  multi1d<Real> wloopr(lengtht);
  for(r = 0; r < lengthr; ++r)
  {  
    for(t = 0; t < lengtht; ++t)
    {
      fuz_wlp1[r][t] = fuz_wlp1[r][t] * ddummy;
      wloopr[t]      = fuz_wlp1[r][t];
    } 
    write(xml, "r", r);       
    write(xml, "wloopr", wloopr);       // write out fuz_wlp1
  } // end for r 
  pop(xml);				// XML end tag for fuz_wlp1

   ddummy = 1.0 / double(Layout::vol()*8*Nc*(Nd-1)*(Nd-2) ) ;

  for(t = 0;t  < ( lengtht); ++t )
  {
    n = -1;
    for(r = 0;r  < ( lengthr); ++r )
      for(s = 0;s  <= ( r); ++s )
      {
	n = n + 1;
	fuz_wlp2[n][t] = fuz_wlp2[n][t] * ddummy;
      }
  }
  multi1d<Real> wlooprs(lengtht);
  push(xml,"fuz_wlp2");			// XML tag for fuz_wlp2
  write(xml, "lengthr", lengthr);

  n = -1;
  for(r = 0; r < lengthr; ++r)
  {
    write(xml, "r", r);       
    for(s = 0;s  <= ( r); ++s )
      {
	n = n + 1;
    	write(xml, "s", s);       
        for(t = 0; t < lengtht; ++t)
          wlooprs[t]      = fuz_wlp2[n][t];
        write(xml, "wlooprs", wlooprs);       // write out fuz_wlp2
      }  // end for s
  } // end for r 
  pop(xml);				// XML end tag for fuz_wlp2
  QDPIO::cout << "fuz_wlp1 and fuz_wlp2 written to .xml file " << endl;
}

};
