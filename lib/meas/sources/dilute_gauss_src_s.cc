// $Id: dilute_gauss_src_s.cc,v 3.0 2006-04-03 04:59:06 edwards Exp $
/*! \file
 *  \brief Variety of Z2 noise sources
 */

#include "chromabase.h"
#include "meas/hadron/stag_propShift_s.h"

namespace Chroma {

//! Volume source of complex Z2 noise
/*!
 * \ingroup hadron
 *
 * This routine is specific to Wilson fermions! 
 *
 * \param a      Source fermion
 *
 *
 *  This type of source is required to compute disconnected
 *  diagrams. The source is complex Z2 noise, hence there
 *  is an additional normalization factor of 1/sqrt(2).
 *
 *
 */



//! Timeslice source of complex gaussian noise
/*!
 * \ingroup hadron
 *
 *
 * \param a      Source fermion
 * \param slice        time slice
 * \param mu           direction of slice
 *
 *
 *  This type of source is useful for computing hadronic
 *  decay like diagrams. 
 *
 *  QUARK MASS DEPENDENCE OF HADRON MASSES FROM LATTICE QCD.
 * By UKQCD Collaboration (M. Foster et al.)
 * Published in Phys.Rev.D59:074503,1999
 * e-Print Archive: hep-lat/9810021 
 *
 */

/*************************************************************************/
void gaussian_on_timeslice(LatticeStaggeredFermion& a, int slice, int mu){

  // compute a source of z2 noise on slice x_mu = slice
  LatticeStaggeredFermion tmp; 
  gaussian(tmp) ;
  printf("slice=%d\n",slice);fflush(stdout);
  a = where(Layout::latticeCoordinate(mu) == slice, tmp, LatticeStaggeredFermion(zero));
}


/*************************************************************************/
void gaussian_on_parity(LatticeStaggeredFermion& a, int parity){

  //"parity" should be passed in as either 1 ---> odd sites
  //                                 or    0 ---> even sites
  

  // compute a volume source of z2 noise
  LatticeStaggeredFermion tmp; 

  gaussian(tmp) ;       // fill tmp with gaussian complexs  

  a = where((((Layout::latticeCoordinate(0) + Layout::latticeCoordinate(1)+
	       Layout::latticeCoordinate(2) + Layout::latticeCoordinate(3))%2)
	     == parity), 
	    tmp, LatticeStaggeredFermion(zero));

}

/*************************************************************************/

void gaussian_color_src(LatticeStaggeredFermion& a, int color_index){

  LatticeComplex          lat_rand;
  LatticeColorVector      latcolor = zero;
  //  LatticeStaggeredFermion latfield = zero;
  const int               spin_index = 0 ;


  a=zero;

  gaussian(lat_rand) ;   // fill c with gaussian rand numbers (both components)

  pokeSpin(a, pokeColor(latcolor, lat_rand, color_index), spin_index);

}



/*************************************************************************/
void gaussian_color_src_on_slice(LatticeStaggeredFermion& a, int color_index, 
				 int slice, int mu){

  
  LatticeStaggeredFermion tmp; 
  LatticeComplex          lat_rand;
  LatticeColorVector      latcolor = zero;
  const int               spin_index = 0 ;


  a=zero;                    // zero out the fermion source, use zero, not 0

  gaussian(lat_rand) ;       // put a gaussian rand on all sites of lat_rand

  pokeSpin(tmp, pokeColor(latcolor, lat_rand, color_index), spin_index);

  //  pokeColor(tmp, lat_rand, color); // put lat_rand on color component of tmp

  //copy tmp into a on all sites in timeslice x_mu=slice
  a = where(Layout::latticeCoordinate(mu) == slice, tmp, LatticeStaggeredFermion(zero));

}
/*************************************************************************/
void gaussian_color_src_on_parity(LatticeStaggeredFermion& a, int color_index,
				  int parity){

  //"parity" should be passed in as either 1 ---> odd sites
  //                                 or    0 ---> even sites
  

  LatticeStaggeredFermion tmp; 
  LatticeComplex          lat_rand;
  LatticeColorVector      latcolor = zero;
  const int               spin_index = 0 ;

  a=zero;                    // zero out the fermion source, use zero, not 0

  gaussian(lat_rand) ;       // put a gaussian rand on all sites of lat_rand

  pokeSpin(tmp, pokeColor(latcolor, lat_rand, color_index), spin_index);

 //  pokeColor(tmp, lat_rand, color); // put lat_rand on color component of tmp

  //copy tmp into a on all sites with given parity
  a = where((((Layout::latticeCoordinate(0) + Layout::latticeCoordinate(1) +
	       Layout::latticeCoordinate(2)+Layout::latticeCoordinate(3))%2)
	     == parity), 
	    tmp, LatticeStaggeredFermion(zero));

}
/*************************************************************************/
void gaussian_parity_src_on_slice(LatticeStaggeredFermion& a,
			       int parity, int slice, int mu){

  //"parity" should be passed in as either 1 ---> odd sites
  //                                 or    0 ---> even sites
  

  LatticeStaggeredFermion tmp; 
  LatticeComplex          lat_rand;
  LatticeColorVector      latcolor = zero;
  const int               spin_index = 0 ;

  a=zero;                    // zero out the fermion source, use zero, not 0

  gaussian(tmp) ;       // put a gaussian rand on all sites of tmp


 //  pokeColor(tmp, lat_rand, color); // put lat_rand on color component of tmp

  //copy tmp into a on all sites with given parity
  a = where(((((Layout::latticeCoordinate(0) + Layout::latticeCoordinate(1) +
	       Layout::latticeCoordinate(2)+Layout::latticeCoordinate(3))%2)
	      == parity) && (Layout::latticeCoordinate(mu) == slice)), 
	    tmp, LatticeStaggeredFermion(zero));

}
/*************************************************************************/
void gaussian_on_mod_timeslice(LatticeStaggeredFermion& a, int slice, int mu,
			       int seperation){

  // compute a source of z2 noise on slice x_mu % seperation = slice
  LatticeStaggeredFermion tmp; 
  gaussian(tmp) ;
  printf("slice=%d\n",slice);fflush(stdout);
  //  a = where((Layout::latticeCoordinate(mu))%seperation == slice, tmp, 
  //	    LatticeStaggeredFermion(zero));

  a = where(((Layout::latticeCoordinate(mu))-slice)%seperation == 0, tmp, 
	    LatticeStaggeredFermion(zero));
}

/*************************************************************************/
void gaussian_on_corner(LatticeStaggeredFermion& a, int corner_index){

  // define a gaussian source on one corner of all hypercubes

  multi1d<int> coord(Nd);

  PropIndexTodelta(corner_index, coord) ;

  

  // compute a volume source of z2 noise
  LatticeStaggeredFermion tmp; 

  gaussian(tmp) ;       // fill tmp with gaussian complexs  

  //  a = where((
  //	     (
  //	      ((Layout::latticeCoordinate(0))%2==coord[0]) && 
  //	      ((Layout::latticeCoordinate(1))%2==coord[1])) &&
  //	     (
  //	      ((Layout::latticeCoordinate(2))%2==coord[2]) &&
  //	      ((Layout::latticeCoordinate(3))%2==coord[3]))),
  //	    tmp, LatticeStaggeredFermion(zero));

  a = where((
	     (
	      (((Layout::latticeCoordinate(0))-coord[0])%2==0) && 
	      (((Layout::latticeCoordinate(1))-coord[1])%2==0)) &&
	     (
	      (((Layout::latticeCoordinate(2))-coord[2])%2==0) &&
	      (((Layout::latticeCoordinate(3))-coord[3])%2==0))),
	    tmp, LatticeStaggeredFermion(zero));

}

/*************************************************************************/
void gaussian_corner_on_dbl_slice(LatticeStaggeredFermion& a, 
				   int corner_index,
				   int slice, int mu){

  // define a gaussian source on one corner of all hypercubes on 
  // double-timeslice==slice.
  // double-timeslice = timeslice/2

  multi1d<int> coord(Nd);

  PropIndexTodelta(corner_index, coord) ;

  

  // compute a volume source of z2 noise
  LatticeStaggeredFermion tmp; 

  gaussian(tmp) ;       // fill tmp with gaussian complexs  

  //  a = where(((
  //	      (
  //	       ((Layout::latticeCoordinate(0))%2==coord[0]) && 
  //	       ((Layout::latticeCoordinate(1))%2==coord[1])) &&
  //	      (
  //	       ((Layout::latticeCoordinate(2))%2==coord[2]) &&
  //	       ((Layout::latticeCoordinate(3))%2==coord[3]))) &&
  //	     (Layout::latticeCoordinate(mu))/2== slice),
  //	    tmp, LatticeStaggeredFermion(zero));

  a = where(((
	      (
	       (((Layout::latticeCoordinate(0))-coord[0])%2==0) && 
	       (((Layout::latticeCoordinate(1))-coord[1])%2==0)) &&
	      (
	       (((Layout::latticeCoordinate(2))-coord[2])%2==0) &&
	       (((Layout::latticeCoordinate(3))-coord[3])%2==0))) &&
	     (Layout::latticeCoordinate(mu))/2== slice),
	    tmp, LatticeStaggeredFermion(zero));

}
/*************************************************************************/
void gaussian_corner_on_mod_dbl_slice(LatticeStaggeredFermion& a, 
				   int corner_index,
				   int slice, int mu, int seperation){

  // define a gaussian source on one corner of all hypercubes on 
  // double-timeslice==slice.
  // double-timeslice is 2 timeslicesw

  multi1d<int> coord(Nd);

  PropIndexTodelta(corner_index, coord) ;
  
  // compute a volume source of z2 noise
  LatticeStaggeredFermion tmp; 

  gaussian(tmp) ;       // fill tmp with gaussian complexs  

  //  a = where(((
  //	      (
  //	       ((Layout::latticeCoordinate(0))%2==coord[0]) && 
  //	       ((Layout::latticeCoordinate(1))%2==coord[1])) &&
  //	      (
  //	       ((Layout::latticeCoordinate(2))%2==coord[2]) &&
  //	       ((Layout::latticeCoordinate(3))%2==coord[3]))) &&
  //	     ((Layout::latticeCoordinate(mu))/2)%seperation == slice),
  //	    tmp, LatticeStaggeredFermion(zero));

  a = where(((
	      (
	       (((Layout::latticeCoordinate(0))-coord[0])%2==0) && 
	       (((Layout::latticeCoordinate(1))-coord[1])%2==0)) &&
	      (
	       (((Layout::latticeCoordinate(2))-coord[2])%2==0) &&
	       (((Layout::latticeCoordinate(3))-coord[3])%2==0))) &&
	     (((Layout::latticeCoordinate(mu))/2-slice)%(seperation) == 0)),
	    tmp, LatticeStaggeredFermion(zero));

}

/*************************************************************************/
void gaussian_color_src_on_mod_slice(LatticeStaggeredFermion& a, 
				     int color_index, int slice, int mu,
				     int seperation){

  
  LatticeStaggeredFermion tmp; 
  LatticeComplex          lat_rand;
  LatticeColorVector      latcolor = zero;
  const int               spin_index = 0 ;


  a=zero;                    // zero out the fermion source, use zero, not 0

  gaussian(lat_rand) ;       // put a gaussian rand on all sites of lat_rand

  pokeSpin(tmp, pokeColor(latcolor, lat_rand, color_index), spin_index);

  //copy tmp into a on all sites in timeslice x_mu%seperation=slice

  //  a = where((Layout::latticeCoordinate(mu)%seperation) == slice, 
  //	    tmp, LatticeStaggeredFermion(zero));

  a = where(((Layout::latticeCoordinate(mu)-slice)%seperation) == 0, 
	    tmp, LatticeStaggeredFermion(zero));

}
/*************************************************************************/

}  // end namespace Chroma

