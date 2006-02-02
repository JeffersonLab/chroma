/*! This function, converts a set of shifts stored in delta, into
 *  an index for the array of propagators. 
 */

#include "chromabase.h"
#include "meas/hadron/stag_propShift_s.h"
#include "meas/smear/displacement.h"

namespace Chroma {

int deltaToPropIndex(multi1d<int>& delta)
{

  if( Nd != 4 ) {
    QDPIO::cerr << "This routine only works for Nd=4, not Nd="<<Nd << endl;
    QDP_abort(1);
  }

  /* Paranoia -- Delta has to have Nd elements */
  if( delta.size() != Nd ) {
    QDPIO::cerr << "Delta has to have Nd elements as opposed to " << delta.size();
    QDP_abort(1);
  }

  /* Paranoia2 -- Delta can contain only 0 (no shift) and 1s (shift) */
  for( int mu = 0; mu < Nd; mu++) {
    if( delta[mu] != 0 && delta[mu] != 1 ) {
      QDPIO::cerr << "Delta must be made up of either zeros or ones. Element " << mu << " has value " << delta[mu] << endl;
      QDP_abort(1);
    }
  }

  /* The delta array is turned into a lexicographic index within
     the 4d hypercube */
  return(delta[0] + 2*( delta[1] + 2*(delta[2] + 2*delta[3])));
}




void PropIndexTodelta(int src_index, multi1d<int>& delta)
{

  if( Nd != 4 ) {
    QDPIO::cerr << "This routine only works for Nd=4, not Nd="<<Nd << endl;
    QDP_abort(1);
  }

  /* Paranoia -- Delta has to have Nd elements */
  if( delta.size() != Nd ) {
    QDPIO::cerr << "Delta has to have Nd elements as opposed to " << delta.size();
    QDP_abort(1);
  }

  /* Paranoia1 -- src_index must lie in range */
  if( src_index <  0 || src_index  > 15  ) {
    QDPIO::cerr << "src_index out of range " << src_index << "\n" ;
    QDP_abort(1);
  }

  //
  delta[0] = src_index % 2 ; 

  int tmp = ( src_index - delta[0] ) / 2 ; 
  delta[1] = tmp % 2 ; 

  int tmpA = ( tmp - delta[1] ) / 2 ; 
  delta[2] = tmpA % 2 ; 

  tmp = ( tmpA - delta[2] ) / 2 ; 
  delta[3] = tmp % 2 ; 


  /* Paranoia3 -- Delta can contain only 0 (no shift) and 1s (shift) */
  for( int mu = 0; mu < Nd; mu++) {
    if( delta[mu] != 0 && delta[mu] != 1 ) {
      QDPIO::cerr << "Delta must be made up of either zeros or ones. Element " << mu << " has value " << delta[mu] << endl;
      QDP_abort(1);
    }
  }

  return ;
}


/*! Given an array of forward shifts (up to 1 in each dimension) and a
 *  propagator. This routine will carry out up to 4 shifts on the input,
 *  possibly 1 in each dimension
 */

LatticeStaggeredPropagator shiftDeltaProp(multi1d<int>& delta,
                                 const LatticeStaggeredPropagator& src)
{

  int mu;

  if( delta.size() != Nd ) {
    QDPIO::cerr << "Delta has to have Nd elements as opposed to " << delta.size();
    QDP_abort(1);
  }

  for( mu = 0; mu < Nd; mu++) {
    if( delta[mu] != 0 && delta[mu] != 1 ) {
      QDPIO::cerr << "Delta must be made up of either zeros or ones. Element " << mu << " has value " << delta[mu] << endl;
      QDP_abort(1);
    }
  }

  LatticeStaggeredPropagator ret_val = src;
  LatticeStaggeredPropagator tmp;

  for( mu = 0; mu < Nd; mu++) {
    if( delta[mu] == 1 ) {
      // This at the moment cannot occur without a temporary
      tmp = shift(ret_val, FORWARD, mu);
      ret_val = tmp;
    }
  }

  return ret_val;
}




/**********************************************************************/

template<typename T>
T shiftDeltaPropCov_t(multi1d<int>& delta,
		    const T & src,
		    multi1d<LatticeColorMatrix> u, bool sym_flag){
 
 /* This variant of  shiftDeltaProp takes as an argument:
         sym_flag=true for symmetric shifting;
         sym_flag=false for no symmetric shifting;

     This is the gauge covariant shiftDeltaProp for gauge fixed configs.

     Shift averages over all paths between start and end locations.
  */

 
  int mu;                               /* shift direction */
  int num_perms;                        /* number of permutations of shifts */
  int i;                                /* loop index */
  int j;                                /* loop index */
  int fac;                              /* index for computing factorial */
  const int length = 1 ;                /* shift only one lat unit */
  double inv_num_perms;                  /* = 1/num_perms */

  T ret_val=zero; /* returned value */
  T tmp1;      /* temp for holding shifted prop */
  T tmp2;      /* temp for holding shifted prop */

  const int shift_index[24][4]=
    {{0,1,2,3},{1,0,2,3},{0,2,1,3},{2,0,1,3},{1,2,0,3},{2,1,0,3},
     {0,1,3,2},{1,0,3,2},{0,3,1,2},{3,0,1,2},{1,3,0,2},{3,1,0,2},
     {0,2,3,1},{0,3,2,1},{2,0,3,1},{2,3,0,1},{3,2,0,1},{3,0,2,1},
     {1,2,3,0},{1,3,2,0},{2,1,3,0},{2,3,1,0},{3,1,2,0},{3,2,1,0}};

  /* 
     Some notes about the shift_index. This 24x4 matrix is all of the
     24 permutations of 4 different links, that is, they make up
     the 24 different paths connecting opposite corners of a 4-d
     hypercube.

     The ordering is significant!
     The "upper-left" 6x3 matrix contained within it is the collection
     of all 6 permutations of 3 paths:
     {{0,1,2},{1,0,2},{0,2,1},{2,0,1},{1,2,0},{2,1,0}}

     Likewise the "upper-left" 2x2 matrix contained within is
     the collection of both 2-link paths connecting opposite corners of 
     a square:
     {{0,1},{1,0}}

     Of course the first element (a 1x1) is the one 1-link path between 
     two adjacent vertices: 
     {{0}}

     Of course this little hack crashes and burns when Nd>4

  */

  /*debug*/ 
  /*
      if(sym_flag){
	printf("SYM SHIFT!\n");
      }else{
	printf("NOT SYM SHIFT!\n");
      }
  */
  if( Nd != 4 ) {
    QDPIO::cerr << "This only works for Nd=4 , not Nd=" << Nd;
    QDP_abort(1);
  }

  if( delta.size() != Nd ) {
    QDPIO::cerr << "Delta has to have Nd elements as opposed to " << delta.size();
    QDP_abort(1);
  }

  for( mu = 0; mu < Nd; mu++) {
    if( delta[mu] != 0 && delta[mu] != 1 ) {
      QDPIO::cerr << "Delta must be made up of either zeros or ones. Element " << mu << " has value " << delta[mu] << endl;
      QDP_abort(1);
    }
  }

  int num_shifts = 0 ; 
  multi1d<int>  delta_order(4);
  delta_order = 0 ;

  for( mu = 0; mu < Nd; mu++) {
    if(delta[mu] ==1 ){
      delta_order[num_shifts] = mu ;
      ++num_shifts ;
    }
  }


  /* find the number of permutations = factorial(num_shifts) */
  num_perms=1;
  fac=num_shifts;

  while(fac>1){
    num_perms*=fac;
    fac--;
  }
  inv_num_perms= 1.0/((double)num_perms);

  ret_val=zero;                          /* initialize returned value */

  for(i=0; i<num_perms; i++){
    tmp1=src;
    tmp2=src;
    for(j=0; j<num_shifts; j++){

      displacement(u,tmp1,length, delta_order[shift_index[i][j]]);

      if(sym_flag){

	displacement(u,tmp2,-length, delta_order[shift_index[i][j]]);

	tmp1+=tmp2;
	tmp1*=0.5;
	tmp2=tmp1;       
      }
    }

    ret_val+=tmp1;

  }
  /* normalize by the number of permutations */
  ret_val*=inv_num_perms;

  return ret_val;
}

LatticeStaggeredPropagator shiftDeltaPropCov(multi1d<int>& delta,
			   const LatticeStaggeredPropagator& src,
			   multi1d<LatticeColorMatrix> u, bool sym_flag)
{

  return shiftDeltaPropCov_t<LatticeStaggeredPropagator>(delta,src,u,
						       sym_flag);
}



LatticeStaggeredFermion shiftDeltaPropCov(multi1d<int>& delta,
			   const LatticeStaggeredFermion& src,
			   multi1d<LatticeColorMatrix> u, bool sym_flag)
{

  LatticeStaggeredFermion tmp = 
  shiftDeltaPropCov_t<LatticeStaggeredFermion>(delta,src,u,
					       sym_flag);

  //  LatticeStaggeredFermion tmp =  zero ; // DEBUG

 return tmp ; 
}





/*************************************************************************/
LatticeStaggeredPropagator shiftDeltaProp(multi1d<int>& delta,
                                 const LatticeStaggeredPropagator& src, 
					  bool sym_flag){

  /* This variant of  shiftDeltaProp takes as an argument:
         sym_flag=true for symmetric shifting;
         sym_flag=false for no symmetric shifting;

     This is the non-covariant shiftDeltaProp for gauge fixed configs.

     Maybe should have a test to make sure gauge is fixed???

  */

  int mu;

  if( delta.size() != Nd ) {
    QDPIO::cerr << "Delta has to have Nd elements as opposed to " << delta.size();
    QDP_abort(1);
  }

  for( mu = 0; mu < Nd; mu++) {
    if( delta[mu] != 0 && delta[mu] != 1 ) {
      QDPIO::cerr << "Delta must be made up of either zeros or ones. Element " << mu << " has value " << delta[mu] << endl;
      QDP_abort(1);
    }
  }

  LatticeStaggeredPropagator ret_val = src;
  LatticeStaggeredPropagator tmp1;
  LatticeStaggeredPropagator tmp2;

  for( mu = 0; mu < Nd; mu++) {
    if( delta[mu] == 1 ) {
      // This at the moment cannot occur without a temporary
      tmp1 = shift(ret_val, FORWARD, mu);
      if(sym_flag){
	tmp2 = shift(ret_val, BACKWARD, mu);
      }
      if(sym_flag){
	ret_val = 0.5*(tmp1+tmp2);
      }else{
	ret_val = tmp1;
      }
    }
  }

  return ret_val;
}

}  // end namespace Chroma
