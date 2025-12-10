#ifndef CPP_CLOVER_SITE_APPLY_64_BIT_SSE_H
#define CPP_CLOVER_SITE_APPLY_64_BIT_SSE_H

#include <cpp_clover_types.h>
#include <cpp_dslash_types.h>


#undef CMADD
#undef CONJMADD

#include<xmmintrin.h>
#ifndef DSLASH_USE_SSE3

/* SSE 2 Headers */

#define CMADD(z,x,y)				\
  { \
    __m128d t1,t2,t3,t4; \
    t1 = _mm_mul_pd(x,y); \
    t2 = _mm_shuffle_pd(t1,t1,0x1); \
    t3 = _mm_shuffle_pd(y,y,0x1);\
    t4 = _mm_sub_pd(t1,t2); \
    t2 = _mm_mul_pd(x,t3); \
    t3 = _mm_shuffle_pd(t2,t2,0x1); \
    t3 = _mm_add_pd(t2,t3); \
    t4= _mm_shuffle_pd(t4,t3,0x2); \
    z = _mm_add_pd(z,t4); \
  }
/* z += x*conj(y)    z, x, y are SSE registers containing complex numbers
   ordered with the real part in the low half, imag part 
   in the upper half */
#define CONJMADD(z,x,y)				\
  { \
    __m128d t1,t2,t3,t4; \
    t1 = _mm_mul_pd(x,y); \
    t2 = _mm_shuffle_pd(t1,t1,0x1); \
    t3 = _mm_shuffle_pd(y,y,0x1);\
    t4 = _mm_add_pd(t1,t2); \
    t2 = _mm_mul_pd(x,t3); \
    t3 = _mm_shuffle_pd(t2,t2,0x1); \
    t3 = _mm_sub_pd(t2,t3); \
    t4= _mm_shuffle_pd(t4,t3,0x2); \
    z = _mm_add_pd(z,t4); \
  }

#else 
/* SSE 3 */
#include <pmmintrin.h>
#define CMADD(z,x,y)				\
  { \
    __m128d t1,t2;	      \
    t1 = _mm_mul_pd((x),(y)); \
    t1 = _mm_hsub_pd(t1,t1); \
    t2 = _mm_shuffle_pd((y),(y),0x1);\
    t2 = _mm_mul_pd((x),t2); \
    t2 = _mm_hadd_pd(t2,t2); \
    t1= _mm_shuffle_pd(t1,t2,0x2);		\
    (z) = _mm_add_pd((z),t1);			\
  }


#define CONJMADD(z,x,y)				\
  { \
    __m128d t1,t2; \
    t1 = _mm_mul_pd((x),(y)); \
    t1 = _mm_hadd_pd(t1,t1); \
    t2 = _mm_shuffle_pd((x),(x),0x1);\
    t2 = _mm_mul_pd((y),t2); \
    t2 = _mm_hsub_pd(t2,t2); \
    t1= _mm_shuffle_pd(t1,t2,0x2);		\
    (z) = _mm_add_pd((z),t1);			\
  }

#endif

namespace CPlusPlusClover { 

  // Use 64bit types
  using namespace Clover64BitTypes;
  using namespace Dslash64BitTypes;

  namespace CPlusPlusClover64Bit { 

    inline 
    void cloverSiteApply(FourSpinor result, const CloverTerm clov, const FourSpinor src) 
    {

    __m128d chi0;
    __m128d chi1;
    
    __m128d psi0;
    __m128d psi1;
    __m128d psi2;
    __m128d psi3;
    __m128d psi4;
    __m128d psi5;

    __m128d tmp0;
    __m128d tmp1;



    double* chi_p=(double *)&result[0][0][0];
    double* psi_p=(double *)&src[0][0][0];

    // Block 0
    double* offd_p = (double *)&(clov[0].off_diag[0][0]);
    double* diag_p = (double *)&(clov[0].diag[0]);


    // First things first: Registerize half of psi
    psi0 = _mm_load_pd(psi_p);    // ppsi[0]
    psi1 = _mm_load_pd(psi_p+2);  // ppsi[1]
    psi2 = _mm_load_pd(psi_p+4);  // ppsi[2]
    psi3 = _mm_load_pd(psi_p+6);  // ppsi[3]
    psi4 = _mm_load_pd(psi_p+8);  // ppsi[4]
    psi5 = _mm_load_pd(psi_p+10); // ppsi[5]
    
    // cchi[ 0]  = tri_diag[site][0][ 0]  * ppsi[ 0]
    // cchi[ 1]  = tri_diag[site][0][ 1]  * ppsi[ 1]
    
    // tmp0 = [ diag[1] | diag[0] ]
    tmp0 = _mm_load_pd(diag_p);  //  reads tri_diag[site][0][0]
    //  AND   tri_diat[site][0][1]
    // tmp1 = [ diag[0] | diag[1] ]
    tmp1 = _mm_shuffle_pd(tmp0,tmp0,0x1);
    // tmp0 = [ diag[0] | diag[0] ]
    tmp0 = _mm_shuffle_pd(tmp0,tmp0, 0x0);
    // tmp1 = [ diag[1] | diag[1] ]
    tmp1 = _mm_shuffle_pd(tmp1,tmp1, 0x0);
    
    // cchi[ 0] = tri_diag[site][0][ 0]  * ppsi[ 0]
    chi0 = _mm_mul_pd(tmp0,psi0);
    
    //  cchi[ 1] = tri_diag[site][0][ 1]  * ppsi[ 1]
    chi1 = _mm_mul_pd(tmp1,psi1);
    
    // cchi[ 0] +=  conj(tri_off_diag[site][0][ 0]) * ppsi[ 1]
    // cchi[ 1] +=       tri_off_diag[site][0][ 0]  * ppsi[ 0]
    tmp0 = _mm_load_pd(offd_p); // tri_off_diag[site][0][ 0]
    CONJMADD(chi0,psi1,tmp0);
    CMADD(chi1, psi0, tmp0);
    
    // cchi[ 0] +=  conj(tri_off_diag[site][0][ 1]) * ppsi[ 2]
    // cchi[ 1] +=  conj(tri_off_diag[site][0][ 2]) * ppsi[ 2]
    tmp0 = _mm_load_pd(offd_p+2); 
    tmp1 = _mm_load_pd(offd_p+4); 
    CONJMADD(chi0, psi2, tmp0); 
    CONJMADD(chi1, psi2, tmp1); 
    
    // cchi[ 0] +=  conj(tri_off_diag[site][0][ 3]) * ppsi[ 3]
    // cchi[ 1] +=  conj(tri_off_diag[site][0][ 4]) * ppsi[ 3]
    tmp0 = _mm_load_pd(offd_p+6); // tri_off_diag[site][0][3]
    tmp1 = _mm_load_pd(offd_p+8); // tri_off_diag[site][0][4]
    CONJMADD(chi0, psi3, tmp0); 
    CONJMADD(chi1, psi3, tmp1); 
    
    // cchi[ 0] +=  conj(tri_off_diag[site][0][ 6]) * ppsi[ 4]
    // cchi[ 1] +=  conj(tri_off_diag[site][0][ 7]) * ppsi[ 4]
    tmp0 = _mm_load_pd(offd_p+12); // tri_off_diag[site][0][6]
    tmp1 = _mm_load_pd(offd_p+14); // tri_off_diag[site][0][7]
    CONJMADD(chi0, psi4,tmp0);
    CONJMADD(chi1, psi4,tmp1); 
    
    // cchi[ 0] +=  conj(tri_off_diag[site][0][10]) * ppsi[ 5];
    // cchi[ 1] +=  conj(tri_off_diag[site][0][11]) * ppsi[ 5];
    tmp0 = _mm_load_pd(offd_p+20); // tri_off_diag[site][0][10]
    tmp1 = _mm_load_pd(offd_p+22); // tri_off_diag[site][0][11]
    CONJMADD(chi0, psi5,tmp0);
    CONJMADD(chi1, psi5,tmp1);
    
    _mm_store_pd(chi_p, chi0);
    _mm_store_pd(chi_p+2, chi1);
    
    
    //    cchi[ 2]  =  tri_diag[site][0][ 2]  * ppsi[ 2]
    //    cchi[ 3]  =  tri_diag[site][0][ 3]  * ppsi[ 3]
    
    tmp0 = _mm_load_pd(diag_p+2); //  reads tri_diag[site][0][2]
    //  AND   tri_diat[site][0][3]
    
    // tmp1 = [ diag[2] | diag[3] ]
    tmp1 = _mm_shuffle_pd(tmp0,tmp0,0x1);
    // tmp0 = [ diag[2] | diag[2] ]
    tmp0 = _mm_shuffle_pd(tmp0,tmp0, 0x0);
    // tmp1 = [ diag[3] | diag[3] ]
    tmp1 = _mm_shuffle_pd(tmp1,tmp1, 0x0);
    
    // cchi[ 2] = tri_diag[site][0][ 2]  * ppsi[2]
    chi0 = _mm_mul_pd(tmp0,psi2);
    
    // cchi[ 3] = tri_diag[site][0][ 3]  * ppsi[ 3]
    chi1 = _mm_mul_pd(tmp1,psi3);    
    
    //    cchi[ 2] +=  tri_off_diag[site][0][ 1]  * ppsi[ 0]
    //    cchi[ 2] +=  tri_off_diag[site][0][ 2]  * ppsi[ 1]
    tmp0 = _mm_load_pd(offd_p+2); // tri_off_diag[site][0][ 1]
    tmp1 = _mm_load_pd(offd_p+4); // tri_off_diag[site][0][ 2]
    CMADD(chi0,psi0,tmp0);
    CMADD(chi0,psi1,tmp1);
    
    //    cchi[ 3] +=  tri_off_diag[site][0][ 3]  * ppsi[ 0]
    //    cchi[ 3] +=  tri_off_diag[site][0][ 4]  * ppsi[ 1]
    tmp0 = _mm_load_pd(offd_p+6);
    tmp1 = _mm_load_pd(offd_p+8);
    CMADD(chi1,psi0,tmp0);
    CMADD(chi1,psi1,tmp1);
    
    //    cchi[ 2] +=  conj(tri_off_diag[site][0][ 5]) * ppsi[ 3]
    //    cchi[ 3] +=       tri_off_diag[site][0][ 5]  * ppsi[ 2]
    tmp0 = _mm_load_pd(offd_p+10); // tri_off_diag[site][0][ 5]
    CONJMADD(chi0,psi3,tmp0);
    CMADD(chi1,psi2, tmp0);
    
    //    cchi[ 2] +=  conj(tri_off_diag[site][0][ 8]) * ppsi[ 4]
    //    cchi[ 3] +=  conj(tri_off_diag[site][0][ 9]) * ppsi[ 4]
    tmp0 = _mm_load_pd(offd_p+16); // tri_off_diag[site][0][ 8]
    tmp1 = _mm_load_pd(offd_p+18); // tri_off_diag[site][0][ 9]
    CONJMADD(chi0,psi4,tmp0);
    CONJMADD(chi1,psi4,tmp1);
    
    //    cchi[ 2] +=  conj(tri_off_diag[site][0][12]) * ppsi[ 5];
	//    cchi[ 3] +=  conj(tri_off_diag[site][0][13]) * ppsi[ 5];
    tmp0 = _mm_load_pd(offd_p+24); // tri_off_diag[site][0][12]
    tmp1 = _mm_load_pd(offd_p+26); // tri_off_diag[site][0][13]
    CONJMADD(chi0,psi5,tmp0);
    CONJMADD(chi1,psi5,tmp1);
    
    _mm_store_pd(chi_p+4, chi0);
    _mm_store_pd(chi_p+6, chi1);	       
    
    
    // cchi[ 4] = tri_diag[site][0][ 4]  * ppsi[ 4]
    // cchi[ 5] = tri_diag[site][0][ 5]  * ppsi[ 5]
    
    // tmp 0 = [ diag[5] | diag[4] ]_low
    tmp0 = _mm_load_pd(diag_p+4); //  reads tri_diag[site][0][4]
    //  AND   tri_diat[site][0][5]
    
    // tmp1 = [ diag[4] | diag[5] ]
    tmp1 = _mm_shuffle_pd(tmp0,tmp0,0x1);
    // tmp0 = [ diag[4] | diag[4] ]
    tmp0 = _mm_shuffle_pd(tmp0,tmp0, 0x0);
    // tmp1 = [ diag[5] | diag[5] ]
    tmp1 = _mm_shuffle_pd(tmp1,tmp1, 0x0);
    
    // cchi[ 4] = tri_diag[site][0][4]  * ppsi[4]
    chi0 = _mm_mul_pd(tmp0,psi4);
    
    //  cchi[ 5] = tri_diag[site][0][5]  * ppsi[5]
    chi1 = _mm_mul_pd(tmp1,psi5);
    
    // cchi[ 4] +=        tri_off_diag[site][0][ 6]  * ppsi[ 0]
    // cchi[ 4] +=        tri_off_diag[site][0][ 7]  * ppsi[ 1]
    tmp0 = _mm_load_pd(offd_p+12);
    tmp1 = _mm_load_pd(offd_p+14);
    CMADD(chi0,psi0,tmp0);
    CMADD(chi0,psi1,tmp1);
    
    // cchi[ 4] +=        tri_off_diag[site][0][ 8]  * ppsi[ 2]
    // cchi[ 4] +=        tri_off_diag[site][0][ 9]  * ppsi[ 3]
    tmp0 = _mm_load_pd(offd_p+16);
    tmp1 = _mm_load_pd(offd_p+18);
    CMADD(chi0,psi2,tmp0);
    CMADD(chi0,psi3,tmp1);
    
    // cchi[ 5] +=        tri_off_diag[site][0][10]  * ppsi[ 0]
    // cchi[ 5] +=        tri_off_diag[site][0][11]  * ppsi[ 1]
    tmp0 = _mm_load_pd(offd_p+20);
    tmp1 = _mm_load_pd(offd_p+22);
    CMADD(chi1,psi0,tmp0);
    CMADD(chi1,psi1,tmp1);
    
    // cchi[ 5] +=        tri_off_diag[site][0][12]  * ppsi[ 2]
    // cchi[ 5] +=        tri_off_diag[site][0][13]  * ppsi[ 3]
    tmp0 = _mm_load_pd(offd_p+24);
    tmp1 = _mm_load_pd(offd_p+26);
    CMADD(chi1,psi2,tmp0);
    CMADD(chi1,psi3,tmp1);
    
    // cchi[ 4] +=   conj(tri_off_diag[site][0][14]) * ppsi[ 5];
    // cchi[ 5] +=        tri_off_diag[site][0][14]  * ppsi[ 4];
    tmp0 = _mm_load_pd(offd_p+28);
    CONJMADD(chi0,psi5,tmp0);
    CMADD(chi1,psi4,tmp0);
    
    _mm_store_pd(chi_p+8, chi0);
    _mm_store_pd(chi_p+10,chi1);


    chi_p=(double *)&result[2][0][0];
    psi_p=(double *)&src[2][0][0];
    offd_p = (double *)&(clov[1].off_diag[0][0]);
    diag_p = (double *)&(clov[1].diag[0]);
	
    

    // First things first: Registerize half of psi
    psi0 = _mm_load_pd(psi_p);    // ppsi[0]
    psi1 = _mm_load_pd(psi_p+2);  // ppsi[1]
    psi2 = _mm_load_pd(psi_p+4);  // ppsi[2]
    psi3 = _mm_load_pd(psi_p+6);  // ppsi[3]
    psi4 = _mm_load_pd(psi_p+8);  // ppsi[4]
    psi5 = _mm_load_pd(psi_p+10); // ppsi[5]
    
    // cchi[ 0]  = tri_diag[site][0][ 0]  * ppsi[ 0]
    // cchi[ 1]  = tri_diag[site][0][ 1]  * ppsi[ 1]
    
    // tmp0 = [ diag[1] | diag[0] ]
    tmp0 = _mm_load_pd(diag_p);  //  reads tri_diag[site][0][0]
    //  AND   tri_diat[site][0][1]
    // tmp1 = [ diag[0] | diag[1] ]
    tmp1 = _mm_shuffle_pd(tmp0,tmp0,0x1);
    // tmp0 = [ diag[0] | diag[0] ]
    tmp0 = _mm_shuffle_pd(tmp0,tmp0, 0x0);
    // tmp1 = [ diag[1] | diag[1] ]
    tmp1 = _mm_shuffle_pd(tmp1,tmp1, 0x0);
    
    // cchi[ 0] = tri_diag[site][0][ 0]  * ppsi[ 0]
    chi0 = _mm_mul_pd(tmp0,psi0);
    
    //  cchi[ 1] = tri_diag[site][0][ 1]  * ppsi[ 1]
    chi1 = _mm_mul_pd(tmp1,psi1);
    
    // cchi[ 0] +=  conj(tri_off_diag[site][0][ 0]) * ppsi[ 1]
    // cchi[ 1] +=       tri_off_diag[site][0][ 0]  * ppsi[ 0]
    tmp0 = _mm_load_pd(offd_p); // tri_off_diag[site][0][ 0]
    CONJMADD(chi0,psi1,tmp0);
    CMADD(chi1, psi0, tmp0);
    
    // cchi[ 0] +=  conj(tri_off_diag[site][0][ 1]) * ppsi[ 2]
    // cchi[ 1] +=  conj(tri_off_diag[site][0][ 2]) * ppsi[ 2]
    tmp0 = _mm_load_pd(offd_p+2); 
    tmp1 = _mm_load_pd(offd_p+4); 
    CONJMADD(chi0, psi2, tmp0); 
    CONJMADD(chi1, psi2, tmp1); 
    
    // cchi[ 0] +=  conj(tri_off_diag[site][0][ 3]) * ppsi[ 3]
    // cchi[ 1] +=  conj(tri_off_diag[site][0][ 4]) * ppsi[ 3]
    tmp0 = _mm_load_pd(offd_p+6); // tri_off_diag[site][0][3]
    tmp1 = _mm_load_pd(offd_p+8); // tri_off_diag[site][0][4]
    CONJMADD(chi0, psi3, tmp0); 
    CONJMADD(chi1, psi3, tmp1); 
    
    // cchi[ 0] +=  conj(tri_off_diag[site][0][ 6]) * ppsi[ 4]
    // cchi[ 1] +=  conj(tri_off_diag[site][0][ 7]) * ppsi[ 4]
    tmp0 = _mm_load_pd(offd_p+12); // tri_off_diag[site][0][6]
    tmp1 = _mm_load_pd(offd_p+14); // tri_off_diag[site][0][7]
    CONJMADD(chi0, psi4,tmp0);
    CONJMADD(chi1, psi4,tmp1); 
    
    // cchi[ 0] +=  conj(tri_off_diag[site][0][10]) * ppsi[ 5];
    // cchi[ 1] +=  conj(tri_off_diag[site][0][11]) * ppsi[ 5];
    tmp0 = _mm_load_pd(offd_p+20); // tri_off_diag[site][0][10]
    tmp1 = _mm_load_pd(offd_p+22); // tri_off_diag[site][0][11]
    CONJMADD(chi0, psi5,tmp0);
    CONJMADD(chi1, psi5,tmp1);
    
    _mm_store_pd(chi_p, chi0);
    _mm_store_pd(chi_p+2, chi1);
    
    
    //    cchi[ 2]  =  tri_diag[site][0][ 2]  * ppsi[ 2]
    //    cchi[ 3]  =  tri_diag[site][0][ 3]  * ppsi[ 3]
    
    tmp0 = _mm_load_pd(diag_p+2); //  reads tri_diag[site][0][2]
    //  AND   tri_diat[site][0][3]
    
    // tmp1 = [ diag[2] | diag[3] ]
    tmp1 = _mm_shuffle_pd(tmp0,tmp0,0x1);
    // tmp0 = [ diag[2] | diag[2] ]
    tmp0 = _mm_shuffle_pd(tmp0,tmp0, 0x0);
    // tmp1 = [ diag[3] | diag[3] ]
    tmp1 = _mm_shuffle_pd(tmp1,tmp1, 0x0);
    
    // cchi[ 2] = tri_diag[site][0][ 2]  * ppsi[2]
    chi0 = _mm_mul_pd(tmp0,psi2);
    
    // cchi[ 3] = tri_diag[site][0][ 3]  * ppsi[ 3]
    chi1 = _mm_mul_pd(tmp1,psi3);    
    
    //    cchi[ 2] +=  tri_off_diag[site][0][ 1]  * ppsi[ 0]
    //    cchi[ 2] +=  tri_off_diag[site][0][ 2]  * ppsi[ 1]
    tmp0 = _mm_load_pd(offd_p+2); // tri_off_diag[site][0][ 1]
    tmp1 = _mm_load_pd(offd_p+4); // tri_off_diag[site][0][ 2]
    CMADD(chi0,psi0,tmp0);
    CMADD(chi0,psi1,tmp1);
    
    //    cchi[ 3] +=  tri_off_diag[site][0][ 3]  * ppsi[ 0]
    //    cchi[ 3] +=  tri_off_diag[site][0][ 4]  * ppsi[ 1]
    tmp0 = _mm_load_pd(offd_p+6);
    tmp1 = _mm_load_pd(offd_p+8);
    CMADD(chi1,psi0,tmp0);
    CMADD(chi1,psi1,tmp1);
    
    //    cchi[ 2] +=  conj(tri_off_diag[site][0][ 5]) * ppsi[ 3]
    //    cchi[ 3] +=       tri_off_diag[site][0][ 5]  * ppsi[ 2]
    tmp0 = _mm_load_pd(offd_p+10); // tri_off_diag[site][0][ 5]
    CONJMADD(chi0,psi3,tmp0);
    CMADD(chi1,psi2, tmp0);
    
    //    cchi[ 2] +=  conj(tri_off_diag[site][0][ 8]) * ppsi[ 4]
    //    cchi[ 3] +=  conj(tri_off_diag[site][0][ 9]) * ppsi[ 4]
    tmp0 = _mm_load_pd(offd_p+16); // tri_off_diag[site][0][ 8]
    tmp1 = _mm_load_pd(offd_p+18); // tri_off_diag[site][0][ 9]
    CONJMADD(chi0,psi4,tmp0);
    CONJMADD(chi1,psi4,tmp1);
    
    //    cchi[ 2] +=  conj(tri_off_diag[site][0][12]) * ppsi[ 5];
	//    cchi[ 3] +=  conj(tri_off_diag[site][0][13]) * ppsi[ 5];
    tmp0 = _mm_load_pd(offd_p+24); // tri_off_diag[site][0][12]
    tmp1 = _mm_load_pd(offd_p+26); // tri_off_diag[site][0][13]
    CONJMADD(chi0,psi5,tmp0);
    CONJMADD(chi1,psi5,tmp1);
    
    _mm_store_pd(chi_p+4, chi0);
    _mm_store_pd(chi_p+6, chi1);	       
    
    
    // cchi[ 4] = tri_diag[site][0][ 4]  * ppsi[ 4]
    // cchi[ 5] = tri_diag[site][0][ 5]  * ppsi[ 5]
    
    // tmp 0 = [ diag[5] | diag[4] ]_low
    tmp0 = _mm_load_pd(diag_p+4); //  reads tri_diag[site][0][4]
    //  AND   tri_diat[site][0][5]
    
    // tmp1 = [ diag[4] | diag[5] ]
    tmp1 = _mm_shuffle_pd(tmp0,tmp0,0x1);
    // tmp0 = [ diag[4] | diag[4] ]
    tmp0 = _mm_shuffle_pd(tmp0,tmp0, 0x0);
    // tmp1 = [ diag[5] | diag[5] ]
    tmp1 = _mm_shuffle_pd(tmp1,tmp1, 0x0);
    
    // cchi[ 4] = tri_diag[site][0][4]  * ppsi[4]
    chi0 = _mm_mul_pd(tmp0,psi4);
    
    //  cchi[ 5] = tri_diag[site][0][5]  * ppsi[5]
    chi1 = _mm_mul_pd(tmp1,psi5);
    
    // cchi[ 4] +=        tri_off_diag[site][0][ 6]  * ppsi[ 0]
    // cchi[ 4] +=        tri_off_diag[site][0][ 7]  * ppsi[ 1]
    tmp0 = _mm_load_pd(offd_p+12);
    tmp1 = _mm_load_pd(offd_p+14);
    CMADD(chi0,psi0,tmp0);
    CMADD(chi0,psi1,tmp1);
    
    // cchi[ 4] +=        tri_off_diag[site][0][ 8]  * ppsi[ 2]
    // cchi[ 4] +=        tri_off_diag[site][0][ 9]  * ppsi[ 3]
    tmp0 = _mm_load_pd(offd_p+16);
    tmp1 = _mm_load_pd(offd_p+18);
    CMADD(chi0,psi2,tmp0);
    CMADD(chi0,psi3,tmp1);
    
    // cchi[ 5] +=        tri_off_diag[site][0][10]  * ppsi[ 0]
    // cchi[ 5] +=        tri_off_diag[site][0][11]  * ppsi[ 1]
    tmp0 = _mm_load_pd(offd_p+20);
    tmp1 = _mm_load_pd(offd_p+22);
    CMADD(chi1,psi0,tmp0);
    CMADD(chi1,psi1,tmp1);
    
    // cchi[ 5] +=        tri_off_diag[site][0][12]  * ppsi[ 2]
    // cchi[ 5] +=        tri_off_diag[site][0][13]  * ppsi[ 3]
    tmp0 = _mm_load_pd(offd_p+24);
    tmp1 = _mm_load_pd(offd_p+26);
    CMADD(chi1,psi2,tmp0);
    CMADD(chi1,psi3,tmp1);
    
    // cchi[ 4] +=   conj(tri_off_diag[site][0][14]) * ppsi[ 5];
    // cchi[ 5] +=        tri_off_diag[site][0][14]  * ppsi[ 4];
    tmp0 = _mm_load_pd(offd_p+28);
    CONJMADD(chi0,psi5,tmp0);
    CMADD(chi1,psi4,tmp0);
    
    _mm_store_pd(chi_p+8, chi0);
    _mm_store_pd(chi_p+10,chi1);


      }

  } // Namespace
}  // Namespace

#undef CMADD
#undef CONJMADD


#endif
