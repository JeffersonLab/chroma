#include <cstdlib>
#include <iostream>

#include <cpp_dslash_scalar.h>
#include <shift_table_scalar.h>
#include <dispatch_scalar.h>

#include <cpp_clover_scalar.h>
#include <cpp_dslash_scalar_64bit.h>
#include <cpp_clover_site_apply_64bit.h>
#include "allocate.h"
using namespace CPlusPlusWilsonDslash::DslashScalar64Bit;
using namespace CPlusPlusWilsonDslash::Dslash64BitTypes;
using namespace CPlusPlusClover::Clover64BitTypes;

namespace CPlusPlusClover {
 

  /* Constructor */
  CloverSchur4D<double>::CloverSchur4D(const int latt_size[],      
				       void (*getSiteCoords)(int coord[], int node, int linearsite),
				       int (*getLinearSiteIndex)(const int coord[]),
				       int (*nodeNum)(const int coord[])
				       ) 
  {
    s = new ShiftTable(latt_size,
		       getSiteCoords,
		       getLinearSiteIndex,
		       nodeNum);

    xt_spinor = (FourSpinor *)CPlusPlusWilsonDslash::alloc(2*s->totalVolCB()*sizeof(Dslash64BitTypes::FourSpinor)
				    +Cache::CacheLineSize);
    if( xt_spinor == (FourSpinor *)NULL ) { 
      std::cerr << "Unable to allocate temporary" << std::endl;
      exit(1);
    }
    unsigned long pad = 0;
    if ( (unsigned long)xt_spinor % Cache::CacheLineSize != 0 ) { 
      pad = Cache::CacheLineSize - ( (unsigned long)xt_spinor % Cache::CacheLineSize );
    }

    t_spinor = (FourSpinor*)((unsigned char *)xt_spinor + pad);
#if 0
    std::cerr << "Sizeof FourSpinor=" << sizeof(FourSpinor) << std::endl;
    std::cerr << "t_spinor = " << t_spinor << std::endl;
    std::cerr << "xt_spinor = " << xt_spinor << std::endl;
    std::cerr << "Cache::CacheLineSize = " << Cache::CacheLineSize << std::endl;
    std::cerr << " xt_spinor % Cache::CacheLineSize " << (unsigned long)xt_spinor % (unsigned long)Cache::CacheLineSize << std::endl;
    std::cerr << "pad=" << pad << std::endl;
#endif
  }

  CloverSchur4D<double>::~CloverSchur4D() { delete s; CPlusPlusWilsonDslash::dealloc(xt_spinor); }

  namespace CloverScalar64Bit {

   
    void ClovDPsiPlus(size_t lo, size_t hi, int id, const void *ptr);
    void ClovDPsiMinus(size_t lo, size_t hi, int id, const void *ptr);

  }

  using  namespace CloverScalar64Bit;
  // The operator 
  void CloverSchur4D<double>::operator() (double* res, 
					  const double* psi, 
					  const double *u, /* Gauge field suitably packed */
				 
					  const double *clov_oo,
					  const double *invclov_ee,
				   int isign)

  {
    if (isign == 1) {  

      CPlusPlusClover::dispatchToThreads((void (*)(size_t, size_t, int, const void *))&ClovDPsiPlus, 
					 (void*)psi,
					 (void*)res,
					 (void *)u,
					 (void *)invclov_ee, 
					 (void *)clov_oo, 
					 (void *)t_spinor, 
					 (void *)s,
					 s->totalVolCB());
      

    }

    if( isign == -1) {

      CPlusPlusClover::dispatchToThreads((void (*)(size_t, size_t, int, const void *))&ClovDPsiMinus, 
					 (void *)psi,
					 (void *)res,
					 (void *)u,
					 (void *)invclov_ee,
					 (void *)clov_oo,
					 (void *)t_spinor, 
					 (void *)s, 
					 (int)s->totalVolCB());
      
    }    
  }

  using namespace DslashScalar64Bit;
  using namespace CPlusPlusClover64Bit;
  namespace CloverScalar64Bit {

    void ClovDPsiPlus(size_t lo, size_t hi, int id, const void *ptr)
    {
      int ix,iy,iz;                                      /* Ix is corrent site */
      /* Iy is a neighbour */
      /* Iz is next site in loop */
      
      const CloverThreadWorkerArgs *a = (const CloverThreadWorkerArgs*)ptr;              

      ShiftTable *shift = (ShiftTable *)a->s;
      int total_vol_cb = shift->totalVolCB();

      GaugeMatrix (*gauge_field)[4] ALIGN = (GaugeMatrix(*)[4])a->u;        /* My gauge field */
      FourSpinor *psi = (FourSpinor *)a->psi;                        /* Source */
      FourSpinor *res = (FourSpinor *)a->res;                        /* Result */
      
      FourSpinor *t_spinor = (FourSpinor *)a->t_spinor;   /* Temporary spinor */
      CloverTerm* invclov = (CloverTerm *)a->invclov_ee;  /* Clover Inverse */
      CloverTerm* clov = (CloverTerm *)a->clov_oo;        /* Clover term */

      GaugeMatrix *up,*um;                               /* Pointer to FORWARD neighbour */
      FourSpinor *s,*sp,*sm,*rn;                       /* Pointer to BACKWARD neighbour */
      
      FourSpinor d_oe_psi __attribute__((aligned(16)));
				       

      int thissite;
	 
      
      /*********** loop over all lattice sites *************************/
      for (ix=lo;ix<hi;ix++) {
	/************** direction +0 ***********/
	thissite=shift->siteTable(ix);
	
	/* This is like a prefetch 
	   - we peel it off the loop */
	
	/* Get 4 spinor from forward direction */
	sp=&psi[shift->forwardNeighbor(ix,0) ];
	
	
	/* Get Gauge Field */
	up=&(gauge_field[thissite][0]);
	
	
	/******************************   Direction 0  ********************* */
	/* Prefetch back spinor for next dir: -0 */
	iy=shift->backwardNeighbor(ix,0);
	sm=&psi[iy];
	
	
	/* Prefetch back gauge field */
	um=&(gauge_field[iy][0]);
	
	dslash_plus_dir0_forward(*sp,*up,d_oe_psi);
	
	
	/* sm and um should already be prefetched */
	
	/* Now prefetch forward neighbor for next dir (1+) */
	/* And gauge field */
	sp=&psi[ shift->forwardNeighbor(ix,1) ];
	
	up =&(gauge_field[thissite][1]);
	
	dslash_plus_dir0_backward_add(*sm,*um,d_oe_psi);
	
	
	/********************** Direction 1 ************************ */
	/* up and sp should be prefetched */
	
	/* Prefetch backwards spinor and gauge field */
	iy=shift->backwardNeighbor(ix,1);
	sm=&psi[iy];
	
	um=&(gauge_field[iy][1]);
	
	
	
	dslash_plus_dir1_forward_add(*sp,*up,d_oe_psi);
	
	
	
	/* Prefetch forwards spinor and gauge field for next dir: 2+ */
	iy=shift->forwardNeighbor(ix,2);
	sp=&psi[iy];
	
	up = &(gauge_field[thissite][2]);
	
	
	
	dslash_plus_dir1_backward_add(*sm,*um,d_oe_psi);
	
	
	
	/********************** Direction 2 ************************* */
	/* Prefetch back spinor and gauge field  */
	
	iy=shift->backwardNeighbor(ix,2);
	sm=&psi[iy];
	
	um=&(gauge_field[iy][2]);
	
	
	
	dslash_plus_dir2_forward_add(*sp,*up,d_oe_psi);
	
	
	/* Prefetch forward spinor and gauge field for next direction: 3+ */
	iy=shift->forwardNeighbor(ix,3);
	sp=&psi[iy];
	
	up = &(gauge_field[thissite][3]);
	
	
	dslash_plus_dir2_backward_add(*sm,*um,d_oe_psi);
	
	
	/* ******************* Direction +3 **************************** */
	/* Prefetch back spinor and gauge field  3- */
	
	iy=shift->backwardNeighbor(ix,3);
	sm=&psi[iy];
	
	um=&(gauge_field[iy][3]);
	
	
	dslash_plus_dir3_forward_add(*sp,*up,d_oe_psi);      
	
	dslash_plus_dir3_backward_add_store(*sm,*um,d_oe_psi); 
	cloverSiteApply(t_spinor[thissite], invclov[thissite], d_oe_psi);
      }
      

      
      // Other checkerboard
      int low = lo + total_vol_cb;
      int high = hi + total_vol_cb;
      
	 
      /*********** loop over all lattice sites *************************/

      for (int ix = low; ix < high; ix++) { 
	// NOTE: low and high are lo and hi on opposite cb
	thissite=shift->siteTable(ix);
      
	/* Get 4 spinor from forward direction */
	sp=&t_spinor[ shift->forwardNeighbor(ix,0) ];
	
	/* Get Gauge Field */
	up=&(gauge_field[thissite][0]);
	
	/************** direction +0 ***********/
	//	     rn=&D_psi_oe[thissite];
	
	/******************************   Direction 0  ********************* */
	/* Prefetch back spinor for next dir: -0 */
	iy=shift->backwardNeighbor(ix,0);
	sm=&t_spinor[iy];

	     
	/* Prefetch back gauge field */
	um=&(gauge_field[iy][0]);
	dslash_plus_dir0_forward(*sp,*up,d_oe_psi);
	     
	     
	/* sm and um should already be prefetched */
	/* Now prefetch forward neighbor for next dir (1+) */
	/* And gauge field */
	sp=&t_spinor[ shift->forwardNeighbor(ix,1) ];
	
	up =&(gauge_field[thissite][1]);
	     
	dslash_plus_dir0_backward_add(*sm,*um,d_oe_psi);
	     
	     
	/********************** Direction 1 ************************ */
	/* up and sp should be prefetched */
	
	/* Prefetch backwards spinor and gauge field */
	iy=shift->backwardNeighbor(ix,1);
	sm=&t_spinor[iy];
	
	um=&(gauge_field[iy][1]);
	
	
	     
	dslash_plus_dir1_forward_add(*sp,*up,d_oe_psi);
	
	
	
	/* Prefetch forwards spinor and gauge field for next dir: 2+ */
	iy=shift->forwardNeighbor(ix,2);
	sp=&t_spinor[iy];
	
	up = &(gauge_field[thissite][2]);
	
	
	
	dslash_plus_dir1_backward_add(*sm,*um,d_oe_psi);
	
	
	
	/********************** Direction 2 ************************* */
	/* Prefetch back spinor and gauge field  */
	
	iy=shift->backwardNeighbor(ix,2);
	sm=&t_spinor[iy];
	
	um=&(gauge_field[iy][2]);
	
	
	
	dslash_plus_dir2_forward_add(*sp,*up,d_oe_psi);
	
	
	/* Prefetch forward spinor and gauge field for next direction: 3+ */
	iy=shift->forwardNeighbor(ix,3);
	sp=&t_spinor[iy];
	
	up = &(gauge_field[thissite][3]);
	
	
	dslash_plus_dir2_backward_add(*sm,*um,d_oe_psi);
	
	
	/* ******************* Direction +3 **************************** */
	/* Prefetch back spinor and gauge field  3- */
	
	iy=shift->backwardNeighbor(ix,3);
	sm=&t_spinor[iy];
	
	um=&(gauge_field[iy][3]);
	
	
	dslash_plus_dir3_forward_add(*sp,*up,d_oe_psi);      
	
	cloverSiteApply(res[ thissite ], clov[thissite], psi[ thissite ]);	
	dslash_plus_dir3_backward_add_store(*sm,*um,d_oe_psi); 
	



	// Final result
	res[thissite][0][0][0] -= d_oe_psi[0][0][0];
	res[thissite][0][0][1] -= d_oe_psi[0][0][1];
	res[thissite][0][1][0] -= d_oe_psi[0][1][0];
	res[thissite][0][1][1] -= d_oe_psi[0][1][1];
	res[thissite][0][2][0] -= d_oe_psi[0][2][0];
	res[thissite][0][2][1] -= d_oe_psi[0][2][1];

	res[thissite][1][0][0] -= d_oe_psi[1][0][0];
	res[thissite][1][0][1] -= d_oe_psi[1][0][1];
	res[thissite][1][1][0] -= d_oe_psi[1][1][0];
	res[thissite][1][1][1] -= d_oe_psi[1][1][1];
	res[thissite][1][2][0] -= d_oe_psi[1][2][0];
	res[thissite][1][2][1] -= d_oe_psi[1][2][1];

	res[thissite][2][0][0] -= d_oe_psi[2][0][0];
	res[thissite][2][0][1] -= d_oe_psi[2][0][1];
	res[thissite][2][1][0] -= d_oe_psi[2][1][0];
	res[thissite][2][1][1] -= d_oe_psi[2][1][1];
	res[thissite][2][2][0] -= d_oe_psi[2][2][0];
	res[thissite][2][2][1] -= d_oe_psi[2][2][1];


	res[thissite][3][0][0] -= d_oe_psi[3][0][0];
	res[thissite][3][0][1] -= d_oe_psi[3][0][1];
	res[thissite][3][1][0] -= d_oe_psi[3][1][0];
	res[thissite][3][1][1] -= d_oe_psi[3][1][1];
	res[thissite][3][2][0] -= d_oe_psi[3][2][0];
	res[thissite][3][2][1] -= d_oe_psi[3][2][1];

      }
    }



    void ClovDPsiMinus(size_t lo, size_t hi, int id, const void *ptr )
    {
      int ix,iy,iz;                          /* ix is the current site */
      /* iy is the neighbour for prefetching */
      /* iz is the prefetch site for the 
	 next loop iteration */
      
      const CloverThreadWorkerArgs *a = (const CloverThreadWorkerArgs*)ptr;    /* Cast the void args pointer */
      ShiftTable *shift = (ShiftTable *)a->s;
      int total_vol_cb = shift->totalVolCB();
      GaugeMatrix ALIGN (*gauge_field)[4] = (GaugeMatrix (*)[4])a->u; /* Gauge field */
      FourSpinor *psi = (FourSpinor *)a->psi;                 /* Source spinor */
      FourSpinor *res = (FourSpinor *)a->res;                 /* Result spinor */
      FourSpinor *t_spinor = (FourSpinor *)a->t_spinor;
      const CloverTerm *clov = (const CloverTerm *)a->clov_oo;
      const CloverTerm *invclov = (const CloverTerm *)a->invclov_ee;
  
      GaugeMatrix *up,*um;                        /* us for multiply (PLUS/MINUS) */
      FourSpinor *sp,*sm,*rn;                   /* spinor pointers sp sm are the 
						   neighbours, rn is the result */
      
      FourSpinor d_oe_psi  __attribute__((aligned(16)));
      int thissite;
      
      
      /************************ loop over all lattice sites *************************/
      
      for (ix=lo;ix<hi;ix++) {
	thissite=shift->siteTable(ix);
      
	/* 'peel this off' to allow prefetching */
	iy=shift->forwardNeighbor(ix,0);
	sp=&psi[iy];
	up=&(gauge_field[thissite][0]);
	
	/* Prefetch back spinor and gauge field */
	iy=shift->backwardNeighbor(ix,0);
	sm=&psi[iy];

	um=&(gauge_field[iy][0]);

	
	dslash_minus_dir0_forward(*sp,*up,d_oe_psi);
	
	
	/* Prefetch spinor and gauge field for next direction (1) */
	iy=shift->forwardNeighbor(ix,1);
	sp=&psi[iy];

	up = &(gauge_field[thissite][1]);

	
	dslash_minus_dir0_backward_add(*sm,*um,d_oe_psi);    
	
	/* ******************* direction +1 ******************* */
	/* Prefetch backward spinor and gauge field */
	
	iy=shift->backwardNeighbor(ix,1);
	sm=&psi[iy];

	um=&(gauge_field[iy][1]);

	
	dslash_minus_dir1_forward_add(*sp,*up,d_oe_psi);
	
        
	/* Prefetch forward spinor and gauge field for next direction: 2 */
	iy=shift->forwardNeighbor(ix,2);
	sp=&psi[iy];

	up =&(gauge_field[thissite][2]);

	
	dslash_minus_dir1_backward_add(*sm,*um,d_oe_psi);    
	
	
	/* ******************* direction +2 **************************** */
	/* Prefetch backward spinor and gauge field */
	
	iy=shift->backwardNeighbor(ix,2); 
	sm=&psi[iy]; 

	um=&(gauge_field[iy][2]);

	
	dslash_minus_dir2_forward_add(*sp,*up,d_oe_psi);    
        
	
	/* Prefetch spinor and gauge field for next direction: 3+ */
	iy=shift->forwardNeighbor(ix,3);
	sp=&psi[iy];

	
	/* Prefetch gauge field for nex dir: 3+ */
	up = &(gauge_field[thissite][3]);

	
	dslash_minus_dir2_backward_add(*sm,*um,d_oe_psi);    
	
        
	
	/*******************  direction +3 ************************** */
	/* Prefetch backward spinor and gauge field */
	iy=shift->backwardNeighbor(ix,3);
	sm=&psi[iy];

	um=&(gauge_field[iy][3]);

	
	dslash_minus_dir3_forward_add(*sp,*up,d_oe_psi);    
	
	

	
	dslash_minus_dir3_backward_add_store(*sm,*um,d_oe_psi);    
	cloverSiteApply(t_spinor[thissite], invclov[thissite], d_oe_psi);       
      
      }

      /* Get forward neighbour of low in the x direction */
      const int low  =  total_vol_cb+lo;                 /* First site for this thread */
      const int high  = total_vol_cb+hi;                /* Last site for this thread */

      
      /************************ loop over all lattice sites *************************/
      
      for (int ix2=lo;ix2<hi;ix2++) {
	ix=ix2+total_vol_cb;
	thissite=shift->siteTable(ix);
      
	iy=shift->forwardNeighbor(ix,0);
	sp=&t_spinor[iy];
	up=&(gauge_field[thissite][0]);

	/* Prefetch back spinor and gauge field */
	iy=shift->backwardNeighbor(ix,0);
	sm=&t_spinor[iy];

	um=&(gauge_field[iy][0]);

	
	dslash_minus_dir0_forward(*sp,*up,d_oe_psi);
	
	
	/* Prefetch spinor and gauge field for next direction (1) */
	iy=shift->forwardNeighbor(ix,1);
	sp=&t_spinor[iy];

	up = &(gauge_field[thissite][1]);

	
	dslash_minus_dir0_backward_add(*sm,*um,d_oe_psi);    
	
	/* ******************* direction +1 ******************* */
	/* Prefetch backward spinor and gauge field */
	
	iy=shift->backwardNeighbor(ix,1);
	sm=&t_spinor[iy];

	um=&(gauge_field[iy][1]);

	
	dslash_minus_dir1_forward_add(*sp,*up,d_oe_psi);
	
        
	/* Prefetch forward spinor and gauge field for next direction: 2 */
	iy=shift->forwardNeighbor(ix,2);
	sp=&t_spinor[iy];

	up =&(gauge_field[thissite][2]);

	
	dslash_minus_dir1_backward_add(*sm,*um,d_oe_psi);    
	
	
	/* ******************* direction +2 **************************** */
	/* Prefetch backward spinor and gauge field */
	
	iy=shift->backwardNeighbor(ix,2); 
	sm=&t_spinor[iy]; 

	um=&(gauge_field[iy][2]);

	
	dslash_minus_dir2_forward_add(*sp,*up,d_oe_psi);    
        
	
	/* Prefetch spinor and gauge field for next direction: 3+ */
	iy=shift->forwardNeighbor(ix,3);
	sp=&t_spinor[iy];

	
	/* Prefetch gauge field for nex dir: 3+ */
	up = &(gauge_field[thissite][3]);

	
	dslash_minus_dir2_backward_add(*sm,*um,d_oe_psi);    
	
        
	
	/*******************  direction +3 ************************** */
	/* Prefetch backward spinor and gauge field */
	iy=shift->backwardNeighbor(ix,3);
	sm=&t_spinor[iy];

	um=&(gauge_field[iy][3]);

	
	dslash_minus_dir3_forward_add(*sp,*up,d_oe_psi);    
	

	
	dslash_minus_dir3_backward_add_store(*sm,*um,d_oe_psi);    
	// res = A_oo psi 
	cloverSiteApply(res[ thissite ], clov[thissite], psi[ thissite ]);

	// Subtract d_oe_psi from a_oo
	res[thissite][0][0][0] -= d_oe_psi[0][0][0];
	res[thissite][0][0][1] -= d_oe_psi[0][0][1];
	res[thissite][0][1][0] -= d_oe_psi[0][1][0];
	res[thissite][0][1][1] -= d_oe_psi[0][1][1];
	res[thissite][0][2][0] -= d_oe_psi[0][2][0];
	res[thissite][0][2][1] -= d_oe_psi[0][2][1];

	res[thissite][1][0][0] -= d_oe_psi[1][0][0];
	res[thissite][1][0][1] -= d_oe_psi[1][0][1];
	res[thissite][1][1][0] -= d_oe_psi[1][1][0];
	res[thissite][1][1][1] -= d_oe_psi[1][1][1];
	res[thissite][1][2][0] -= d_oe_psi[1][2][0];
	res[thissite][1][2][1] -= d_oe_psi[1][2][1];

	res[thissite][2][0][0] -= d_oe_psi[2][0][0];
	res[thissite][2][0][1] -= d_oe_psi[2][0][1];
	res[thissite][2][1][0] -= d_oe_psi[2][1][0];
	res[thissite][2][1][1] -= d_oe_psi[2][1][1];
	res[thissite][2][2][0] -= d_oe_psi[2][2][0];
	res[thissite][2][2][1] -= d_oe_psi[2][2][1];


	res[thissite][3][0][0] -= d_oe_psi[3][0][0];
	res[thissite][3][0][1] -= d_oe_psi[3][0][1];
	res[thissite][3][1][0] -= d_oe_psi[3][1][0];
	res[thissite][3][1][1] -= d_oe_psi[3][1][1];
	res[thissite][3][2][0] -= d_oe_psi[3][2][0];
	res[thissite][3][2][1] -= d_oe_psi[3][2][1];
      
      }
    }

  }

} // Namespace
