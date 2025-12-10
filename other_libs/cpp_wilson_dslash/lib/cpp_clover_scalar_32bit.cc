#include <cstdlib>
#include <iostream>

#include <cpp_dslash_scalar.h>
#include <shift_table_scalar.h>
#include <dispatch_scalar.h>

#include <cpp_clover_scalar.h>
#include <cpp_dslash_scalar_32bit.h>
#include <cpp_clover_site_apply_32bit.h>
#include "allocate.h"
using namespace CPlusPlusWilsonDslash::DslashScalar32Bit;
using namespace CPlusPlusWilsonDslash::Dslash32BitTypes;
using namespace CPlusPlusWilsonDslash::Cache;

namespace CPlusPlusClover {
 

  /* Constructor */
  CloverSchur4D<float>::CloverSchur4D(const int latt_size[],      
				      void (*getSiteCoords)(int coord[], int node, int linearsite),
				      int (*getLinearSiteIndex)(const int coord[]),
				      int (*nodeNum)(const int coord[])
		) 
  {
    s = new CPlusPlusWilsonDslash::ShiftTable(latt_size,
					      getSiteCoords,
					      getLinearSiteIndex,
					      nodeNum);

    xt_spinor = (FourSpinor *)CPlusPlusWilsonDslash::alloc(2*s->totalVolCB()*sizeof(FourSpinor)
				    +Cache::CacheLineSize);
    if( xt_spinor == (FourSpinor *)NULL ) { 
      std::cerr << "Unable to allocate temporary" << std::endl;
      exit(1);
    }
    unsigned long pad = 0;
    if ( (unsigned long)xt_spinor % Cache::CacheLineSize != 0 ) { 
      pad = Cache::CacheLineSize - (unsigned long)xt_spinor % Cache::CacheLineSize;
    }
    t_spinor = (FourSpinor *)((unsigned char *)xt_spinor + pad);
  }
      

  CloverSchur4D<float>::~CloverSchur4D() { delete s; CPlusPlusWilsonDslash::dealloc(xt_spinor); }



  namespace CloverScalar32Bit {

   
    void DClovPsiPlus(size_t lo, size_t hi, int id, const void *ptr);
    void DClovPsiMinus(size_t lo, size_t hi, int id, const void *ptr);

  }

  using namespace CloverScalar32Bit;

  // The operator 
  void CloverSchur4D<float>::operator() (float* res, 
					 const float* psi, 
					 const float *u, /* Gauge field suitably packed */
					 
					 const float *clov_oo,
					 const float *invclov_ee,
					 int isign)

  {
    if (isign == 1) {
      

      CPlusPlusClover::dispatchToThreads((void (*)(size_t, size_t, int, const void *))&DClovPsiPlus, 
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

      CPlusPlusClover::dispatchToThreads((void (*)(size_t, size_t, int, const void *))&DClovPsiMinus, 
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

  using namespace DslashScalar32Bit;
  using namespace CPlusPlusClover32Bit;
  namespace CloverScalar32Bit {

    
    void DClovPsiPlus(size_t lo, size_t hi, int id, const void *ptr)
    {
      const CloverThreadWorkerArgs *a  
	= (const CloverThreadWorkerArgs*)ptr;  /* Cast the (void *) to an (ThreadWorkerArgs*) */

      int ix1;                              /* Index of current site */
      int iy1,iy2;                          /* Index of neighbour of current site (iy1) 
					       and current site+1 (iy2) */
      
      int iz1;                              /* Index of next site in loop */
      
      GaugeMatrix (*gauge_field)[4]  =  (GaugeMatrix (*)[4])a->u; /* Packed Gauge fields */
      FourSpinor *psi  =  (FourSpinor *)a->psi;           /* Source spinor */
      FourSpinor *res  =  (FourSpinor *)a->res;           /* Result spinor */
      FourSpinor *t_spinor = (FourSpinor *)a->t_spinor;   /* Temporary spinor */
      
      ShiftTable *s = (ShiftTable *)(a->s); 
      int total_vol_cb = s->totalVolCB();   
      
      CloverTerm* invclov = (CloverTerm *)a->invclov_ee;
      CloverTerm* clov = (CloverTerm *)a->clov_oo;
      
      /* Pointers to the neighboring u-s */
      GaugeMatrix *up1 ALIGN;                  /* U[ x  ] */
      GaugeMatrix *um1 ALIGN;                  /* U[ x - mu ] */
      
    /* 4 - Spinor pointers */
      FourSpinor *sp1 ALIGN;
      FourSpinor *sm1 ALIGN;
      FourSpinor *sn1 ALIGN;
      
      /* Half Vectors */
      HalfSpinor r12_1 ALIGN; /* Site 1 upper */
      HalfSpinor r34_1 ALIGN; /* Site 1 lower */
      
      FourSpinor d_oe_psi ALIGN;

      /* note that we want the spinors from the checkerboard opposite the one we are writing to */
      /* We are doing things in bunches of two sites */
      

 
    

      // First loop: Target is the even subset (from the odd)
      // So cb = 0 for even
      // since ix1 = cb*total_vol_cb + i and cb = 0 => ix1 = i
      // can use ix1 directly as a loop counter
      for(int ix1 = lo; ix1 < hi; ix1++) { 

	// This is an even site
	int thissite = s->siteTable( ix1 );
	
	// These are odd sites
	int fsite = s->forwardNeighbor(ix1,0);
	int bsite = s->backwardNeighbor(ix1,0);
	
	/******************************* direction +0 *********************************/
	
	/* ...(1-isign*gamma(0)) psi(x + \hat{0}) */
	
	/* Prefetch the backward neighbours for the following 1 + isign gamma(0) case */
	sp1 = &psi[ fsite  ];
	up1 = &(gauge_field[thissite][0]);
	dslash_plus_dir0_forward(*sp1, *up1, r12_1, r34_1);
	
	/* Now prefetch for the  1 + \gamma_0 U^\dagger case */
	sm1 = &psi[ bsite ];
	um1 = &(gauge_field[bsite][0]);
	dslash_plus_dir0_backward_add(*sm1, *um1, r12_1, r34_1);
      
	/* Prefetch gauge field for next direction */
	fsite = s->forwardNeighbor(ix1,1);
	bsite = s->backwardNeighbor(ix1,1);
	
	up1 = &(gauge_field[thissite][1]);
	sp1 = &psi[ fsite ];
	dslash_plus_dir1_forward_add(*sp1, *up1, r12_1, r34_1);
	
	/* Prefetch spinors for the -1 direction */
	sm1 = &psi[bsite];
	um1 = &(gauge_field[bsite][1]);
	dslash_plus_dir1_backward_add(*sm1, *um1, r12_1, r34_1);
	
	/* Prefetch forward neighbour for direction 2+ */
	fsite = s->forwardNeighbor(ix1,2);
	bsite = s->backwardNeighbor(ix1,2);
	
	up1 = &(gauge_field[thissite][2]);
	sp1 = &psi[fsite];
	dslash_plus_dir2_forward_add(*sp1, *up1, r12_1, r34_1);
	
	/* Prefetch sm1 & sm2 for -ve direction */
	sm1 = &psi[bsite];
	um1 = &(gauge_field[bsite][2]);
	dslash_plus_dir2_backward_add(*sm1, *um1, r12_1, r34_1);
	
      
	/* Prefetch spinors for direction 3+ */
	fsite = s->forwardNeighbor(ix1,3);
	bsite = s->backwardNeighbor(ix1,3);
	
	sp1 = &psi[ fsite ];
	up1 = &(gauge_field[thissite][3]);
	dslash_plus_dir3_forward_add(*sp1, *up1, r12_1, r34_1);
	
      
	sm1 = &psi[bsite]; 
	um1 = &(gauge_field[bsite][3]); 
	dslash_plus_dir3_backward_add_store(*sm1, *um1, r12_1, r34_1, d_oe_psi);
	// Apply Clover Inverse	  
	cloverSiteApply(t_spinor[thissite], invclov[thissite], d_oe_psi);


      }

      // Now the opposite checkerboard...
      for(int ix2= lo; ix2 < hi; ix2++) { 
	int ix1 = total_vol_cb + ix2;
	
	int thissite = s->siteTable( ix1 );

	int fsite = s->forwardNeighbor(ix1,0);
	int bsite = s->backwardNeighbor(ix1,0);
	
	/******************************* direction +0 *********************************/
	
	/* ...(1-isign*gamma(0)) psi(x + \hat{0}) */
	
	/* Prefetch the backward neighbours for the following 1 + isign gamma(0) case */
	sp1 = &t_spinor[ fsite  ];
	up1 = &(gauge_field[thissite][0]);
	dslash_plus_dir0_forward(*sp1, *up1, r12_1, r34_1);
      
	/* Now prefetch for the  1 + \gamma_0 U^\dagger case */
	sm1 = &t_spinor[ bsite ];
	um1 = &(gauge_field[bsite][0]);
	dslash_plus_dir0_backward_add(*sm1, *um1, r12_1, r34_1);
      
	/* Prefetch gauge field for next direction */
	fsite = s->forwardNeighbor(ix1,1);
	bsite = s->backwardNeighbor(ix1,1);
	
	up1 = &(gauge_field[thissite][1]);
	sp1 = &t_spinor[ fsite ];
	dslash_plus_dir1_forward_add(*sp1, *up1, r12_1, r34_1);
	
	/* Prefetch spinors for the -1 direction */
	sm1 = &t_spinor[bsite];
	um1 = &(gauge_field[bsite][1]);
	dslash_plus_dir1_backward_add(*sm1, *um1, r12_1, r34_1);
	
	/* Prefetch forward neighbour for direction 2+ */
	fsite = s->forwardNeighbor(ix1,2);
	bsite = s->backwardNeighbor(ix1,2);
	
	up1 = &(gauge_field[thissite][2]);
	sp1 = &t_spinor[fsite];
	dslash_plus_dir2_forward_add(*sp1, *up1, r12_1, r34_1);
	
	/* Prefetch sm1 & sm2 for -ve direction */
	sm1 = &t_spinor[bsite];
	um1 = &(gauge_field[bsite][2]);
	dslash_plus_dir2_backward_add(*sm1, *um1, r12_1, r34_1);
	
      
	/* Prefetch spinors for direction 3+ */
	fsite = s->forwardNeighbor(ix1,3);
	bsite = s->backwardNeighbor(ix1,3);
	
	sp1 = &t_spinor[ fsite ];
	up1 = &(gauge_field[thissite][3]);
	dslash_plus_dir3_forward_add(*sp1, *up1, r12_1, r34_1);
	
      
	sm1 = &t_spinor[bsite]; 
	um1 = &(gauge_field[bsite][3]); 
	//	sn1 = &res[thissite];  /*we always walk across the result lexicographically */
	
	// d_oe_psi now holds:  D_oe A^{-1}_ee D_eo
	dslash_plus_dir3_backward_add_store(*sm1, *um1, r12_1, r34_1, d_oe_psi);

	// res = A_oo psi 
	cloverSiteApply(res[ thissite ], clov[thissite], psi[ thissite ]);


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
    
    void DClovPsiMinus(size_t lo, size_t hi, int id, const void *ptr)
    {
      const CloverThreadWorkerArgs *a  = (const CloverThreadWorkerArgs*)ptr;   /* Downcast to args */
      int ix1,iy1,iy2,iz1;                   /* Coordinates ix1 - current
					      iy1 - index of neighbour
					      iy1 - index of neighbour of ix1+1 
					      iz1 - index of first of next pair (ix+2) */
      //      const int cb = a->cb;
    
    //  const int low  =  icolor_start[cb]+lo;                 /* First site for this thread */
    // const int high  = icolor_start[cb]+ hi;                /* Last site for this thread */
      
      GaugeMatrix (*gauge_field)[4]  = (GaugeMatrix (*)[4]) a->u;  /* Gauge field */
      
      FourSpinor *psi  =  (FourSpinor *)a->psi;            /* Source 4-spinor */
      FourSpinor *res  =  (FourSpinor *)a->res;            /* Result 4-spinor */
      FourSpinor *t_spinor = (FourSpinor *)a->t_spinor;
      const CloverTerm *clov = (const CloverTerm *)a->clov_oo;
      const CloverTerm *invclov = (const CloverTerm *)a->invclov_ee;
      
      ShiftTable *s = (ShiftTable *)a->s;
      int total_vol_cb = s->totalVolCB();
      
      GaugeMatrix *up1 ALIGN;                  /* U[ ix ] */
      GaugeMatrix *um1 ALIGN;                  /* U[ ix - mu ] */
      
      FourSpinor *sp1 ALIGN;                 /* 4 spinor psi[ix1+mu] */
      FourSpinor *sm1 ALIGN;                 /* 4 spinor psi[ix1-mu] */
      FourSpinor *sn1 ALIGN;                 /* 4 spinor result */
      
      HalfSpinor r12_1;                         /* site 1 halfspinor top half */
      HalfSpinor r34_1;                         /* site 1 halfspinor bottom half */
      FourSpinor d_oe_psi ALIGN;

#if 0
      /* Get forward neighbour of low in the x direction */
      const int low  =  cb*total_vol_cb+lo;                 /* First site for this thread */
      const int high  = cb*total_vol_cb+hi;                /* Last site for this thread */
     
#endif 
      /************************ loop over all lattice sites *************************/
      
      for (ix1 = lo;ix1<hi;ix1++) {
	
	// ix1 loops through even sites (cb=0)
	int thissite = s->siteTable( ix1 );

	// Odd Sites 
	int fsite = s->forwardNeighbor(ix1,0);
	
	sp1 = &psi[ fsite  ];
	up1 = &(gauge_field[thissite][0]);
	dslash_minus_dir0_forward(*sp1, *up1, r12_1, r34_1);
      
	/* Now prefetch for the  1 + \gamma_0 U^\dagger case */
	iy1 =  s->backwardNeighbor(ix1,0);
	sm1 = &psi[ iy1 ];
	um1 = &(gauge_field[iy1][0]);
	dslash_minus_dir0_backward_add(*sm1, *um1, r12_1, r34_1);
      
	/* Prefetch gauge field for next direction */
	fsite = s->forwardNeighbor(ix1,1);
	up1 = &(gauge_field[thissite][1]);
	sp1 = &psi[ fsite ];
	dslash_minus_dir1_forward_add(*sp1, *up1, r12_1, r34_1);
	
	/* Prefetch spinors for the -1 direction */
	iy1 = s->backwardNeighbor(ix1,1);
	sm1 = &psi[iy1];
	um1 = &(gauge_field[iy1][1]);
	dslash_minus_dir1_backward_add(*sm1, *um1, r12_1, r34_1);
	
	/* Prefetch forward neighbour for direction 2+ */
	fsite = s->forwardNeighbor(ix1,2);
	up1 = &(gauge_field[thissite][2]);
	sp1 = &psi[fsite];
	dslash_minus_dir2_forward_add(*sp1, *up1, r12_1, r34_1);
	
	/* Prefetch sm1 & sm2 for -ve direction */
	iy1 = s->backwardNeighbor(ix1,2);
	sm1 = &psi[iy1];
	um1 = &(gauge_field[iy1][2]);
	dslash_minus_dir2_backward_add(*sm1, *um1, r12_1, r34_1);
	
	
	/* Prefetch spinors for direction 3+ */
	fsite = s->forwardNeighbor(ix1,3);
	
	up1 = &(gauge_field[thissite][3]);
	sp1 = &psi[ fsite  ];
	dslash_minus_dir3_forward_add(*sp1, *up1, r12_1, r34_1);
	
	
	iy1 = s->backwardNeighbor(ix1,3);
	sm1 = &psi[iy1]; 
	um1 = &(gauge_field[iy1][3]); 

	dslash_minus_dir3_backward_add_store(*sm1, *um1, r12_1, r34_1, d_oe_psi);
	cloverSiteApply(t_spinor[thissite], invclov[thissite], d_oe_psi);       
      }

      // These are now even sites
      for (int ix2 = lo;ix2<hi;ix2++) {

	// Odd sites
	ix1 = ix2 + total_vol_cb;
	int thissite = s->siteTable( ix1 );
       
	// Even site
	int fsite = s->forwardNeighbor(ix1,0);
	

	sp1 = &t_spinor[ fsite  ];
	up1 = &(gauge_field[thissite][0]);
	dslash_minus_dir0_forward(*sp1, *up1, r12_1, r34_1);
      
	/* Now prefetch for the  1 + \gamma_0 U^\dagger case */
	iy1 =  s->backwardNeighbor(ix1,0);
	sm1 = &t_spinor[ iy1 ];
	um1 = &(gauge_field[iy1][0]);
	dslash_minus_dir0_backward_add(*sm1, *um1, r12_1, r34_1);
      
	/* Prefetch gauge field for next direction */
	fsite = s->forwardNeighbor(ix1,1);
	up1 = &(gauge_field[thissite][1]);
	sp1 = &t_spinor[ fsite ];
	dslash_minus_dir1_forward_add(*sp1, *up1, r12_1, r34_1);
	
	/* Prefetch spinors for the -1 direction */
	iy1 = s->backwardNeighbor(ix1,1);
	sm1 = &t_spinor[iy1];
	um1 = &(gauge_field[iy1][1]);
	dslash_minus_dir1_backward_add(*sm1, *um1, r12_1, r34_1);
	
	/* Prefetch forward neighbour for direction 2+ */
	fsite = s->forwardNeighbor(ix1,2);
	up1 = &(gauge_field[thissite][2]);
	sp1 = &t_spinor[fsite];
	dslash_minus_dir2_forward_add(*sp1, *up1, r12_1, r34_1);
	
	/* Prefetch sm1 & sm2 for -ve direction */
	iy1 = s->backwardNeighbor(ix1,2);
	sm1 = &t_spinor[iy1];
	um1 = &(gauge_field[iy1][2]);
	dslash_minus_dir2_backward_add(*sm1, *um1, r12_1, r34_1);
	
	
	/* Prefetch spinors for direction 3+ */
	fsite = s->forwardNeighbor(ix1,3);
	
	up1 = &(gauge_field[thissite][3]);
	sp1 = &t_spinor[ fsite  ];
	dslash_minus_dir3_forward_add(*sp1, *up1, r12_1, r34_1);
	
	
	iy1 = s->backwardNeighbor(ix1,3);
	sm1 = &t_spinor[iy1]; 
	um1 = &(gauge_field[iy1][3]); 

	dslash_minus_dir3_backward_add_store(*sm1, *um1, r12_1, r34_1, d_oe_psi);


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
}

