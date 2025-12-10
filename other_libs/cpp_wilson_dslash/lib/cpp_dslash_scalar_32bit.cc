#include <cpp_dslash_scalar.h>
#include <shift_table_scalar.h>
#include <dispatch_scalar.h>

#include <cpp_dslash_scalar_32bit.h>

using namespace CPlusPlusWilsonDslash::DslashScalar32Bit;
using namespace CPlusPlusWilsonDslash::Dslash32BitTypes;

namespace CPlusPlusWilsonDslash {
 

  /* Constructor */
  Dslash<float>::Dslash(const int latt_size[],      
		void (*getSiteCoords)(int coord[], int node, int linearsite),
		int (*getLinearSiteIndex)(const int coord[]),
		int (*nodeNum)(const int coord[])
		) 
  {
    s = new ShiftTable(latt_size,
		       getSiteCoords,
		       getLinearSiteIndex,
		       nodeNum);

  }
      
  Dslash<float>::~Dslash() { delete s; }

  //  int Dslash<float>::getPathSite(int site) const
  //  { 
  //    return s->getPathSite(site);
  //  }
 
  // The operator 
  void Dslash<float>::operator() (float* res, 
				  float* psi, 
				  float *u, /* Gauge field suitably packed */
				  int isign,
				  int cb) 
  {
    if (isign == 1) {  

      CPlusPlusWilsonDslash::dispatchToThreads((void (*)(size_t, size_t, int, const void *))&DPsiPlus, 
					       (void*)psi,
					       (void*)res,
					       (void *)u,
					       (void *)s,
					       1-cb,
					       s->totalVolCB());


    }

    if( isign == -1) {

      CPlusPlusWilsonDslash::dispatchToThreads((void (*)(size_t, size_t, int, const void *))&DPsiMinus, 
					       (void *)psi,
					       (void *)res,
					       (void *)u,
					       (void *)s,
					       1-cb,
					       (int)s->totalVolCB());
      
    }    
  }


  namespace DslashScalar32Bit {
    
    void DPsiPlus(size_t lo, size_t hi, int id, const void *ptr)
    {
      const ThreadWorkerArgs *a  
	= (const ThreadWorkerArgs*)ptr;  /* Cast the (void *) to an (ThreadWorkerArgs*) */

      int ix1;                              /* Index of current site */
      int iy1,iy2;                          /* Index of neighbour of current site (iy1) 
					     and current site+1 (iy2) */
    
      int iz1;                              /* Index of next site in loop */
      
      GaugeMatrix (*gauge_field)[4]  =  (GaugeMatrix (*)[4])a->u; /* Packed Gauge fields */
      FourSpinor *psi  =  (FourSpinor *)a->psi;           /* Source spinor */
      FourSpinor *res  =  (FourSpinor *)a->res;           /* Result spinor */
    
      const int cb = a->cb;
      ShiftTable *s = (ShiftTable *)(a->s); 
      int total_vol_cb = s->totalVolCB();   
    
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
      
      
      /* note that we want the spinors from the checkerboard opposite the one we are writing to */
      /* We are doing things in bunches of two sites */
      
      /* Get forward neighbour of low in the x direction */
      /* Get forward neighbour of low in the x direction */
      const int low  =  cb*total_vol_cb+lo;                 /* First site for this thread */
      const int high  = cb*total_vol_cb+hi;                /* Last site for this thread */
    

      for (int ix1 = low; ix1 < high; ix1++)  {
	//int thissite = s->getPathSite(ix1);

	int thissite = s->siteTable( ix1 );
	int fsite = s->forwardNeighbor(ix1,0);
	int bsite = s->backwardNeighbor(ix1,0);
	
	/******************************* direction +0 *********************************/
	
	/* ...(1-isign*gamma(0)) psi(x + \hat{0}) */
	
	/* Prefetch the backward neighbours for the following 1 + isign gamma(0) case */
	sp1 = &psi[ fsite  ];
	up1 = &(gauge_field[ix1][0]);
	dslash_plus_dir0_forward(*sp1, *up1, r12_1, r34_1);
      
	/* Now prefetch for the  1 + \gamma_0 U^\dagger case */
	sm1 = &psi[ bsite ];
	um1 = &(gauge_field[bsite][0]);
	dslash_plus_dir0_backward_add(*sm1, *um1, r12_1, r34_1);
      
	/* Prefetch gauge field for next direction */
	fsite = s->forwardNeighbor(ix1,1);
	bsite = s->backwardNeighbor(ix1,1);
	
	up1 = &(gauge_field[ix1][1]);
	sp1 = &psi[ fsite ];
	dslash_plus_dir1_forward_add(*sp1, *up1, r12_1, r34_1);
	
	/* Prefetch spinors for the -1 direction */
	sm1 = &psi[bsite];
	um1 = &(gauge_field[bsite][1]);
	dslash_plus_dir1_backward_add(*sm1, *um1, r12_1, r34_1);
	
	/* Prefetch forward neighbour for direction 2+ */
	fsite = s->forwardNeighbor(ix1,2);
	bsite = s->backwardNeighbor(ix1,2);
	
	up1 = &(gauge_field[ix1][2]);
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
	up1 = &(gauge_field[ix1][3]);
	dslash_plus_dir3_forward_add(*sp1, *up1, r12_1, r34_1);
	
      
	sm1 = &psi[bsite]; 
	um1 = &(gauge_field[bsite][3]); 
	sn1 = &res[ix1];  /*we always walk across the result lexicographically */
	dslash_plus_dir3_backward_add_store(*sm1, *um1, r12_1, r34_1, *sn1);
	
      }
    }
    
    void DPsiMinus(size_t lo, size_t hi, int id, const void *ptr)
    {
      const ThreadWorkerArgs *a  = (const ThreadWorkerArgs*)ptr;   /* Downcast to args */
      int ix1,iy1,iy2,iz1;                   /* Coordinates ix1 - current
					      iy1 - index of neighbour
					      iy1 - index of neighbour of ix1+1 
					      iz1 - index of first of next pair (ix+2) */
      const int cb = a->cb;
    
    //  const int low  =  icolor_start[cb]+lo;                 /* First site for this thread */
    // const int high  = icolor_start[cb]+ hi;                /* Last site for this thread */
      
      GaugeMatrix (*gauge_field)[4]  = (GaugeMatrix (*)[4]) a->u;  /* Gauge field */
      
      FourSpinor *psi  =  (FourSpinor *)a->psi;            /* Source 4-spinor */
      FourSpinor *res  =  (FourSpinor *)a->res;            /* Result 4-spinor */
      ShiftTable *s = (ShiftTable *)a->s;
      int total_vol_cb = s->totalVolCB();
      
      GaugeMatrix *up1 ALIGN;                  /* U[ ix ] */
      GaugeMatrix *um1 ALIGN;                  /* U[ ix - mu ] */
      
      FourSpinor *sp1 ALIGN;                 /* 4 spinor psi[ix1+mu] */
      FourSpinor *sm1 ALIGN;                 /* 4 spinor psi[ix1-mu] */
      FourSpinor *sn1 ALIGN;                 /* 4 spinor result */
      
      HalfSpinor r12_1;                         /* site 1 halfspinor top half */
      HalfSpinor r34_1;                         /* site 1 halfspinor bottom half */
      
      /* Get forward neighbour of low in the x direction */
      const int low  =  cb*total_vol_cb+lo;                 /* First site for this thread */
      const int high  = cb*total_vol_cb+hi;                /* Last site for this thread */
      
      /************************ loop over all lattice sites *************************/
      for (int ix1 = low ; ix1 < high; ix1++) {
	//	ix1 = s->getPathSite(cb, site_iter);

	//	int thissite = s->getPathSite( ix1 );
	int thissite = s->siteTable(ix1);
	int fsite = s->forwardNeighbor(ix1,0);
	
	sp1 = &psi[ fsite  ];
	up1 = &(gauge_field[ix1][0]);
	dslash_minus_dir0_forward(*sp1, *up1, r12_1, r34_1);
      
	/* Now prefetch for the  1 + \gamma_0 U^\dagger case */
	iy1 =  s->backwardNeighbor(ix1,0);
	sm1 = &psi[ iy1 ];
	um1 = &(gauge_field[iy1][0]);
	dslash_minus_dir0_backward_add(*sm1, *um1, r12_1, r34_1);
      
	/* Prefetch gauge field for next direction */
	fsite = s->forwardNeighbor(ix1,1);
	up1 = &(gauge_field[ix1][1]);
	sp1 = &psi[ fsite ];
	dslash_minus_dir1_forward_add(*sp1, *up1, r12_1, r34_1);
	
	/* Prefetch spinors for the -1 direction */
	iy1 = s->backwardNeighbor(ix1,1);
	sm1 = &psi[iy1];
	um1 = &(gauge_field[iy1][1]);
	dslash_minus_dir1_backward_add(*sm1, *um1, r12_1, r34_1);
	
	/* Prefetch forward neighbour for direction 2+ */
	fsite = s->forwardNeighbor(ix1,2);
	up1 = &(gauge_field[ix1][2]);
	sp1 = &psi[fsite];
	dslash_minus_dir2_forward_add(*sp1, *up1, r12_1, r34_1);
	
	/* Prefetch sm1 & sm2 for -ve direction */
	iy1 = s->backwardNeighbor(ix1,2);
	sm1 = &psi[iy1];
	um1 = &(gauge_field[iy1][2]);
	dslash_minus_dir2_backward_add(*sm1, *um1, r12_1, r34_1);
	
	
	/* Prefetch spinors for direction 3+ */
	fsite = s->forwardNeighbor(ix1,3);
	
	up1 = &(gauge_field[ix1][3]);
	sp1 = &psi[ fsite  ];
	dslash_minus_dir3_forward_add(*sp1, *up1, r12_1, r34_1);
	
	
	iy1 = s->backwardNeighbor(ix1,3);
	sm1 = &psi[iy1]; 
	um1 = &(gauge_field[iy1][3]); 
	sn1 = &res[ix1];  /*we always walk across the result lexicographically */
	dslash_minus_dir3_backward_add_store(*sm1, *um1, r12_1, r34_1, *sn1);
	
      }
      
    }
    
  }

} // Namespace
