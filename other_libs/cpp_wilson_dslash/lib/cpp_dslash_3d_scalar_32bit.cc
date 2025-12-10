#include <cpp_dslash_scalar.h>
#include <shift_table_3d_scalar.h>
#include <dispatch_scalar.h>
#include <cpp_dslash_scalar_32bit.h>

using namespace CPlusPlusWilsonDslash::DslashScalar32Bit;
using namespace CPlusPlusWilsonDslash::Dslash32BitTypes;

namespace CPlusPlusWilsonDslash {
   /* Constructor */
  Dslash3D<float>::Dslash3D(const int latt_size[],      
			    void (*getSiteCoords)(int coord[], int node, int linearsite),
			    int (*getLinearSiteIndex)(const int coord[]),
			    int (*nodeNum)(const int coord[])
			    ) 
  {
    s = new ShiftTable3D(latt_size,
			 getSiteCoords,
			 getLinearSiteIndex,
			 nodeNum);

  }
      
  Dslash3D<float>::~Dslash3D() { delete s; }

  // The operator 
  void Dslash3D<float>::operator() (float* res, 
				    float* psi, 
				    float *u, /* Gauge field suitably packed */
				    int isign,
				    int cb) 
  {
    if (isign == 1) {  
      
      CPlusPlusWilsonDslash::dispatchToThreads((void (*)(size_t, size_t, int, const void *))&DPsiPlus3D, 
					       (void*)psi,
					       (void*)res,
					       (void *)u,
					       (void *)s,
					       1-cb,
					       s->totalVolCB());


    }

    if( isign == -1) {

      CPlusPlusWilsonDslash::dispatchToThreads((void (*)(size_t, size_t, int, const void *))&DPsiMinus3D, 
					       (void *)psi,
					       (void *)res,
					       (void *)u,
					       (void *)s,
					       1-cb,
					       (int)s->totalVolCB());
      
    }    
  }


  namespace DslashScalar32Bit {
    
    void DPsiPlus3D(size_t lo, size_t hi, int id, const void *ptr)
    {
      const ThreadWorkerArgs *a  = (const ThreadWorkerArgs*)ptr;  /* Cast the (void *) to an (ThreadWorkerArgs*) */
      int ix1;                              /* Index of current site */
      int iy1,iy2;                          /* Index of neighbour of current site (iy1) 
					       and current site+1 (iy2) */
      
      int iz1;                              /* Index of next site in loop */
      
      GaugeMatrix (*gauge_field)[4]  =  (GaugeMatrix (*)[4])a->u; /* Packed Gauge fields */
      FourSpinor *psi  =  (FourSpinor *)a->psi;           /* Source spinor */
      FourSpinor *res  =  (FourSpinor *)a->res;           /* Result spinor */
      
      const int cb = a->cb;
      ShiftTable3D *s = (ShiftTable3D *)(a->s); 
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
      const int low  =  cb*total_vol_cb+lo;                 /* First site for this thread */
      const int high  = cb*total_vol_cb+hi;                /* Last site for this thread */
      
      for (ix1 = low; ix1< high; ix1++) {
	
	int thissite = s->siteTable(ix1);
	int fsite = s->forwardNeighbor(ix1,0);
	int bsite = s->backwardNeighbor(ix1,0);

	/* Dir = 0 */
	sp1 = &psi[ fsite ];
	up1 = &(gauge_field[thissite][0]);
	dslash_plus_dir0_forward(*sp1, *up1, r12_1, r34_1);

	sm1 = &psi[bsite];	
	um1 = &(gauge_field[bsite][0]);       
	dslash_plus_dir0_backward_add(*sm1, *um1, r12_1, r34_1);
	
	/* Dir = 1 */
	fsite = s->forwardNeighbor(ix1,1);
	bsite = s->backwardNeighbor(ix1,1);
	
	sp1 = &psi[ fsite ];
	up1 = &(gauge_field[thissite][1]);
	dslash_plus_dir1_forward_add(*sp1, *up1, r12_1, r34_1);

	sm1 = &psi[bsite];	
	um1 = &(gauge_field[bsite][1]);
	dslash_plus_dir1_backward_add(*sm1, *um1, r12_1, r34_1);
	
	
	/* Dir = 2 */
	fsite = s->forwardNeighbor(ix1,2);
	bsite = s->backwardNeighbor(ix1,2);
	
	sp1 = &psi[fsite ];
	up1 = &(gauge_field[thissite][2]);
	dslash_plus_dir2_forward_add(*sp1, *up1, r12_1, r34_1);
	
	sm1 = &psi[bsite];
	um1 = &(gauge_field[bsite][2]);
	sn1 = &res[thissite];
	dslash_plus_dir2_backward_add_store(*sm1, *um1, r12_1, r34_1, *sn1);
	
      }
    }
      
    void DPsiMinus3D(size_t lo, size_t hi, int id, const void *ptr)
    {
      const ThreadWorkerArgs *a  = (const ThreadWorkerArgs*)ptr;   /* Downcast to args */
      int ix1,iy1,iy2,iz1;                   /* Coordinates ix1 - current
						iy1 - index of neighbour
						iy1 - index of neighbour of ix1+1 
						iz1 - index of first of next pair (ix+2) */
      const int cb = a->cb;
      ShiftTable3D *s = (ShiftTable3D *)(a->s); 
      int total_vol_cb = s->totalVolCB();         
      GaugeMatrix (*gauge_field)[4]  = (GaugeMatrix (*)[4])a->u;  /* Gauge field */
      
      FourSpinor *psi  = (FourSpinor *)a->psi;            /* Source 4-spinor */
      FourSpinor *res  = (FourSpinor *)a->res;            /* Result 4-spinor */
      
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
      
      for (ix1 = low; ix1<high; ix1++)  {
	
	int thissite = s->siteTable( ix1 );
	int x_plus_mu = s->forwardNeighbor(ix1,0);
	/* Dir = 0 */
	sp1 = &psi[ s->forwardNeighbor(ix1,0)  ];
	up1 = &(gauge_field[thissite][0]);
	dslash_minus_dir0_forward(*sp1, *up1, r12_1, r34_1);
	
	iy1 =  s->backwardNeighbor(ix1,0);
	um1 = &(gauge_field[iy1][0]);       
	sm1 = &psi[ iy1 ];
	dslash_minus_dir0_backward_add(*sm1, *um1, r12_1, r34_1);
	
	/* Dir = 1 */
	sp1 = &psi[ s->forwardNeighbor(ix1,1) ];
	up1 = &(gauge_field[thissite][1]);
	dslash_minus_dir1_forward_add(*sp1, *up1, r12_1, r34_1);
	
	iy1 = s->backwardNeighbor(ix1,1);
	um1 = &(gauge_field[iy1][1]);
	sm1 = &psi[iy1];
	dslash_minus_dir1_backward_add(*sm1, *um1, r12_1, r34_1);
	
	
	/* Dir = 2 */
	sp1 = &psi[s->forwardNeighbor(ix1,2)];
	up1 = &(gauge_field[thissite][2]);
	dslash_minus_dir2_forward_add(*sp1, *up1, r12_1, r34_1);
	
	iy1 = s->backwardNeighbor(ix1,2);
	sm1 = &psi[iy1];
	um1 = &(gauge_field[iy1][2]);
	sn1 = &res[thissite];
	dslash_minus_dir2_backward_add_store(*sm1, *um1, r12_1, r34_1, *sn1);
	
      }
    }
  }

} // Namespace
