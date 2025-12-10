#include <cpp_dslash_scalar.h>
#include <shift_table_3d_scalar.h>
#include <dispatch_scalar.h>
#include <cpp_dslash_scalar_64bit.h>

using namespace CPlusPlusWilsonDslash::DslashScalar64Bit;
using namespace CPlusPlusWilsonDslash::Dslash64BitTypes;

namespace CPlusPlusWilsonDslash {
 

  /* Constructor */
  Dslash3D<double>::Dslash3D(const int latt_size[],      
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
      
  Dslash3D<double>::~Dslash3D() { delete s; }

  // The operator 
  void Dslash3D<double>::operator() (double* res, 
				    double* psi, 
				    double *u, /* Gauge field suitably packed */
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


  namespace DslashScalar64Bit {
    
    void DPsiPlus3D(size_t lo, size_t hi, int id, const void *ptr)
    {
      int ix,iy,iz;                                   /* Ix is current site */
      /* Iy is a neighbour */
      /* Iz is next site in loop */

      const ThreadWorkerArgs *a = (const ThreadWorkerArgs*)ptr;
      
      ShiftTable3D *shift = (ShiftTable3D*)a->s;
      int total_vol_cb = shift->totalVolCB();
      int cb = a->cb;

      GaugeMatrix (*gauge_field)[4] ALIGN =(GaugeMatrix(*)[4])a->u;/* My gauge field */
      FourSpinor *psi = (FourSpinor*)a->psi;                        /* Source */
      FourSpinor *res = (FourSpinor*)a->res;         /* Result */
      
      GaugeMatrix *up,*um;                               /* Pointer to FORWARD neighbour */
      FourSpinor *s,*sp,*sm,*rn;                       /* Pointer to BACKWARD neighbour */

      FourSpinor temp;
      int fsite, bsite, thissite;

   
      /******** loop over all lattice sites *************************/  
  
      /* Get forward neighbour of low in the x direction */
      const int low  =  cb*total_vol_cb+lo;                 /* First site for this thread */
      const int high  = cb*total_vol_cb+hi;                /* Last site for this thread */
      
      for (ix=low;ix<high;ix++) {
	thissite = shift->siteTable( ix );
	fsite = shift->forwardNeighbor( ix, 0);
	bsite = shift->backwardNeighbor( ix, 0);
	
	/* Result... Why do I need an RN? */
	rn=&res[thissite];
	
	/* Dir = 0 */
	sp=&psi[ fsite ];
	up=&(gauge_field[thissite][0]);
	dslash_plus_dir0_forward(*sp,*up,*rn);
	
	sm=&psi[ bsite ];
	um=&gauge_field[ bsite ][0];
	dslash_plus_dir0_backward_add(*sm,*um,*rn);
	
	fsite = shift->forwardNeighbor(  ix, 1);
	bsite = shift->backwardNeighbor(  ix, 1);
	
	sp=&psi[ fsite ];
	up=&(gauge_field[thissite][1]);
	dslash_plus_dir1_forward_add(*sp,*up,*rn);
	
	sm=&psi[ bsite ];
	um=&gauge_field[ bsite ][1];
	dslash_plus_dir1_backward_add(*sm,*um,*rn);
	
	fsite = shift->forwardNeighbor(  ix, 2);
	bsite = shift->backwardNeighbor(  ix, 2);
	
	sp=&psi[ fsite ];
	up=&(gauge_field[thissite][2]);
	dslash_plus_dir2_forward_add(*sp,*up,*rn);
	
	sm=&psi[ bsite ];
	um=&gauge_field[ bsite ][2];
	dslash_plus_dir2_backward_add(*sm,*um,*rn);
      }
    }
      
    void DPsiMinus3D(size_t lo, size_t hi, int id, const void *ptr)
    {
      int ix,iy,iz;                          /* ix is the current site */
      /* iy is the neighbour for prefetching */
      /* iz is the prefetch site for the 
	 next loop iteration */
      
      const ThreadWorkerArgs *a = (const ThreadWorkerArgs*)ptr;    /* Cast the void args pointer */

      ShiftTable3D *shift = (ShiftTable3D*)a->s;
      int total_vol_cb = shift->totalVolCB();
      const int cb = a->cb;
      
      
      GaugeMatrix ALIGN (*gauge_field)[4] =(GaugeMatrix (*)[4])a->u; /* Gauge field */
      FourSpinor *psi = (FourSpinor*)a->psi;                 /* Source spinor */
      FourSpinor *res = (FourSpinor*)a->res;                 /* Result spinor */
      GaugeMatrix *up,*um;                        /* us for multiply (PLUS/MINUS) */
      FourSpinor *sp,*sm,*rn;                   /* spinor pointers sp sm are the 
						     neighbours, rn is the result */
      
      
      /************************ loop over all lattice sites *************************/
      const int low  =  cb*total_vol_cb+lo;                 /* First site for this thread */
      const int high  = cb*total_vol_cb+hi;                /* Last site for this thread */
      int fsite, bsite, thissite;
      
      for (ix=low;ix<high;ix++) {
	thissite = shift->siteTable( ix );
	fsite = shift->forwardNeighbor(  ix, 0);
	bsite = shift->backwardNeighbor(  ix, 0);
	
	/* Result... Why do I need an RN? */
	rn=&res[thissite];
	
	/* Dir = 0 */
	sp=&psi[ fsite ];
	up=&(gauge_field[thissite][0]);
	dslash_minus_dir0_forward(*sp,*up,*rn);
	
	sm=&psi[ bsite ];
	um=&gauge_field[ bsite ][0];
	dslash_minus_dir0_backward_add(*sm,*um,*rn);
	
	fsite = shift->forwardNeighbor(  ix, 1);
	bsite = shift->backwardNeighbor(  ix, 1);
	
	sp=&psi[ fsite ];
	up=&(gauge_field[thissite][1]);
	dslash_minus_dir1_forward_add(*sp,*up,*rn);
	
	sm=&psi[ bsite ];
	um=&gauge_field[ bsite ][1];
	dslash_minus_dir1_backward_add(*sm,*um,*rn);
	
	fsite = shift->forwardNeighbor( ix, 2);
	bsite = shift->backwardNeighbor( ix, 2);
	
	sp=&psi[ fsite ];
	up=&(gauge_field[thissite][2]);
	dslash_minus_dir2_forward_add(*sp,*up,*rn);
	
	sm=&psi[ bsite ];
	um=&gauge_field[ bsite ][2];
	dslash_minus_dir2_backward_add(*sm,*um,*rn);
	
	
      }
    }


  } // namespace 64 bit 

} // Namespace
