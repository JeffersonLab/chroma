#include <cpp_dslash_scalar.h>
#include <shift_table_scalar.h>
#include <dispatch_scalar.h>

#include <cpp_dslash_scalar_64bit.h>

using namespace CPlusPlusWilsonDslash::DslashScalar64Bit;
using namespace CPlusPlusWilsonDslash::Dslash64BitTypes;

namespace CPlusPlusWilsonDslash {
 

  /* Constructor */
  Dslash<double>::Dslash(const int latt_size[],      
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

  Dslash<double>::~Dslash() { delete s; }
  
  //  int Dslash<double>::getPathSite(int site) const
  // { 
  //  return s->getPathSite(site);
  // }
 
  // The operator 
  void Dslash<double>::operator() (double* res, 
				   double* psi, 
				   double *u, /* Gauge field suitably packed */
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


  namespace DslashScalar64Bit {
       void DPsiPlus(size_t lo, size_t hi, int id, const void *ptr)
       {
	 int ix,iy,iz;                                      /* Ix is corrent site */
	 /* Iy is a neighbour */
	 /* Iz is next site in loop */

	 const ThreadWorkerArgs *a = (const ThreadWorkerArgs*)ptr;                /* Downcast argument */
	 ShiftTable *shift = (ShiftTable *)a->s;
	 int total_vol_cb = shift->totalVolCB();
	 int cb = a->cb;
	 const int low = cb*total_vol_cb+lo;
	 const int high = cb*total_vol_cb+hi;
	 
	 GaugeMatrix ALIGN (*gauge_field)[4] = (GaugeMatrix(*)[4])a->u;        /* My gauge field */
	 FourSpinor *psi = (FourSpinor *)a->psi;                        /* Source */
	 FourSpinor *res = (FourSpinor *)a->res;                        /* Result */

	 GaugeMatrix *up,*um;                               /* Pointer to FORWARD neighbour */
	 FourSpinor *s,*sp,*sm,*rn;                       /* Pointer to BACKWARD neighbour */
	 
	 FourSpinor temp;
	 int thissite;
	 
	 //thissite=shift->getPathSite(low);
	 thissite=shift->siteTable(low);
	 /* This is like a prefetch 
	    - we peel it off the loop */
	 
	 /* Get 4 spinor from forward direction */
	 sp=&psi[shift->forwardNeighbor(low,0) ];
	 
	 
	 /* Get Gauge Field */
	 up=&(gauge_field[low][0]);
	 
	 /************************ loop over all lattice sites *************************/
	 for (ix=low;ix<high;ix++) 
	   {
	     /******************************* direction +0 *********************************/
	     rn=&res[ix];
	     
	     /******************************   Direction 0  ********************* */
	     /* Prefetch back spinor for next dir: -0 */
	     iy=shift->backwardNeighbor(ix,0);
	     sm=&psi[iy];

	     
	     /* Prefetch back gauge field */
	     um=&(gauge_field[iy][0]);
	     
	     dslash_plus_dir0_forward(*sp,*up,*rn);
	     
	     
	     /* sm and um should already be prefetched */
	     
	     /* Now prefetch forward neighbor for next dir (1+) */
	     /* And gauge field */
	     sp=&psi[ shift->forwardNeighbor(ix,1) ];

	     up =&(gauge_field[ix][1]);
	     
	     dslash_plus_dir0_backward_add(*sm,*um,*rn);
	     
	     
	     /********************** Direction 1 ************************ */
	     /* up and sp should be prefetched */
	     
	     /* Prefetch backwards spinor and gauge field */
	     iy=shift->backwardNeighbor(ix,1);
	     sm=&psi[iy];

	     um=&(gauge_field[iy][1]);

	     
	     
	     dslash_plus_dir1_forward_add(*sp,*up,*rn);
	     
	     
	     
	     /* Prefetch forwards spinor and gauge field for next dir: 2+ */
	     iy=shift->forwardNeighbor(ix,2);
	     sp=&psi[iy];

	     up = &(gauge_field[ix][2]);

	     
	     
	     dslash_plus_dir1_backward_add(*sm,*um,*rn);
	     
	     
	     
	     /********************** Direction 2 ************************* */
	     /* Prefetch back spinor and gauge field  */
	     
	     iy=shift->backwardNeighbor(ix,2);
	     sm=&psi[iy];

	     um=&(gauge_field[iy][2]);


	     
	     dslash_plus_dir2_forward_add(*sp,*up,*rn);
	     
	     
	     /* Prefetch forward spinor and gauge field for next direction: 3+ */
	     iy=shift->forwardNeighbor(ix,3);
	     sp=&psi[iy];

	     up = &(gauge_field[ix][3]);

	     
	     dslash_plus_dir2_backward_add(*sm,*um,*rn);
	     
	     
	     /* ******************* Direction +3 **************************** */
	     /* Prefetch back spinor and gauge field  3- */
	     
	     iy=shift->backwardNeighbor(ix,3);
	     sm=&psi[iy];

	     um=&(gauge_field[iy][3]);

	     
	     dslash_plus_dir3_forward_add(*sp,*up,*rn);      
	     
	     /* Next site */
	     iz=ix+1;
	     // thissite = shift->getPathSite(iz);
	     thissite = shift->siteTable(iz);
	     if (iz == high) { /* If we're on the last site, prefetch first site to avoid */
	       iz=0;           /* Running beyond array bounds */
	     }
	     
	     
	     /* Prefetch forward spinor and gauge field for next site, dir 0 */
	     iy=shift->forwardNeighbor(iz,0);
	     sp=&psi[iy];

	     up=&(gauge_field[iz][0]);

	     
	     
	     dslash_plus_dir3_backward_add_store(*sm,*um,*rn);      
	     
	   }
       }



    void DPsiMinus(size_t lo, size_t hi, int id, const void *ptr )
    {
      int ix,iy,iz;                          /* ix is the current site */
      /* iy is the neighbour for prefetching */
      /* iz is the prefetch site for the 
	 next loop iteration */
      
      const ThreadWorkerArgs *a = (const ThreadWorkerArgs*)ptr;    /* Cast the void args pointer */
      const int cb = a->cb;
      ShiftTable *shift = (ShiftTable *)a->s;
      int total_vol_cb = shift->totalVolCB();
      GaugeMatrix (*gauge_field)[4] ALIGN = (GaugeMatrix (*)[4])a->u; /* Gauge field */
      FourSpinor *psi = (FourSpinor *)a->psi;                 /* Source spinor */
      FourSpinor *res = (FourSpinor *)a->res;                 /* Result spinor */
      GaugeMatrix *up,*um;                        /* us for multiply (PLUS/MINUS) */
      FourSpinor *sp,*sm,*rn;                   /* spinor pointers sp sm are the 
						   neighbours, rn is the result */
      
      /* Get forward neighbour of low in the x direction */
      const int low  =  cb*total_vol_cb+lo;                 /* First site for this thread */
      const int high  = cb*total_vol_cb+hi;                /* Last site for this thread */
      int thissite;
      
      //      thissite=shift->getPathSite(low);
      thissite = shift->siteTable(low);
      /* 'peel this off' to allow prefetching */
      iy=shift->forwardNeighbor(low,0);
      sp=&psi[iy];
      up=&(gauge_field[low][0]);
      
      /************************ loop over all lattice sites *************************/
      
      for (ix=low;ix<high;ix++) {
	rn=&res[ix]; /* Pouinter to result */
	
	
	/* Prefetch back spinor and gauge field */
	iy=shift->backwardNeighbor(ix,0);
	sm=&psi[iy];

	um=&(gauge_field[iy][0]);

	
	dslash_minus_dir0_forward(*sp,*up,*rn);
	
	
	/* Prefetch spinor and gauge field for next direction (1) */
	iy=shift->forwardNeighbor(ix,1);
	sp=&psi[iy];

	up = &(gauge_field[ix][1]);

	
	dslash_minus_dir0_backward_add(*sm,*um,*rn);    
	
	/* ******************* direction +1 ******************* */
	/* Prefetch backward spinor and gauge field */
	
	iy=shift->backwardNeighbor(ix,1);
	sm=&psi[iy];

	um=&(gauge_field[iy][1]);

	
	dslash_minus_dir1_forward_add(*sp,*up,*rn);
	
        
	/* Prefetch forward spinor and gauge field for next direction: 2 */
	iy=shift->forwardNeighbor(ix,2);
	sp=&psi[iy];

	up =&(gauge_field[ix][2]);

	
	dslash_minus_dir1_backward_add(*sm,*um,*rn);    
	
	
	/* ******************* direction +2 **************************** */
	/* Prefetch backward spinor and gauge field */
	
	iy=shift->backwardNeighbor(ix,2); 
	sm=&psi[iy]; 

	um=&(gauge_field[iy][2]);

	
	dslash_minus_dir2_forward_add(*sp,*up,*rn);    
        
	
	/* Prefetch spinor and gauge field for next direction: 3+ */
	iy=shift->forwardNeighbor(ix,3);
	sp=&psi[iy];

	
	/* Prefetch gauge field for nex dir: 3+ */
	up = &(gauge_field[ix][3]);

	
	dslash_minus_dir2_backward_add(*sm,*um,*rn);    
	
        
	
	/*******************  direction +3 ************************** */
	/* Prefetch backward spinor and gauge field */
	iy=shift->backwardNeighbor(ix,3);
	sm=&psi[iy];

	um=&(gauge_field[iy][3]);

	
	dslash_minus_dir3_forward_add(*sp,*up,*rn);    
	
	
	/* Next site in loop. We peeled this off the loop to start with so we can prefetch... */
	iz=ix+1; 
	if (iz==high) { /* If we are on the last site, we should prefetch the first element, to
			   avoid running past the array bounds */
	  iz=0;
	}
	//	thissite = shift->getPathSite(iz);
	thissite = shift->siteTable(iz);

	/* Prefetch the spinor and gauge field for next site, dir 0+ */
	iy=shift->forwardNeighbor(iz,0);
	sp=&psi[iy];

	/* Prefetch the gauge field for next site dir  0+ */
	up=&(gauge_field[iz][0]);

	
	dslash_minus_dir3_backward_add_store(*sm,*um,*rn);    
	
      }
    }

  }

} // Namespace
