/*******************************************************************************
 * $Id: sse_su3dslash_64bit_scalar_3d.c,v 1.3 2008-03-05 19:45:13 bjoo Exp $
 * 
 * Action of the 64bit single-node Wilson-Dirac operator D_w on a given spinor field
 *
 * The externally accessible function is   sse_su3dslash_wilson
 * 
 * void sse_su3dslash_wilson(float *u, float *psi, float *res, int isign, int cb)
 *
 * NOTE:
 * u: base pointer to gauge field
 * psi: base pointer to input spinor field on FULL lattice
 * res: base pointer to output spinor field on FULL lattice
 * isign: 1-->normal, -1--> swaps 1 - gamma(mu^) for 1 + gamma(mu^)
 * cb: checkerboard (0/1) of input fields
 *
 * Oringal author: Chris McClendon <cmcclend@jlab.org>
 * Acknowledgements to: Martin Luescher <luscher@mail.desy.de>
 * Date: 9/15/2001
 *
 *******************************************************************************/


#include <sse_config.h>

#ifdef __cplusplus
extern "C" {
#endif


#include <stdlib.h>
#include <stdio.h>
#include <math.h>


/* externally callable routine: tnxtsu3dslash */
/* requires gauge fields packed by no_funnystuff_pack_gauge_field of packer64.c */
/* This routine applies the operator D' to Psi, putting the result in Res. */

/*	       Nd-1 */
/*	       --- */
/*	       \ */
/*   res(x)  :=  >  U  (x) (1 - isign gamma  ) psi(x+mu) */
/*	       /    mu			  mu */
/*	       --- */
/*	       mu=0 */

/*	             Nd-1 */
/*	             --- */
/*	             \    + */
/*                +    >  U  (x-mu) (1 + isign gamma  ) psi(x-mu) */
/*	             /    mu			   mu */
/*	             --- */
/*	             mu=0 */

/* Arguments: */

/*  U	      Gauge field					(Read) */
/*  Psi	      Pseudofermion field				(Read) */
/*  Res	      Pseudofermion field				(Write) */
/*		      + */
/*  ISign      D' or D'  ( +1 | -1 ) respectively		(Read) */
/*  CB	      Checkerboard of input vector			(Read) */


#include <sse_align.h>
#include <shift_tables_scalar.h>
#include <dispatch_scalar.h>
#include <site_dslash_64bit_scalar.h>


  
  static int initP_3d = 0;          /* Reference count */
  static int *shift_table;       /* Shift Table */

  extern int* site_table_3d;
  extern int total_vol_3d_cb;


/***************** start of initialization routine ***************************************/
  void init_sse_su3dslash_3d(const int latt_size[],
			     void (*getSiteCoords)(int coord[], int node, int linearsite),
			 
			     int (*getLinearSiteIndex)(const int coord[]),
			     int (*nodeNumber)(const int coord[]))
{
  int mu;

  /* If we are already initialised, then increase the refcount and return */
  if (initP_3d > 0) {

    initP_3d++;
    return;
  }

  /* Check problem size - 4D  */
  for(mu=0; mu < 3; mu++)  {
    if ( latt_size[mu] % 2 != 0 ) {
      fprintf(stderr,"This is a Dslash with checkerboarding in 3 dimensions. GLOBAL dimension 0,1,2 (corresponding to x,y,z) must be even,  Your lattice is not like this: latt_size[%d]=%d\n", 
	      mu, latt_size[mu]);
      
      exit(1);
    }  
  }

  /* If make_shift_tables cannot allocate, it will barf */
  shift_table = make_shift_tables_3d(latt_size,
				     getSiteCoords,
				     getLinearSiteIndex);

  /* Set init _flag */
  initP_3d = 1;
}


void free_sse_su3dslash_3d(void)
{
  /* If we are uninitialised just return */
  if (initP_3d == 0) {
    return;
  }

  /* Otherwise decrease the refcount */
  initP_3d--;

  /* If the refcount has now hit 0 then free stuff */
  if( initP_3d == 0 ) {
    free_shift_tables(&shift_table);
  }

}

/***************** end of initialization routine ***************************************/


/* prototypes - for Thread Slaves */

void D_psi_fun_plus_3d(size_t lo,size_t hi, int id, const void *ptr);
void D_psi_fun_minus_3d(size_t lo,size_t hi, int id, const void *ptr);



/* External routine */
void sse_su3dslash_wilson_3d(double *u, double *psi, double *res, int isign, int cb)
{

  if (isign == 1)  {
    dispatch_to_threads(D_psi_fun_plus_3d, 
		(spinor_array*)psi,
		(spinor_array*)res,
		(my_mat_array)u,
		1-cb,
		total_vol_3d_cb);
  }

  if( isign == -1) {
    dispatch_to_threads(D_psi_fun_minus_3d, 
		(spinor_array*)psi,
		(spinor_array*)res,
		(my_mat_array)u,
		1-cb,
		total_vol_3d_cb);
  }
}

 


void D_psi_fun_plus_3d(size_t lo, size_t hi, int id, const void *ptr)
{
  int ix,iy,iz;                                      /* Ix is corrent site */
                                                     /* Iy is a neighbour */
                                                     /* Iz is next site in loop */

  const ThreadWorkerArgs *a = (const ThreadWorkerArgs*)ptr;                /* Downcast argument */
  int cb = a->cb;

  u_mat_array (*gauge_field)[4] ALIGN = a->u;        /* My gauge field */
  spinor_array *psi = a->psi;                        /* Source */
  spinor_array *res = a->res;                        /* Result */

  u_mat_array *up,*um;                               /* Pointer to FORWARD neighbour */
  spinor_array *s,*sp,*sm,*rn;                       /* Pointer to BACKWARD neighbour */

  spinor_array temp;
  int fsite, bsite, thissite;

   
  /************************ loop over all lattice sites *************************/  
  
  /* Get forward neighbour of low in the x direction */
  const int low  =  cb*total_vol_3d_cb+lo;                 /* First site for this thread */
  const int high  = cb*total_vol_3d_cb+hi;                /* Last site for this thread */

  for (ix=low;ix<high;ix++) 
  {
    thissite = site_table_3d[ ix ];
    fsite = forward_neighbor_3d( shift_table, ix, 0);
    bsite = backward_neighbor_3d( shift_table, ix, 0);

    /* Result... Why do I need an RN? */
    rn=&res[thissite];

    /* Dir = 0 */
    sp=&psi[ fsite ];
    up=&(gauge_field[thissite][0]);
    dslash_plus_dir0_forward(*sp,*up,*rn);

    sm=&psi[ bsite ];
    um=&gauge_field[ bsite ][0];
    dslash_plus_dir0_backward_add(*sm,*um,*rn);

    fsite = forward_neighbor_3d( shift_table, ix, 1);
    bsite = backward_neighbor_3d( shift_table, ix, 1);

    sp=&psi[ fsite ];
    up=&(gauge_field[thissite][1]);
    dslash_plus_dir1_forward_add(*sp,*up,*rn);

    sm=&psi[ bsite ];
    um=&gauge_field[ bsite ][1];
    dslash_plus_dir1_backward_add(*sm,*um,*rn);

    fsite = forward_neighbor_3d( shift_table, ix, 2);
    bsite = backward_neighbor_3d( shift_table, ix, 2);

    sp=&psi[ fsite ];
    up=&(gauge_field[thissite][2]);
    dslash_plus_dir2_forward_add(*sp,*up,*rn);

    sm=&psi[ bsite ];
    um=&gauge_field[ bsite ][2];
    dslash_plus_dir2_backward_add(*sm,*um,*rn);
  }
}



void D_psi_fun_minus_3d(size_t lo, size_t hi, int id, const void *ptr )
{
  int ix,iy,iz;                          /* ix is the current site */
                                         /* iy is the neighbour for prefetching */
                                         /* iz is the prefetch site for the 
					    next loop iteration */

  const ThreadWorkerArgs *a = (const ThreadWorkerArgs*)ptr;    /* Cast the void args pointer */
  const int cb = a->cb;


  u_mat_array (*gauge_field)[4] ALIGN = a->u; /* Gauge field */
  spinor_array *psi = a->psi;                 /* Source spinor */
  spinor_array *res = a->res;                 /* Result spinor */
  u_mat_array *up,*um;                        /* us for multiply (PLUS/MINUS) */
  spinor_array *sp,*sm,*rn;                   /* spinor pointers sp sm are the 
						 neighbours, rn is the result */

   
/************************ loop over all lattice sites *************************/
  const int low  =  cb*total_vol_3d_cb+lo;                 /* First site for this thread */
  const int high  = cb*total_vol_3d_cb+hi;                /* Last site for this thread */
  int fsite, bsite, thissite;

  for (ix=low;ix<high;ix++) {
    thissite = site_table_3d[ ix ];
    fsite = forward_neighbor_3d( shift_table, ix, 0);
    bsite = backward_neighbor_3d( shift_table, ix, 0);

    /* Result... Why do I need an RN? */
    rn=&res[thissite];

    /* Dir = 0 */
    sp=&psi[ fsite ];
    up=&(gauge_field[thissite][0]);
    dslash_minus_dir0_forward(*sp,*up,*rn);

    sm=&psi[ bsite ];
    um=&gauge_field[ bsite ][0];
    dslash_minus_dir0_backward_add(*sm,*um,*rn);

    fsite = forward_neighbor_3d( shift_table, ix, 1);
    bsite = backward_neighbor_3d( shift_table, ix, 1);

    sp=&psi[ fsite ];
    up=&(gauge_field[thissite][1]);
    dslash_minus_dir1_forward_add(*sp,*up,*rn);

    sm=&psi[ bsite ];
    um=&gauge_field[ bsite ][1];
    dslash_minus_dir1_backward_add(*sm,*um,*rn);

    fsite = forward_neighbor_3d( shift_table, ix, 2);
    bsite = backward_neighbor_3d( shift_table, ix, 2);

    sp=&psi[ fsite ];
    up=&(gauge_field[thissite][2]);
    dslash_minus_dir2_forward_add(*sp,*up,*rn);

    sm=&psi[ bsite ];
    um=&gauge_field[ bsite ][2];
    dslash_minus_dir2_backward_add(*sm,*um,*rn);


  }
}


#ifdef __cplusplus
}
#endif



