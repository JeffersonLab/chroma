/*******************************************************************************
 * $Id: sse_su3dslash_32bit_scalar_3d.c,v 1.3 2008-03-05 19:45:13 bjoo Exp $
 * 
 * Action of the 32bit single-node Wilson-Dirac operator D_w on a given spinor field
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
 * Date: 9/11/2001
 *
 *******************************************************************************/

#include <sse_config.h>
#include <sse_align.h>
#include <shift_tables_scalar.h>
#include <types32.h>
#include <dispatch_scalar.h>
#include <site_dslash_32bit_scalar.h>
#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>



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

  /* Volume and initialization */
  static int initP_3d = 0;

  /* These are needed for the shift tables */
  static int *shift_table;

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
    if (initP_3d > 0)  {
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

  /* Construct all the shift tables needed */
  /* 4 dimensions * 2 directions { aka FORWARD and BACKWARD } * volume */
    /* shift_table and icolor start are set, latt_size is read */
  shift_table = make_shift_tables_3d(latt_size,
				     getSiteCoords,
				     getLinearSiteIndex);

  /* Set flag to signify initialization */
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

  /* If the refcount has now hit 0 then free the shift tables*/
  if( initP_3d == 0 ) {
    free_shift_tables_3d(&shift_table);
  }
}

/***************** end of initialization routine ***************************************/


/* Apply Dslash */
void D_psi_fun_plus_3d(size_t lo,size_t hi, int id, const void *ptr);

/* Apply Dslash Dagger */
void D_psi_fun_minus_3d(size_t lo,size_t hi, int id, const void *ptr);


/* routine sse_su3dslash_wilson
   u: base pointer to gauge field
   psi: base pointer to input spinor field on full lattice
   res: base pointer to output spinor field on full lattice
   isign: 1-->normal, -1--> swaps 1 - gamma(mu^) for 1 + gamma(mu^)
   cb: checkerboard (0/1) of input fields
*/
void sse_su3dslash_wilson_3d(float *u, float *psi, float *res, int isign, int cb)
{

  if (isign == 1) {
    dispatch_to_threads(D_psi_fun_plus_3d, 
			(spinor_array*)psi,
			(spinor_array*)res,
			(my_mat_array)u,
			1-cb,
			total_vol_3d_cb);
  }

  if( isign == -1) 
  {
    dispatch_to_threads(D_psi_fun_minus_3d, 
			(spinor_array*)psi,
			(spinor_array*)res,
			(my_mat_array)u,
			1-cb,
			total_vol_3d_cb);
  }
}

// #include <sse32.h>


void D_psi_fun_plus_3d(size_t lo,size_t hi, int id, const void *ptr)
{

  const ThreadWorkerArgs *a  = (const ThreadWorkerArgs*)ptr;  /* Cast the (void *) to an (ThreadWorkerArgs*) */
  int ix1;                              /* Index of current site */
  int iy1,iy2;                          /* Index of neighbour of current site (iy1) 
					   and current site+1 (iy2) */

  int iz1;                              /* Index of next site in loop */

  u_mat_array (*gauge_field)[4]  =  a->u; /* Packed Gauge fields */
  spinor_array *psi  =  a->psi;           /* Source spinor */
  spinor_array *res  =  a->res;           /* Result spinor */

  const int cb = a->cb;


  

  /* Pointers to the neighboring u-s */
  u_mat_array *up1 ALIGN;                  /* U[ x  ] */
  u_mat_array *up2 ALIGN;                  /* U[ (x+1)  ] */
  u_mat_array *um1 ALIGN;                  /* U[ x - mu ] */
  u_mat_array *um2 ALIGN;                  /* U[ (x+1)-mu ] */

  /* 4 - Spinor pointers */
  spinor_array *sp1 ALIGN;
  spinor_array *sm1 ALIGN;
  spinor_array *sn1 ALIGN;

  /* Half Vectors */
  halfspinor_array r12_1 ALIGN; /* Site 1 upper */
  halfspinor_array r34_1 ALIGN; /* Site 1 lower */


  
  /* note that we want the spinors from the checkerboard opposite the one we are writing to */
  /* We are doing things in bunches of two sites */

  /* Get forward neighbour of low in the x direction */
  const int low  =  cb*total_vol_3d_cb+lo;                 /* First site for this thread */
  const int high  = cb*total_vol_3d_cb+hi;                /* Last site for this thread */

  for (ix1 = low; ix1< high; ix1++) 
  {
    int thissite = site_table_3d[ ix1 ];
    int fsite = forward_neighbor_3d(shift_table,ix1,0);
    int bsite = backward_neighbor_3d(shift_table,ix1,0);

    /* Dir = 0 */
    sp1 = &psi[ fsite ];
    up1 = &(gauge_field[thissite][0]);
    dslash_plus_dir0_forward(*sp1, *up1, r12_1, r34_1);

    um1 = &(gauge_field[bsite][0]);       
    sm1 = &psi[bsite];
    dslash_plus_dir0_backward_add(*sm1, *um1, r12_1, r34_1);

    /* Dir = 1 */
    fsite = forward_neighbor_3d(shift_table,ix1,1);
    bsite = backward_neighbor_3d(shift_table,ix1,1);

    sp1 = &psi[ fsite ];
    up1 = &(gauge_field[thissite][1]);
    dslash_plus_dir1_forward_add(*sp1, *up1, r12_1, r34_1);

    um1 = &(gauge_field[bsite][1]);
    sm1 = &psi[bsite];
    dslash_plus_dir1_backward_add(*sm1, *um1, r12_1, r34_1);


    /* Dir = 2 */
    fsite = forward_neighbor_3d(shift_table,ix1,2);
    bsite = backward_neighbor_3d(shift_table,ix1,2);

    sp1 = &psi[fsite ];
    up1 = &(gauge_field[thissite][2]);
    dslash_plus_dir2_forward_add(*sp1, *up1, r12_1, r34_1);

    sm1 = &psi[bsite];
    um1 = &(gauge_field[bsite][2]);
    sn1 = &res[thissite];
    dslash_plus_dir2_backward_add_store(*sm1, *um1, r12_1, r34_1, *sn1);

#if 0
    /* Dir = 3 */
    fsite = forward_neighbor_3d(shift_table,ix1,3);
    bsite = backward_neighbor_3d(shift_table,ix1,3);

    sp1 = &psi[ fsite ];
    up1 = &(gauge_field[thissite][3]);
    dslash_plus_dir3_forward_add(*sp1, *up1, r12_1, r34_1);

    sm1 = &psi[bsite]; 
    um1 = &(gauge_field[bsite][3]); 


    dslash_plus_dir3_backward_add_store(*sm1, *um1, r12_1, r34_1, *sn1);
#endif      
 
  }
}

/*ok, now this routine here is just like the one above except isign has a different value, which means that the
signs used for 1 +- gamma(mu^) must be swapped...*/ 
void D_psi_fun_minus_3d(size_t lo,size_t hi, int id, const void *ptr)
{
  const ThreadWorkerArgs *a  = (const ThreadWorkerArgs*)ptr;   /* Downcast to args */
  int ix1,iy1,iy2,iz1;                   /* Coordinates ix1 - current
					    iy1 - index of neighbour
					    iy1 - index of neighbour of ix1+1 
					    iz1 - index of first of next pair (ix+2) */
  const int cb = a->cb;

  u_mat_array (*gauge_field)[4]  =  a->u;  /* Gauge field */

  spinor_array *psi  =  a->psi;            /* Source 4-spinor */
  spinor_array *res  =  a->res;            /* Result 4-spinor */

  u_mat_array *up1 ALIGN;                  /* U[ ix ] */
  u_mat_array *um1 ALIGN;                  /* U[ ix - mu ] */

  spinor_array *sp1 ALIGN;                 /* 4 spinor psi[ix1+mu] */
  spinor_array *sm1 ALIGN;                 /* 4 spinor psi[ix1-mu] */
  spinor_array *sn1 ALIGN;                 /* 4 spinor result */

  halfspinor_array r12_1;                         /* site 1 halfspinor top half */
  halfspinor_array r34_1;                         /* site 1 halfspinor bottom half */

  /* Get forward neighbour of low in the x direction */
  const int low  =  cb*total_vol_3d_cb+lo;                 /* First site for this thread */
  const int high  = cb*total_vol_3d_cb+hi;                /* Last site for this thread */

  for (ix1 = low; ix1<high; ix1++) 
  {
    int thissite = site_table_3d[ ix1 ];
    int x_plus_mu = forward_neighbor_3d(shift_table,ix1,0);
    /* Dir = 0 */
    sp1 = &psi[ forward_neighbor_3d(shift_table,ix1,0)  ];
    up1 = &(gauge_field[thissite][0]);
    dslash_minus_dir0_forward(*sp1, *up1, r12_1, r34_1);

    iy1 =  backward_neighbor_3d(shift_table,ix1,0);
    um1 = &(gauge_field[iy1][0]);       
    sm1 = &psi[ iy1 ];
    dslash_minus_dir0_backward_add(*sm1, *um1, r12_1, r34_1);

    /* Dir = 1 */
    sp1 = &psi[ forward_neighbor_3d(shift_table,ix1,1) ];
    up1 = &(gauge_field[thissite][1]);
    dslash_minus_dir1_forward_add(*sp1, *up1, r12_1, r34_1);

    iy1 = backward_neighbor_3d(shift_table,ix1,1);
    um1 = &(gauge_field[iy1][1]);
    sm1 = &psi[iy1];
    dslash_minus_dir1_backward_add(*sm1, *um1, r12_1, r34_1);


    /* Dir = 2 */
    sp1 = &psi[forward_neighbor_3d(shift_table,ix1,2)];
    up1 = &(gauge_field[thissite][2]);
    dslash_minus_dir2_forward_add(*sp1, *up1, r12_1, r34_1);

    iy1 = backward_neighbor_3d(shift_table,ix1,2);
    sm1 = &psi[iy1];
    um1 = &(gauge_field[iy1][2]);
    sn1 = &res[thissite];
    dslash_minus_dir2_backward_add_store(*sm1, *um1, r12_1, r34_1, *sn1);

#if 0
    sp1 = &psi[ forward_neighbor_3d(shift_table,ix1,3) ];
    up1 = &(gauge_field[thissite][3]);
    dslash_minus_dir3_forward_add(*sp1, *up1, r12_1, r34_1);

    iy1 = backward_neighbor_3d(shift_table,ix1,3);
    sm1 = &psi[iy1]; 
    um1 = &(gauge_field[iy1][3]); 


    dslash_minus_dir3_backward_add_store(*sm1, *um1, r12_1, r34_1, *sn1);
#endif
      
  }
}

#ifdef __cplusplus
}
#endif
