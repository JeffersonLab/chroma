/*******************************************************************************
 * $Id: sse_su3dslash_64bit_parscalar_3d.c,v 1.4 2008-03-05 19:45:13 bjoo Exp $
 * 
 * Action of the 32bit parallel Wilson-Dirac operator D_w on a given spinor field
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

/* requires gauge fields packed by pack_gauge_field of intpar_table.c */

/* externally callable function: wnxtsu3dslash */
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

/* U	     Gauge field					(Read) */
/* Psi	     Pseudofermion field				(Read) */
/* Res	     Pseudofermion field				(Write) */
/*		      + */
/* ISign     D' or D'  ( +1 | -1 ) respectively	                (Read) */
/* CB	     Checkerboard of input vector			(Read) */


#include <sse_config.h>

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <qmp.h>
#include <types64.h>
#include <sse_align.h>
#include <dispatch_parscalar.h>
#include <shift_tables_parscalar.h>
#include <decomp.h>
#include <decomp_hvv.h>
#include <recons.h>
#include <mvv_recons_64bit.h>

  static int initP_3d=0;

  extern int subgrid_vol_3d;
  extern int subgrid_vol_cb_3d;

  extern int* site_table_3d;
  extern int* offset_table_body_3d;

  extern int Nd3;


  static inline
  int halfspinor_buffer_offset(HalfSpinorOffsetType type, int site, int mu)
  {
    int LocalNd3=3;
    return offset_table_body_3d[mu + LocalNd3*( site + subgrid_vol_3d*type) ];
  }






  /* this routine is similar to wnxtsu3dslash, except instead of handling the second site's worth in the same loop, the second spin component's worth must be handled seperately */
void decomp_3d_plus(size_t lo, size_t hi, int id, const void *ptr)
{
  const ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr; 
  const int cb = a->cb; 
  spinor_array *psi ALIGN = a->spinor;
  halfspinor_array *chi ALIGN = a->half_spinor;
  
  int ix1;
  halfspinor_array *s3 ALIGN ;
  spinor_array* sp ALIGN;
  const int low = cb*subgrid_vol_cb_3d + lo; 
  const int high = cb*subgrid_vol_cb_3d + hi;


  for (ix1=low;ix1<high;++ix1) {
    int thissite = site_table_3d[ ix1 ];
    sp=&psi[thissite];
    s3 = chi + halfspinor_buffer_offset(DECOMP_SCATTER,ix1,0);
    decomp_gamma0_minus(*sp, *s3);
    /******************************* direction +1 *********************************/
    s3 = chi + halfspinor_buffer_offset(DECOMP_SCATTER,ix1,1);
    decomp_gamma1_minus(*sp, *s3);
    
    /******************************* direction +2 *********************************/
    s3 = chi + halfspinor_buffer_offset(DECOMP_SCATTER,ix1,2);
    decomp_gamma2_minus(*sp, *s3);

  }
}


void decomp_hvv_3d_plus(size_t lo, size_t hi, int id, const void *ptr)
{
  const ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr; 
  const int cb = a->cb; 
  u_mat_array (*gauge_field)[4] = a->u;
  spinor_array *psi = a->spinor;
  halfspinor_array *chi = a->half_spinor;

  int ix1=0;

  u_mat_array *um ALIGN;

  halfspinor_array *s3 ALIGN;

  spinor_array *sm ALIGN; 
   


  const int low = cb*subgrid_vol_cb_3d + lo; 
  const int high = cb*subgrid_vol_cb_3d + hi;


  for (ix1=low;ix1<high;++ix1) {
    int thissite = site_table_3d[ ix1 ];

    sm=&psi[thissite];
     
    /******************************* direction -0 *********************************/
    um=&gauge_field[thissite][0];
    s3 = chi + halfspinor_buffer_offset(DECOMP_HVV_SCATTER,ix1,0);
    decomp_hvv_gamma0_plus(*sm, *um, *s3);
    
    /******************************* direction -1 *********************************/
    um=&gauge_field[thissite][1];
    s3 = chi + halfspinor_buffer_offset(DECOMP_HVV_SCATTER,ix1,1);
    decomp_hvv_gamma1_plus(*sm, *um, *s3);

    /******************************* direction -2 *********************************/
    um=&gauge_field[thissite][2];
    s3 = chi + halfspinor_buffer_offset(DECOMP_HVV_SCATTER,ix1,2);
    decomp_hvv_gamma2_plus(*sm, *um, *s3);
  }
}


void mvv_recons_3d_plus(size_t lo, size_t hi, int id, const void *ptr)
{
  const ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr; 
  const int cb = a->cb; 
  int ix1=0;

  u_mat_array (*gauge_field)[4] = a->u;
  spinor_array *psi = a->spinor;
  halfspinor_array *chi = a->half_spinor;
  u_mat_array *up ALIGN;
  halfspinor_array *s3 ALIGN, *s4 ALIGN;
  spinor_array part_sum ALIGN, *result ALIGN;
 


  const int low = cb*subgrid_vol_cb_3d + lo; 
  const int high = cb*subgrid_vol_cb_3d + hi;


  for (ix1=low;ix1<high;++ix1) {
    int thissite = site_table_3d[ ix1 ];

	 
    /******************************* direction +0 *********************************/	
    up=&gauge_field[thissite][0];    
    s3 = chi + halfspinor_buffer_offset(RECONS_MVV_GATHER,ix1,0);
    mvv_recons_gamma0_plus(*s3, *up, part_sum);

	   
    /******************************* direction +1 *********************************/
    up=&gauge_field[thissite][1];
    s3 = chi + halfspinor_buffer_offset(RECONS_MVV_GATHER,ix1,1);
    mvv_recons_gamma1_plus_add(*s3, *up, part_sum);


    /******************************* direction +2 *********************************/
    up=&gauge_field[thissite][2];
    s3 = chi + halfspinor_buffer_offset(RECONS_MVV_GATHER,ix1,2);
    mvv_recons_gamma2_plus_add_store(*s3, *up, part_sum, psi[thissite]);
  }
}



/*optimized for SZIN spin basis */
void recons_3d_plus(size_t lo, size_t hi, int id, const void *ptr)
{


  const ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr; 
  const int cb = a->cb; 
  int ix1=0;

  spinor_array *psi = a->spinor;
  halfspinor_array *chi = a->half_spinor;

  halfspinor_array *hs0 ALIGN , *hs1 ALIGN, *hs2 ALIGN;
  spinor_array *s1 ALIGN, *rn ALIGN;

 int low = cb*subgrid_vol_cb_3d + lo;
  int high = cb*subgrid_vol_cb_3d + hi;
  
  
  /************************ loop over all lattice sites *************************/
  for (ix1 =low; ix1 <high ;ix1++) {

    int thissite = site_table_3d[ ix1 ];

    /* first spin component of result */
    hs0 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1,0);
    hs1 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1,1);	  
    hs2 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1,2);
    recons_3dir_plus(*hs0, *hs1, *hs2, psi[thissite]);
  }
 
}


/************now for isign = -1  **********************/


void decomp_3d_minus(size_t lo, size_t hi, int id, const void *ptr)
{
   

  const ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr; 
  const int cb = a->cb; 
  int ix1=0;

  spinor_array *psi = a->spinor;
  halfspinor_array *chi = a->half_spinor;

  halfspinor_array *s3 ALIGN; 
  spinor_array *sp ALIGN;
 
  int low = cb*subgrid_vol_cb_3d + lo;
  int high = cb*subgrid_vol_cb_3d + hi;
  
  
  /************************ loop over all lattice sites *************************/
  for (ix1 =low; ix1 <high ;ix1++) {

    int thissite = site_table_3d[ ix1 ];       
    sp=&psi[thissite];
    
    s3 = chi + halfspinor_buffer_offset(DECOMP_SCATTER,ix1,0);
    decomp_gamma0_plus(*sp, *s3);

    s3 = chi + halfspinor_buffer_offset(DECOMP_SCATTER,ix1,1);
    decomp_gamma1_plus(*sp, *s3);
    
    s3 = chi + halfspinor_buffer_offset(DECOMP_SCATTER,ix1,2);
    decomp_gamma2_plus(*sp, *s3);
    
  }
}


void decomp_hvv_3d_minus(size_t lo, size_t hi, int id, const void *ptr)
{
  
  const ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr; 
  const int cb = a->cb; 
  int ix1=0;
  
  
  u_mat_array (*gauge_field)[4] = a->u;
  
  spinor_array *psi = a->spinor;
  halfspinor_array *chi = a->half_spinor;
  
  u_mat_array *um ALIGN;
  halfspinor_array *s3 ALIGN;
  spinor_array *s1 ALIGN, *sm ALIGN;
  
  
  int low = cb*subgrid_vol_cb_3d + lo;
  int high = cb*subgrid_vol_cb_3d + hi;
  
  for (ix1 =low; ix1 <high ;ix1++) {

    int thissite = site_table_3d[ ix1 ];    

    sm=&psi[thissite];
    um=&gauge_field[thissite][0];

    s3 = chi + halfspinor_buffer_offset(DECOMP_HVV_SCATTER,ix1,0);
    decomp_hvv_gamma0_minus(*sm, *um, *s3);	   

    um=&gauge_field[thissite][1];
    s3 = chi + halfspinor_buffer_offset(DECOMP_HVV_SCATTER,ix1,1);
    decomp_hvv_gamma1_minus(*sm, *um, *s3);	   

    um=&gauge_field[thissite][2];
    s3 = chi + halfspinor_buffer_offset(DECOMP_HVV_SCATTER,ix1,2);
    decomp_hvv_gamma2_minus(*sm, *um, *s3);	   

  }
}


void mvv_recons_3d_minus(size_t lo, size_t hi, int id, const void *ptr)
{
   	
  
  const ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr; 
  const int cb = a->cb; 
  int ix1=0;

  u_mat_array (*gauge_field)[4] = a->u;
  spinor_array *psi = a->spinor;
  halfspinor_array *chi = a->half_spinor;

  u_mat_array *up ALIGN;

  halfspinor_array *s3 ALIGN, *s4 ALIGN;

  spinor_array rs ALIGN,*rn ALIGN;

  
  int low = cb*subgrid_vol_cb_3d + lo;
  int high = cb*subgrid_vol_cb_3d + hi;
  
  for (ix1 =low; ix1 <high ;ix1++) {
    int thissite = site_table_3d[ ix1 ];       

    up=&gauge_field[thissite][0];
    s3 = chi + halfspinor_buffer_offset(RECONS_MVV_GATHER,ix1,0);
    mvv_recons_gamma0_minus(*s3, *up, rs);

    up = &gauge_field[thissite][1];
    s3 = chi + halfspinor_buffer_offset(RECONS_MVV_GATHER,ix1,1);
    mvv_recons_gamma1_minus_add(*s3, *up, rs);
	   
    up = &gauge_field[thissite][2];
    s3 = chi + halfspinor_buffer_offset(RECONS_MVV_GATHER,ix1,2);
    mvv_recons_gamma2_minus_add_store(*s3, *up, rs, psi[thissite]);
  }

}



void recons_3d_minus(size_t lo, size_t hi, int id, const void *ptr)
{
  
  const ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr; 
  const int cb = a->cb; 
  int ix1=0;


  spinor_array *psi = a->spinor;
  halfspinor_array *chi = a->half_spinor;

  halfspinor_array *hs0,*hs1,*hs2,*hs3;   
  spinor_array  *s1 ALIGN,  *rn ALIGN;

  int low = cb*subgrid_vol_cb_3d + lo;
  int high = cb*subgrid_vol_cb_3d + hi;
  
  
  /************************ loop over all lattice sites *************************/
  for (ix1 =low; ix1 <high ;ix1++) {

    int thissite = site_table_3d[ ix1 ];     
    
    rn=&psi[thissite];
    /* first spin component of result */
    hs0 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1,0);
    hs1 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1,1);
    hs2 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1,2);
    recons_3dir_minus(*hs0, *hs1, *hs2, psi[thissite]);
    /* end of loop */
  }
}

static QMP_mem_t* xchi1_3d;
static QMP_mem_t* xchi2_3d;
static halfspinor_array* chi1_3d;
static halfspinor_array* chi2_3d;


/* Nearest neighbor communication channels */
static int total_comm_3d = 0;
static QMP_msgmem_t forw_msg_3d[3][2];
static QMP_msgmem_t back_msg_3d[3][2];
static QMP_msghandle_t forw_mh_3d[3][2];
static QMP_msghandle_t back_mh_3d[3][2];
static QMP_msghandle_t forw_all_mh_3d;
static QMP_msghandle_t back_all_mh_3d;



  void init_sse_su3dslash_3d(const int latt_size[],
			     void (*getSiteCoords)(int coord[], int node, int linearsite),
			     
			     int (*getLinearSiteIndex)(const int coord[]),
			     int (*nodeNumber)(const int coord[]))   // latt_size not used, here for scalar version
{
  const int *machine_size = QMP_get_logical_dimensions();
  int bound[2][4][3];

  int mu, num, nsize;

  /* If we are already initialised, then increase the refcount and return */
  if (initP_3d > 0) 
  {
    initP_3d++;
    return;
  }


  /* Otherwise initialise */
  if (QMP_get_logical_number_of_dimensions() != 4)
  {
    QMP_error("init_sse_su3dslash: number of logical dimensions does not match problem");
    QMP_abort(1);
  }
    
 /* Check problem size - 3D  */
  for(mu=0; mu < 3; mu++)  {
    if ( latt_size[mu] % 2 != 0 ) {
      fprintf(stderr,"This is a Dslash with checkerboarding in 3 dimensions. GLOBAL dimensions 0,1,2 (corresponding to x,y,z) must be even. In addition LOCAL dimension 0 (x) has to be even.  Your lattice does not meet the GLOBAL requirement: latt_size[%d]=%d\n", 
	      mu, latt_size[mu]);
      
      exit(1);
    }
  }
  
  num = latt_size[0] / machine_size[0];
  if ( num % 2 != 0 )
    {
      fprintf(stderr,"This is a Dslash with checkerboarding in 3 dimensions. GLOBAL dimensions 0,1,2 (corresponding to x,y,z) must be even. In addition LOCAL dimension 0 (x) has to be even. Your lattice does not meet the LOCAL requirement: sublattice_size[0]=%d\n", 
	      num);
    QMP_abort(1);
  }

  make_shift_tables_3d(bound, 
		       getSiteCoords,
		       getLinearSiteIndex,
		       nodeNumber);

  /* Allocated space for the floating temps */
  /* Wasteful - allocate 3 times subgrid_vol_cb. Otherwise, need to pack the TAIL{1,2} halfspinor_buffer_offsets */
  nsize = 2*3*2*sizeof(double)*3*subgrid_vol_cb_3d*Nd3;  /* Note 3x4 half-ferm temps */
  if ((xchi1_3d = QMP_allocate_aligned_memory(nsize,64,0)) == 0)
  {
    QMP_error("init_wnxtsu3dslash: could not initialize xchi1_3d");
    QMP_abort(1);
  }
  if ((xchi2_3d = QMP_allocate_aligned_memory(nsize,64,0)) == 0)
  {
    QMP_error("init_wnxtsu3dslash: could not initialize xchi2_3d");
    QMP_abort(1);
  }
    
  chi1_3d = (halfspinor_array*)QMP_get_memory_pointer(xchi1_3d);
  chi2_3d = (halfspinor_array*)QMP_get_memory_pointer(xchi2_3d);
  
  /* Loop over all communicating directions and build up the two message
   * handles. If there is no communications, the message handles will not
   * be initialized 
   */
  num = 0;
    
  for(mu=0; mu < Nd3; ++mu) 
  {
    if(machine_size[mu] > 1) 
    {
      if (bound[0][0][mu] == 0)
      {
	QMP_error("init_sse_dslash: type 0 message size is 0");
	QMP_abort(1);
      }

      forw_msg_3d[num][0] = QMP_declare_msgmem(chi1_3d+subgrid_vol_cb_3d*(1+3*mu), bound[0][0][mu]*sizeof(halfspinor_array));
      forw_msg_3d[num][1] = QMP_declare_msgmem(chi1_3d+subgrid_vol_cb_3d*(2+3*mu), bound[0][0][mu]*sizeof(halfspinor_array));
      forw_mh_3d[num][0]  = QMP_declare_receive_relative(forw_msg_3d[num][1], mu, +1, 0);
      forw_mh_3d[num][1]  = QMP_declare_send_relative(forw_msg_3d[num][0], mu, -1, 0);
	
      if (bound[0][1][mu] == 0)
      {
	QMP_error("init_sse_dslash: type 0 message size is 0");
	QMP_abort(1);
      }

      back_msg_3d[num][0] = QMP_declare_msgmem(chi2_3d+subgrid_vol_cb_3d*(1+3*mu), bound[0][1][mu]*sizeof(halfspinor_array));
      back_msg_3d[num][1] = QMP_declare_msgmem(chi2_3d+subgrid_vol_cb_3d*(2+3*mu), bound[0][1][mu]*sizeof(halfspinor_array));
      back_mh_3d[num][0]  = QMP_declare_receive_relative(back_msg_3d[num][1], mu, -1, 0);
      back_mh_3d[num][1]  = QMP_declare_send_relative(back_msg_3d[num][0], mu, +1, 0);
	
      num++;
    }
  }

  if (num > 0) {
    forw_all_mh_3d = QMP_declare_multiple(&(forw_mh_3d[0][0]), 2*num);
    back_all_mh_3d = QMP_declare_multiple(&(back_mh_3d[0][0]), 2*num);
  }
  
  total_comm_3d = num;
  initP_3d = 1;
}

void free_sse_su3dslash_3d(void)
{
  const int *machine_size = QMP_get_logical_dimensions();
  int mu, num;

  /* If we are uninitialised just return */
  if (initP_3d == 0) {
    return;
  }

  /* Otherwise decrease the refcount */
  initP_3d--;

  /* If the refcount has now hit 0 then free stuff */
  if( initP_3d == 0 ) { 

    /* Free space */
    QMP_free_memory(xchi1_3d);
    QMP_free_memory(xchi2_3d);
    free_shift_tables();
    
    if (total_comm_3d > 0) {
      
      QMP_free_msghandle(forw_all_mh_3d);
      QMP_free_msghandle(back_all_mh_3d);
  
      num = 0;
      
      for(mu=0; mu < Nd3; ++mu) {
	
	if(machine_size[mu] > 1) {
	  QMP_free_msgmem(forw_msg_3d[num][0]);
	  QMP_free_msgmem(forw_msg_3d[num][1]);

	  QMP_free_msgmem(back_msg_3d[num][0]);
	  QMP_free_msgmem(back_msg_3d[num][1]);
	  
	  num++;
	}
      }
    }
    
    total_comm_3d = 0;
  }
}


void sse_su3dslash_wilson_3d(SSEREAL *u, SSEREAL *psi, SSEREAL *res, int isign, int cb)
{

  if (initP_3d == 0) {
    QMP_error("sse_su3dslash_wilson not initialised");
    QMP_abort(1);
  }

  if(isign==1) 
  {
    dispatch_to_threads(decomp_3d_plus,
			(spinor_array*)psi,
			chi1_3d,
			(my_mat_array)u,
			cb,
			subgrid_vol_cb_3d);

    if (total_comm_3d > 0)
      if (QMP_start(forw_all_mh_3d) != QMP_SUCCESS)
      {
	QMP_error("sse_su3dslash_wilson: QMP_start failed in forward direction");
	QMP_abort(1);
      }

    /*other checkerboard's worth */
    dispatch_to_threads(decomp_hvv_3d_plus,
			(spinor_array*)psi,
			chi2_3d,
			(my_mat_array)u,
			cb,
			subgrid_vol_cb_3d);
	
    if (total_comm_3d > 0)
      if (QMP_wait(forw_all_mh_3d) != QMP_SUCCESS)
      {
	QMP_error("sse_su3dslash_wilson: QMP_wait failed in forward direction");
	QMP_abort(1);
      }

   
    if (total_comm_3d > 0)
      if (QMP_start(back_all_mh_3d) != QMP_SUCCESS)
      {
	QMP_error("sse_su3dslash_wilson: QMP_start failed in backward direction");
	QMP_abort(1);
      }

    dispatch_to_threads(mvv_recons_3d_plus,
			(spinor_array*)res,
			chi1_3d,
			(my_mat_array)u,
			1-cb,
			subgrid_vol_cb_3d);
	
    if (total_comm_3d > 0)
      if (QMP_wait(back_all_mh_3d) != QMP_SUCCESS)
      {
	QMP_error("wnxtsu3dslash: QMP_wait failed in backward direction");
	QMP_abort(1);
      }


    dispatch_to_threads(recons_3d_plus,
			(spinor_array*)res, 
			chi2_3d,
			(my_mat_array)u,	
			1-cb,
			subgrid_vol_cb_3d);
  }		
  
  if(isign==-1) 
  {
    dispatch_to_threads(decomp_3d_minus,
			(spinor_array*)psi,
			chi1_3d,
			(my_mat_array)u,
			cb,
			subgrid_vol_cb_3d);
    
   
    if (total_comm_3d > 0)
      if (QMP_start(forw_all_mh_3d) != QMP_SUCCESS)
      {
	QMP_error("sse_su3dslash_wilson: QMP_start failed in forward direction");
	QMP_abort(1);
      }
	
	
    /*other checkerboard's worth */
    dispatch_to_threads(decomp_hvv_3d_minus,
			(spinor_array*)psi,
			chi2_3d,
			(my_mat_array)u,
			cb,
			subgrid_vol_cb_3d);
    
    
    if (total_comm_3d > 0)
      if (QMP_wait(forw_all_mh_3d) != QMP_SUCCESS)
      {
	QMP_error("sse_su3dslash_wilson: QMP_wait failed in forward direction");
	QMP_abort(1);
      }


    if (total_comm_3d > 0)
      if (QMP_start(back_all_mh_3d) != QMP_SUCCESS)
      {
	QMP_error("wnxtsu3dslash: QMP_start failed in backward direction");
	QMP_abort(1);
      }

    /*current cb's u */
    dispatch_to_threads(mvv_recons_3d_minus,
			(spinor_array*)res,
			chi1_3d,
			(my_mat_array)u,
			1-cb,
			subgrid_vol_cb_3d);
    

    if (total_comm_3d > 0)
      if (QMP_wait(back_all_mh_3d) != QMP_SUCCESS)
      {
	QMP_error("wnxtsu3dslash: QMP_wait failed in backward direction"); 
	QMP_abort(1);
      }

    dispatch_to_threads(recons_3d_minus,
			(spinor_array*)res, 
			chi2_3d,
			(my_mat_array)u,	
			1-cb,
			subgrid_vol_cb_3d);
  }		
}


#ifdef __cplusplus
}
#endif
