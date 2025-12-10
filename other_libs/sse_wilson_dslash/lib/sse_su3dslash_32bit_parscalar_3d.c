
#include <sse_config.h>            
#include <shift_tables_parscalar.h>
#include <stdlib.h>
#include <stdio.h>
#include <qmp.h>        /* QMP for the comms */
#include <sse_align.h>  /* Alignment stuff to ensure 16 byte alignments */
#include <types32.h>    /* Types and prefetch macros */
#include <dispatch_parscalar.h>   /* Threads dispatch definition */

#include "decomp.h"
#include "decomp_hvv.h"
#include "mvv_recons_32bit.h"
#include "recons.h" 
#ifdef __cplusplus
extern "C" {
#endif


  static int initP_3d=0;

  extern int subgrid_vol_3d;   /* Number of sites */
  extern int subgrid_vol_cb_3d; /* Number of sites in the checkerboard */



  /* Site table: Places an ordering on the sites. Separate into 2 checkerboards
     Then following a path through the lattice fill out the rest */

  extern int *site_table_3d;    /* Aligned lookup table */

  extern int *offset_table_body_3d;  /* Aligned offset table */

  extern int Nd3;

  static inline 
  int halfspinor_buffer_offset(HalfSpinorOffsetType type, int site, int mu)
  {
    int LocalNd3 = 3;
    return offset_table_body_3d[mu + LocalNd3*( site + subgrid_vol_3d*type) ];
  }

  



/****************************isign corresponding to +1  **************************/

/* the basic operation here is load up a spinor, do the spin decomposition, and store the halfspinor
to a lattice temporary */

  void decomp_3d_plus(size_t lo,size_t hi, int id, const void *ptr) /*need to fix decomp_minus */
  {
    int ix1;                           /* Site index - iz1 used at loop end */
    spinor_array* sp1 ALIGN;                /* Spinor under consideration */
       
    const ThreadWorkerArgs *a = (ThreadWorkerArgs *)ptr;          /* Cast the argument */
    halfspinor_array* chi = a->half_spinor; /* needs to be changed to halfspinor_array and be an array*/
    spinor_array* spinor_field = a->spinor;
    int cb = a->cb;


    halfspinor_array* s3;

    int low;
    int high;

    low = cb*subgrid_vol_cb_3d + lo;
    high = cb*subgrid_vol_cb_3d + hi;

    
    /************************ loop over all lattice sites *************************/
    for (ix1 = low; ix1 < high ; ix1++) {
      
      int thissite = site_table_3d[ix1];

      sp1=&spinor_field[ thissite ];
    
      /******************************* direction +0 *********************************/
      /* first of two sites */
      s3 = chi + halfspinor_buffer_offset(DECOMP_SCATTER, ix1, 0);

      decomp_gamma0_minus(sp1[0], *s3);

      /******************************* direction +1 *********************************/
      s3 = chi + halfspinor_buffer_offset(DECOMP_SCATTER, ix1, 1);

      decomp_gamma1_minus(sp1[0], *s3);
      
      /******************************* direction +2 *********************************/
      s3 = chi + halfspinor_buffer_offset(DECOMP_SCATTER, ix1, 2);

      decomp_gamma2_minus(sp1[0], *s3);

      //      s3 = chi + halfspinor_buffer_offset(DECOMP_SCATTER, ix1, 3);
      // decomp_gamma3_minus(sp1[0], *s3);
      
    }
  }


  /* the basic operations in this routine include loading a spinor, doing 
   * the spin projection, and multiplying the halfspinor by the appropriate 
   * gauge field, and saving the resulting halfspinor to a lattice temporary */
  
  /* need gauge fields on opposite cb */
  void decomp_hvv_3d_plus(size_t lo,size_t hi, int id, const void *ptr)
  {
    
    int ix1;                   /* Site addresses. ix1 = current. 
				  iz1 is for next loop iteration to allow some loop peeling
				  with ix1 */
    
    u_mat_array* um1 ALIGN;    /* Gauge pointer for 1 site */
 
    spinor_array* sm1 ALIGN;   /* spinor */


    const ThreadWorkerArgs *a = (const ThreadWorkerArgs *)ptr;
    spinor_array* spinor_field = a->spinor;
    halfspinor_array* chi = a->half_spinor; /* a 1-d map of a 2-d array */
    my_mat_array gauge_field = a->u;

    int cb = a->cb;

    halfspinor_array* s3;

    int low;
    int high;


    low = cb*subgrid_vol_cb_3d + lo;
    high = cb*subgrid_vol_cb_3d + hi;


    
    /************************ loop over all lattice sites *************************/
    for (ix1 =low; ix1 <high ;ix1++) {

      int thissite = site_table_3d[ ix1 ];

      sm1=&spinor_field[ thissite ];
      um1 = &gauge_field[ thissite ][0];

      s3 = chi + halfspinor_buffer_offset(DECOMP_HVV_SCATTER, ix1, 0);    
      decomp_hvv_gamma0_plus(*sm1,*um1,*s3);


      /******************************* direction +1 *********************************/
      um1 = &gauge_field[ thissite ][1];
      s3 = chi + halfspinor_buffer_offset(DECOMP_HVV_SCATTER, ix1, 1);
      decomp_hvv_gamma1_plus(*sm1,*um1,*s3);
      
      /******************************* direction +2 *********************************/
      um1 = &gauge_field[ thissite ][2];
      s3 = chi + halfspinor_buffer_offset(DECOMP_HVV_SCATTER, ix1, 2);
      decomp_hvv_gamma2_plus(*sm1,*um1,*s3);

      /******************************* direction +3 *********************************/
      // um1 = &gauge_field[ thissite ][3];
      //      s3 = chi + halfspinor_buffer_offset(DECOMP_HVV_SCATTER, ix1, 3);
      // decomp_hvv_gamma3_plus(*sm1,*um1,*s3);

      /******************************** end of loop *********************************/
  }
}
/***************end of decomp_hvv****************/


/* the basic operations in this routine include loading the halfspinor 
 * from memory, multiplying it by the appropriate gauge field, doing the 
 * spin reconstruction, and summing over directions, and saving the partial 
 * sum over directions */

void mvv_recons_3d_plus(size_t lo,size_t hi, int id, const void *ptr)
{

  int ix1;
  int low;
  int high;

    
  u_mat_array* up1 ALIGN;
  spinor_array* sn1 ALIGN;
  halfspinor_array r12_1 ALIGN, r34_1 ALIGN;

  const ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr;
  spinor_array* spinor_field = a->spinor;
  halfspinor_array* chi = a->half_spinor; /* a 1-d map of a 2-d array */

  my_mat_array gauge_field = a->u;
  int cb = a->cb;

  halfspinor_array* s3;
  

  low = cb*subgrid_vol_cb_3d + lo;
  high = cb*subgrid_vol_cb_3d + hi;

  /************************ loop over all lattice sites *************************/
  for (ix1 =low; ix1 <high ;ix1++) {

    //    int thissite = lookup_site(cb,ix1);
    int thissite = site_table_3d[ ix1 ];
  
    up1=&gauge_field[thissite][0];
    s3 = chi + halfspinor_buffer_offset(RECONS_MVV_GATHER, ix1, 0);
    
    mvv_recons_gamma0_plus(*s3, *up1, r12_1, r34_1);


    /***************************** direction +1 ***********************************/
    up1=&gauge_field[thissite][1];
    s3 = chi + halfspinor_buffer_offset(RECONS_MVV_GATHER, ix1, 1);
    mvv_recons_gamma1_plus_add(*s3, *up1, r12_1, r34_1);


    /******************************* direction +2 *********************************/
    up1=&gauge_field[thissite][2];
    s3 = chi + halfspinor_buffer_offset(RECONS_MVV_GATHER, ix1, 2); 
    sn1=&spinor_field[thissite];
    mvv_recons_gamma2_plus_add_store(*s3, *up1, r12_1, r34_1, *sn1);

    /******************************* direction +2 *********************************/
    //    up1=&gauge_field[thissite][3];
    // s3 = chi + halfspinor_buffer_offset(RECONS_MVV_GATHER, ix1, 3); 
    // sn1=&spinor_field[thissite];
    // mvv_recons_gamma3_plus_add_store(*s3, *up1, r12_1, r34_1, *sn1);


  }
}


   
/* this routine takes the partial sum from mvv_recons() and loops 
 * over the output spin components, 2 at a time doing a sum over directions 
 * for each set, accumulating in xmm0-2 and loading the halfspinor 
 * temporaries into xmm3-5 */

void recons_3d_plus(size_t lo,size_t hi, int id, const void *ptr )	
{
  int ix1;
  spinor_array* sn1 ALIGN;
  

  const ThreadWorkerArgs *a = (ThreadWorkerArgs *)ptr;
  spinor_array* spinor_field = a->spinor;
  halfspinor_array* chi = a->half_spinor;
  int cb = a->cb;


 
  halfspinor_array *hs0, *hs1, *hs2;  
  int low;
  int high;

  low = cb*subgrid_vol_cb_3d + lo;
  high = cb*subgrid_vol_cb_3d + hi;
  
  
  /************************ loop over all lattice sites *************************/
  for (ix1 =low; ix1 <high ;ix1++) {

    //    int thissite=lookup_site(cb,ix1);
    int thissite = site_table_3d[ ix1 ];

    hs0 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1,0); 
    hs1 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1,1);
    hs2 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1,2);
    //    hs3 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1,3);
    sn1=&spinor_field[thissite];   
    recons_3dir_plus(*hs0, *hs1, *hs2, *sn1);
  
    /*************************end of loop ****************************/
  }
}
/*****************end of recons**************/





/*************** now for isign corresponding to -1  ****************************************/

void decomp_3d_minus(size_t lo,size_t hi, int id, const void *ptr ) /*need to fix decomp_minus */
{

  int ix1;
  spinor_array* sp1 ALIGN;

  const ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr;
  halfspinor_array* chi = a->half_spinor; /* needs to be changed to halfspinor_array and be an array*/

  halfspinor_array* s3;

  int cb = a->cb;
  spinor_array* spinor_field= a->spinor;
   
  /************************ loop over all lattice sites *************************/
  int low = cb*subgrid_vol_cb_3d + lo;
  int high = cb*subgrid_vol_cb_3d + hi;
  
  
  /************************ loop over all lattice sites *************************/
  for (ix1 =low; ix1 <high ;ix1++) {

    int thissite = site_table_3d[ ix1 ];       
    sp1=&spinor_field[thissite];

   /******************************* direction +0 *********************************/
    /* ...(1-gamma(0))... */
    s3 = chi + halfspinor_buffer_offset(DECOMP_SCATTER,ix1,0);
    decomp_gamma0_plus(sp1[0], *s3);

    /******************************* direction +1 *********************************/
    s3 = chi + halfspinor_buffer_offset(DECOMP_SCATTER,ix1,1);    
    decomp_gamma1_plus(sp1[0], *s3);

    /******************************* direction +2 *********************************/
    s3 = chi + halfspinor_buffer_offset(DECOMP_SCATTER,ix1,2);
    decomp_gamma2_plus(sp1[0], *s3);

    /******************************* direction +2 *********************************/
    //    s3 = chi + halfspinor_buffer_offset(DECOMP_SCATTER,ix1,3);
    //    decomp_gamma3_plus(sp1[0], *s3);

 
  }
}


/* need gauge fields on opposite cb */
void decomp_hvv_3d_minus(size_t lo,size_t hi, int id, const void *ptr )
{

  int ix1;
  u_mat_array* um1 ALIGN;
  spinor_array* sm1 ALIGN;


  const ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr;
  
  spinor_array* spinor_field = a->spinor;

  halfspinor_array* chi = a->half_spinor; /* a 1-d map of a 2-d array */
  
  my_mat_array gauge_field = a->u;
  
  int cb = a->cb;
  
  halfspinor_array* s3;


  /************************ loop over all lattice sites *************************/
  int low = cb*subgrid_vol_cb_3d + lo;
  int high = cb*subgrid_vol_cb_3d + hi;
  
  
  /************************ loop over all lattice sites *************************/
  for (ix1 =low; ix1 <high ;ix1++) {

    int thissite = site_table_3d[ ix1 ];       

    um1=&gauge_field[thissite][0]; 
    sm1=&spinor_field[thissite];

    s3 = chi + halfspinor_buffer_offset(DECOMP_HVV_SCATTER,ix1,0);
    decomp_hvv_gamma0_minus(*sm1, *um1, *s3);


    /******************************* direction -1 *********************************/
    um1=&gauge_field[thissite][1]; 
    s3 = chi + halfspinor_buffer_offset(DECOMP_HVV_SCATTER,ix1,1);
    decomp_hvv_gamma1_minus(*sm1, *um1, *s3);


    /******************************* direction -2 *********************************/    
    um1=&gauge_field[thissite][2]; 
    s3 = chi + halfspinor_buffer_offset(DECOMP_HVV_SCATTER,ix1,2);
    decomp_hvv_gamma2_minus(*sm1, *um1, *s3);

    /******************************* direction -2 *********************************/    
    // um1=&gauge_field[thissite][3]; 
    // s3 = chi + halfspinor_buffer_offset(DECOMP_HVV_SCATTER,ix1,3);
    // decomp_hvv_gamma3_minus(*sm1, *um1, *s3);

  }
}


void mvv_recons_3d_minus(size_t lo,size_t hi, int id, const void *ptr )
{
  int ix1;
  u_mat_array* up1 ALIGN;   /* Gauge pointers for site x and x+1 */
  spinor_array* sn1 ALIGN;  /* The spinor to store to */


  /* Temporaries for the top and bottom parts of spinors. */
  halfspinor_array r12_1 ALIGN,r34_1 ALIGN;

  /* if going to support unpacked gauge fields, need to treat site ix1 and site ix1+1 separately */
  /* to support unpacked gauge fields the prefetches will need to be changed */
  const ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr;
  spinor_array* spinor_field = a->spinor;
  halfspinor_array* chi = a->half_spinor; /* a 1-d map of a 2-d array */
  my_mat_array gauge_field = a->u;
  int cb = a->cb;

  halfspinor_array* s3;
  

  /************************ loop over all lattice sites *************************/
  int low = cb*subgrid_vol_cb_3d + lo;
  int high = cb*subgrid_vol_cb_3d + hi;
  
  
  /************************ loop over all lattice sites *************************/
  for (ix1 =low; ix1 <high ;ix1++) {

    int thissite = site_table_3d[ ix1 ];       

    up1=&gauge_field[thissite][0];
    /******************************* direction +0 *********************************/
    s3 = chi + halfspinor_buffer_offset(RECONS_MVV_GATHER,ix1,0);
    mvv_recons_gamma0_minus(*s3, *up1, r12_1, r34_1);


    /******************************* direction +1 *********************************/
    up1=&gauge_field[thissite][1];
    s3 = chi + halfspinor_buffer_offset(RECONS_MVV_GATHER,ix1,1);
    mvv_recons_gamma1_minus_add(*s3, *up1, r12_1, r34_1);

    /******************************* direction +2 *********************************/
    up1=&gauge_field[thissite][2];
    s3 = chi + halfspinor_buffer_offset(RECONS_MVV_GATHER,ix1,2);
    sn1=&spinor_field[thissite];
    mvv_recons_gamma2_minus_add_store(*s3, *up1, r12_1, r34_1, *sn1);

    //    up1=&gauge_field[thissite][3];
    // s3 = chi + halfspinor_buffer_offset(RECONS_MVV_GATHER,ix1,3);
    // mvv_recons_gamma3_minus_add_store(*s3, *up1, r12_1, r34_1, *sn1);

    /******************************** end of loop *********************************/
  }
}
/******************end of mvv_recons*************************/


void recons_3d_minus(size_t lo,size_t hi, int id, const void *ptr )	
{
  int ix1;
  spinor_array* sn1 ALIGN;


  const ThreadWorkerArgs *a = (ThreadWorkerArgs *)ptr;
  spinor_array* spinor_field = a->spinor;
  halfspinor_array* chi = a->half_spinor; /* a 1-d map of a 2-d array */
  int cb = a->cb;
  halfspinor_array* hs0, *hs1, *hs2;

  /************************ loop over all lattice sites *************************/
  int low = cb*subgrid_vol_cb_3d + lo;
  int high = cb*subgrid_vol_cb_3d + hi;
  
  
  /************************ loop over all lattice sites *************************/
  for (ix1 =low; ix1 <high ;ix1++) {

    int thissite = site_table_3d[ ix1 ];       

    sn1=&spinor_field[thissite];   

    hs0 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1,0); 
    hs1 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1,1);
    hs2 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1,2);
    //    hs3 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1,3);

    recons_3dir_minus(*hs0, *hs1, *hs2, *sn1);
  }

}
/*****************end of isign corresponding to -1 **************************************/


/***************** start of initialization routine ***************************************/


static QMP_mem_t* xchi1_3d;               /* QMP Memory Structures for halfspinor arrays */
static QMP_mem_t* xchi2_3d;               /* xchi1 => FORWARD, xchi2 => BACKWARD         */

static halfspinor_array* chi1_3d;         /* These are the aligned pointers from the QMP Memory structures */
static halfspinor_array* chi2_3d;         /* xchi1 <=> chi1    xchi2 <=> chi2 */



/* Nearest neighbor communication channels */
static int total_comm_3d = 0;
static QMP_msgmem_t forw_msg_3d[3][2];
static QMP_msgmem_t back_msg_3d[3][2];
static QMP_msghandle_t forw_mh_3d_send[3];
static QMP_msghandle_t forw_mh_3d_recv[3];
static QMP_msghandle_t back_mh_3d_send[3];
static QMP_msghandle_t back_mh_3d_recv[3];

#define DSLASH_COLLAPSE_MESSAGES
#ifdef DSLASH_COLLAPSE_MESSAGES
#warning Using Collapsed QMP Message Handles
static QMP_msghandle_t forw_all_mh_3d_send;
static QMP_msghandle_t forw_all_mh_3d_recv;
static QMP_msghandle_t back_all_mh_3d_send;
static QMP_msghandle_t back_all_mh_3d_recv;
#else
#warning Using Individual Message Handles
#endif

/* Initialize the Dslash */
void init_sse_su3dslash_3d(const int latt_size[],
			     void (*getSiteCoords)(int coord[], int node, int linearsite),
			    
			     int (*getLinearSiteIndex)(const int coord[]),
			   int (*nodeNumber)(const int coord[]))
{


  const int *machine_size = QMP_get_logical_dimensions();

  /* bound[cb][forward=0/backward=1][dir] */
  int bound[2][4][3];

  int mu, num, nsize;


  /* If we are already initialised, then increase the refcount and return */
  if (initP_3d > 0) 
  {
    initP_3d++;
    return;
  }


  /* Otherwise initialise */
  /* Check we are in 4D */
  if (QMP_get_logical_number_of_dimensions() != 4) {
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


  /* Make the shift table  -- this sets the vol and vol_cb so we can call getSugridVolCB() after it */
  make_shift_tables_3d(bound,
		       getSiteCoords,
		       getLinearSiteIndex,
		       nodeNumber);


  /* Is Icolor relevant here? I think we'll need to deal with explicit indices */
  /*
  icolor_end[0] = icolor_start[0] + getSubgridVolCB();
  icolor_end[1] = icolor_start[1] + getSubgridVolCB();
  */

  /* Allocate the halfspinor temporaries */
  /* The halfspinor is really the following:
     
        halfspinor[0][Nd-1][vol_cb] - half spinors from the body
	halfspinor[1][Nd-1][vol_cb] - half spinors from the tail to be sent forward
	halfspinor[2][Nd-1][vol_cb] - half spinors from the tail to be sent backward

	The Nd-1 stands for the 3 decomposition directions (gamma_matrices)

	The body and tails are overallocated to their maximum size: vol_cb
	Hence the size is half_spinor_size * vol_cb*3*(Nd-1)
	where Nd=4, 3 is the number of buffers (body, forward tail, backward tail)
	and   half spinor size is 2spin comp * 3 color * 2 floats (=1complex)

  */
  
  /* 2x3x2xsizeof(float) = halfspinor. (Result of decomp)
     3 buffers per direction (body, send, receive 
     3 directions
  */
  nsize = 2*3*2*sizeof(float)*3*subgrid_vol_cb_3d*Nd3;  /* Note 3x3 half-fermions */

 
  /* xchi1 and xchi2 will hold projected spinors memory handles from QMP */
  /* Align it to be 64 bytes aligned - ie cache line size aligned */
  if ((xchi1_3d = QMP_allocate_aligned_memory(nsize,64,0)) == 0)
  {
    fprintf(stderr, "nsize=%d subgrid_vol_cb3d=%d\n", nsize, subgrid_vol_cb_3d);
    fflush(stderr);
    QMP_error("init_wnxtsu3dslash: could not initialize xchi1");
    QMP_abort(1);
  }


  if((xchi2_3d = QMP_allocate_aligned_memory(nsize,64,0)) == 0)
  {
    QMP_error("init_wnxtsu3dslash: could not initialize xchi2");
    QMP_abort(1);
  }


  /* Unwrap the half spinor pointers from the QMP Structures. BTW: This 2 step technique susks so bad! */
  chi1_3d = (halfspinor_array*)QMP_get_memory_pointer(xchi1_3d);

  chi2_3d = (halfspinor_array*)QMP_get_memory_pointer(xchi2_3d); 


  /* Loop over all communicating directions and build up the two message
   * handles. If there is no communications, the message handles will not
   * be initialized 
   */
  num = 0;

  /* Loop over directions */
  for(mu=0; mu < Nd3; ++mu) {

    if(machine_size[mu] > 1) { /* If the machine is not a scalar  in this dimensio */

      if (bound[0][0][mu] == 0) { /* Consistency: Check the boundary in this direction is 0 */
 	QMP_error("init_sse_dslash: type 0 message size is 0 %d %d", mu, bound[0][0][mu]);
	QMP_abort(1);
      }


      /* 
	 Boundary is indexed as            boundary[cb][forward=0/back=1][direction] 
      */	 

      /* cb = 0 */
      forw_msg_3d[num][0] = QMP_declare_msgmem( chi1_3d+subgrid_vol_cb_3d*(1+3*mu),
					     bound[0][0][mu]*sizeof(halfspinor_array));

      /* cb = 1 */
      forw_msg_3d[num][1] = QMP_declare_msgmem( chi1_3d +subgrid_vol_cb_3d*(2+3*mu),
					     bound[1][0][mu]*sizeof(halfspinor_array));

      /* cb = 0: Receive from cb = 1 */
      forw_mh_3d_recv[num]  = QMP_declare_receive_relative(forw_msg_3d[num][1], mu, +1, 0);

      /* cb = 1: send to cb = 0 */
      forw_mh_3d_send[num]  = QMP_declare_send_relative(forw_msg_3d[num][0], mu, -1, 0);
	
      if (bound[0][1][mu] == 0)
      {
	QMP_error("init_sse_dslash: type 1 message size is 0, %d, %d",mu, bound[0][1][mu]);
	QMP_abort(1);
      }

      /* cb = 0 */
      back_msg_3d[num][0] = QMP_declare_msgmem( chi2_3d+subgrid_vol_cb_3d*(1+3*mu),
					     bound[0][1][mu]*sizeof(halfspinor_array));

      /* cb = 1 */
      back_msg_3d[num][1] = QMP_declare_msgmem( chi2_3d+subgrid_vol_cb_3d*(2+3*mu), 
					     bound[1][1][mu]*sizeof(halfspinor_array));

      /* cb = 0: Receive from cb=1 */
      back_mh_3d_recv[num]  = QMP_declare_receive_relative(back_msg_3d[num][1], mu, -1, 0);

      /* cb = 1: Receive from cb=0 */
      back_mh_3d_send[num]  = QMP_declare_send_relative(back_msg_3d[num][0], mu, +1, 0);


	
      num++;
    }
  }

  /* Combine the messages */

  /* Dont combine messages */

#ifdef DSLASH_COLLAPSE_MESSAGES
  if (num > 0) {
    forw_all_mh_3d_send = QMP_declare_multiple(&(forw_mh_3d_send[0]), num);
    forw_all_mh_3d_recv = QMP_declare_multiple(&(forw_mh_3d_recv[0]), num);
    back_all_mh_3d_send = QMP_declare_multiple(&(back_mh_3d_send[0]), num);
    back_all_mh_3d_recv = QMP_declare_multiple(&(back_mh_3d_recv[0]), num);
  }
#endif

  total_comm_3d = num;
  initP_3d = 1;

}


void free_sse_su3dslash_3d(void)
{
  const int *machine_size = QMP_get_logical_dimensions();
  int mu, num;

  /***************HACK********************/
  /* There is a gigE QMP_free bug  (2/6/04). For now, turn off ever free-ing !! */
  return;


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

    /* Free shift table */
    free_shift_tables();
    

    if (total_comm_3d > 0) {
#ifdef DSLASH_COLLAPSE_MESSAGES
      QMP_free_msghandle(forw_all_mh_3d_send);
      QMP_free_msghandle(forw_all_mh_3d_recv);
      QMP_free_msghandle(back_all_mh_3d_send);
      QMP_free_msghandle(back_all_mh_3d_recv);
#else
      for(num=0; num < total_comm_3d; num++) { 
        QMP_free_msghandle(forw_mh_3d_send[num]);
        QMP_free_msghandle(forw_mh_3d_recv[num]);
        QMP_free_msghandle(back_mh_3d_send[num]);
        QMP_free_msghandle(back_mh_3d_recv[num]);
      }
#endif

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

/***************** end of initialization routine ***************************************/


void sse_su3dslash_wilson_3d(float *u, float *psi, float *res, int isign, int cb)
{
  int i; 

  if (initP_3d == 0) {
    QMP_error("sse_su3dslash_wilson not initialized");
    QMP_abort(1);
  }

  if(isign==1) 
  {
    /* Post all receives */
    if (total_comm_3d > 0 ) { 
#ifdef DSLASH_COLLAPSE_MESSAGES
      if(QMP_start(forw_all_mh_3d_recv) != QMP_SUCCESS) { 
	QMP_error("sse_su3dslash_wilson: QMP_start(forw_all_mh_3d_recv)");
	QMP_abort(1);
      }

      if(QMP_start(back_all_mh_3d_recv) != QMP_SUCCESS) { 
	QMP_error("sse_su3dslash_wilson: QMP_start(back_all_mh_3d_recv");
	QMP_abort(1);
      }
#else
      for(i=0; i < total_comm_3d; i++) { 
	if (QMP_start(forw_mh_3d_recv[i]) != QMP_SUCCESS) {
	  QMP_error("sse_su3dslash_wilson: QMP_start failed in forward direction");
	  QMP_abort(1);
	}

	if (QMP_start(back_mh_3d_recv[i]) != QMP_SUCCESS) {
	  QMP_error("sse_su3dslash_wilson: QMP_start failed in backward direction");
	  QMP_abort(1);
	}
      }
#endif
    }

    dispatch_to_threads(decomp_3d_plus,
		(spinor_array*)psi,
		chi1_3d,
		(my_mat_array)u,
		cb,
		subgrid_vol_cb_3d);


    /* Start send forward */
    if (total_comm_3d > 0) {
#ifdef DSLASH_COLLAPSE_MESSAGES
      if( QMP_start(forw_all_mh_3d_send) != QMP_SUCCESS) { 
	QMP_error("sse_su3dslash_wilson_3d: QMP_start(forw_all_mh_3d) failed");
	QMP_abort(1);
      }
#else
      for(i=0; i < total_comm_3d;i++) { 
	if (QMP_start(forw_mh_3d_send[i]) != QMP_SUCCESS) {
      
	  QMP_error("sse_su3dslash_wilson: QMP_start failed in forward direction");
	  QMP_abort(1);
	}
      }
#endif
    }



    dispatch_to_threads(decomp_hvv_3d_plus,
		(spinor_array*)psi,
		chi2_3d,
		(my_mat_array)u,
		cb,
		subgrid_vol_cb_3d);
	


    if (total_comm_3d > 0) {
#ifdef DSLASH_COLLAPSE_MESSAGES
      if( QMP_wait(forw_all_mh_3d_send) != QMP_SUCCESS) {
	QMP_error("sse_su3ddslash_wilson: QMP_wait(forw_all_mh_3d_send)");
	QMP_abort(1);
      }
      if( QMP_start(back_all_mh_3d_send) != QMP_SUCCESS) {
	QMP_error("sse_su3ddslash_wilson: QMP_start(back_all_mh_3d_send)");
	QMP_abort(1);
      }
      if( QMP_wait(forw_all_mh_3d_recv) != QMP_SUCCESS) {
	QMP_error("sse_su3ddslash_wilson: QMP_wait(forw_all_mh_3d_recv)");
	QMP_abort(1);
      }
#else
      /* Finish fwd sends */
      for(i=0; i < total_comm_3d; i++) { 
	if (QMP_wait(forw_mh_3d_send[i]) != QMP_SUCCESS) {
	  QMP_error("sse_su3dslash_wilson: QMP_wait failed in forward direction");
	  QMP_abort(1);
	}

	/* Start bwd sends */
	if (QMP_start(back_mh_3d_send[i]) != QMP_SUCCESS) {
	  QMP_error("sse_su3dslash_wilson: QMP_start failed in backward direction");
	  QMP_abort(1);
	}
      
	/* Finish all forward receives */
	if (QMP_wait(forw_mh_3d_recv[i]) != QMP_SUCCESS) {
	  QMP_error("sse_su3dslash_wilson: QMP_wait failed in forward direction");
	  QMP_abort(1);
	}
      }
#endif
    }

    dispatch_to_threads(mvv_recons_3d_plus,
		(spinor_array*)res,
		chi1_3d,
		(my_mat_array)u,
		1-cb,
		subgrid_vol_cb_3d);



    if (total_comm_3d > 0) { 
#ifdef DSLASH_COLLAPSE_MESSAGES
      if( QMP_wait(back_all_mh_3d_send) != QMP_SUCCESS) {
	QMP_error("sse_su3ddslash_wilson: QMP_wait(back_all_mh_3d_send)");
	QMP_abort(1);
      }
      if( QMP_wait(back_all_mh_3d_recv) != QMP_SUCCESS) {
	QMP_error("sse_su3ddslash_wilson: QMP_start(back_all_mh_3d_recv)");
	QMP_abort(1);
      }
#else
      for(i=0; i < total_comm_3d; i++) { 
	if (QMP_wait(back_mh_3d_send[i]) != QMP_SUCCESS) {

	  QMP_error("wnxtsu3dslash: QMP_wait failed in backward direction");
	  QMP_abort(1);
	}

	if (QMP_wait(back_mh_3d_recv[i]) != QMP_SUCCESS) {
	  
	  QMP_error("wnxtsu3dslash: QMP_wait failed in backward direction");
	  QMP_abort(1);
	}
      }
#endif
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


    /* Start all receives */
    if( total_comm_3d > 0 ) { 
#ifdef DSLASH_COLLAPSE_MESSAGES
      if(QMP_start(forw_all_mh_3d_recv) != QMP_SUCCESS) { 
	QMP_error("sse_su3dslash_wilson: QMP_start(forw_all_mh_3d_recv)");
	QMP_abort(1);
      }

      if(QMP_start(back_all_mh_3d_recv) != QMP_SUCCESS) { 
	QMP_error("sse_su3dslash_wilson: QMP_start(back_all_mh_3d_recv");
	QMP_abort(1);
      }
#else
      for(i=0; i < total_comm_3d;i++) { 
	if (QMP_start(forw_mh_3d_recv[i]) != QMP_SUCCESS) {
	  QMP_error("sse_su3dslash_wilson: QMP_start failed in forward direction");
	  QMP_abort(1);
	}
	if (QMP_start(back_mh_3d_recv[i]) != QMP_SUCCESS) {
	  QMP_error("sse_su3dslash_wilson: QMP_start failed in backward direction");
	  QMP_abort(1);
	}
      }
#endif
    }

    dispatch_to_threads(decomp_3d_minus,
		(spinor_array*)psi,
		chi1_3d,
		(my_mat_array)u,
		cb,
		subgrid_vol_cb_3d);


    if (total_comm_3d > 0) {
#ifdef DSLASH_COLLAPSE_MESSAGES
      if( QMP_start(forw_all_mh_3d_send) != QMP_SUCCESS) { 
	QMP_error("sse_su3dslash_wilson_3d: QMP_start(forw_all_mh_3d) failed");
	QMP_abort(1);
      }
#else
      for(i=0; i < total_comm_3d; i++) { 
	/* Start forward send */
	if (QMP_start(forw_mh_3d_send[i]) != QMP_SUCCESS) {
	  QMP_error("sse_su3dslash_wilson: QMP_start failed in forward direction");
	  QMP_abort(1);
	}
      }
#endif
    }

    dispatch_to_threads(decomp_hvv_3d_minus,
		(spinor_array*)psi,
		chi2_3d,
		(my_mat_array)u,
		cb,
		subgrid_vol_cb_3d);

    if (total_comm_3d > 0) {
#ifdef DSLASH_COLLAPSE_MESSAGES
      if( QMP_wait(forw_all_mh_3d_send) != QMP_SUCCESS) {
	QMP_error("sse_su3ddslash_wilson: QMP_wait(forw_all_mh_3d_send)");
	QMP_abort(1);
      }
      if( QMP_start(back_all_mh_3d_send) != QMP_SUCCESS) {
	QMP_error("sse_su3ddslash_wilson: QMP_start(back_all_mh_3d_send)");
	QMP_abort(1);
      }
      if( QMP_wait(forw_all_mh_3d_recv) != QMP_SUCCESS) {
	QMP_error("sse_su3ddslash_wilson: QMP_wait(forw_all_mh_3d_recv)");
	QMP_abort(1);
      }
#else
      for(i=0; i < total_comm_3d; i++) { 
	/* Finish forward sends */
	if (QMP_wait(forw_mh_3d_send[i]) != QMP_SUCCESS) {
      
	  QMP_error("sse_su3dslash_wilson: QMP_wait failed in forward direction");
	  QMP_abort(1);
	}
      
	/* Start backward sends */
	if (QMP_start(back_mh_3d_send[i]) != QMP_SUCCESS) {
	  QMP_error("sse_su3dslash_wilson: QMP_start failed in backward direction");
	  QMP_abort(1);
	}

	/* Finish forward receives */
	if (QMP_wait(forw_mh_3d_recv[i]) != QMP_SUCCESS) {
	  QMP_error("sse_su3dslash_wilson: QMP_wait failed in forward direction");
	  QMP_abort(1);
	}
      }
#endif
    }

    dispatch_to_threads(mvv_recons_3d_minus,
		(spinor_array*)res,
		chi1_3d,
		(my_mat_array)u,
		1-cb,
		subgrid_vol_cb_3d);
    


    if (total_comm_3d > 0) {
#ifdef DSLASH_COLLAPSE_MESSAGES
      if( QMP_wait(back_all_mh_3d_send) != QMP_SUCCESS) {
	QMP_error("sse_su3ddslash_wilson: QMP_wait(back_all_mh_3d_send)");
	QMP_abort(1);
      }
      if( QMP_wait(back_all_mh_3d_recv) != QMP_SUCCESS) {
	QMP_error("sse_su3ddslash_wilson: QMP_start(back_all_mh_3d_recv)");
	QMP_abort(1);
      }
#else
      /* Finish backward sends */
      for(i=0; i < total_comm_3d; i++) { 
	if (QMP_wait(back_mh_3d_send[i]) != QMP_SUCCESS) {
	  QMP_error("sse_su3dslash_wilson: QMP_wait failed in backward direction");
	  QMP_abort(1);
	}

	/* Finish backwards receives */
	if (QMP_wait(back_mh_3d_recv[i]) != QMP_SUCCESS) {
	  
	  QMP_error("sse_su3dslash_wilson: QMP_wait failed in backward direction");
	  QMP_abort(1);
	}
      }
#endif
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

