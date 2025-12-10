

#include "cpp_dslash_types.h"
#include "shift_table_parscalar.h"
#include "dispatch_parscalar.h"
#include "cpp_dslash_parscalar_utils_64bit.h"



namespace CPlusPlusWilsonDslash {
  using namespace Dslash64BitTypes;

  namespace DslashParscalar64Bit {
  
 void decomp_plus(size_t lo,size_t hi, int id, const void *ptr) /*need to fix decomp_minus */
{
  ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr; 
  ShiftTable<HalfSpinor>* tab = (ShiftTable<HalfSpinor>*)a->s;
  int subgrid_vol_cb = tab->subgridVolCB();

  int cb = a->cb; 
  FourSpinor *psi = (FourSpinor*)a->spinor;
  HalfSpinor *chi = (HalfSpinor*)a->half_spinor;
  
  int ix1=0;
  const int low = cb*subgrid_vol_cb + lo; 
  const int high = cb*subgrid_vol_cb + hi;


  HalfSpinor *s3;
  FourSpinor* sp ALIGN;

  for (ix1=low;ix1<high;ix1++) {
    int thissite = tab->siteTable(ix1);

    sp=&psi[thissite];
    /******************************* direction +0 *********************************/	   
    s3 = tab->halfspinorBufferOffset(DECOMP_SCATTER,ix1,0);
    decomp_gamma0_minus(*sp, *s3);
    /******************************* direction +1 *********************************/
    s3 = tab->halfspinorBufferOffset(DECOMP_SCATTER,ix1,1);
    decomp_gamma1_minus(*sp, *s3);
    
    /******************************* direction +2 *********************************/
    s3 = tab->halfspinorBufferOffset(DECOMP_SCATTER,ix1,2);
    decomp_gamma2_minus(*sp, *s3);

    /******************************* direction +3 *********************************/
    s3 = tab->halfspinorBufferOffset(DECOMP_SCATTER,ix1,3);
    decomp_gamma3_minus(*sp, *s3);

  }
}


/* the basic operations in this routine include loading a spinor, doing 
 * the spin projection, and multiplying the halfspinor by the appropriate 
 * gauge field, and saving the resulting halfspinor to a lattice temporary */

/* need gauge fields on opposite cb */
 void decomp_hvv_plus(size_t lo,size_t hi, int id, const void *ptr)
{
  ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr; 
  ShiftTable<HalfSpinor>* tab = (ShiftTable<HalfSpinor>*)a->s;
  int subgrid_vol_cb = tab->subgridVolCB();

  int cb = a->cb; 
  GaugeMatrix (*gauge_field)[4] = (GaugeMatrix(*)[4])a->u;
  FourSpinor *psi = (FourSpinor*)a->spinor;
  HalfSpinor *chi = (HalfSpinor*)a->half_spinor;
  int ix1=0;
  const int low = cb*subgrid_vol_cb + lo; 
  const int high = cb*subgrid_vol_cb + hi;

  GaugeMatrix *um ALIGN;

  HalfSpinor *s3 ALIGN;

  FourSpinor *sm ALIGN; 

  /*  printf("ID: %d low=%d high=%d: DecompHvvPlus\n", id, low, high);*/

  for (ix1=low;ix1<high;ix1++) 
  {
    int thissite = tab->siteTable(ix1);

    /* Spinor to project*/
    sm=&psi[thissite];

    /******************************* direction -1 *********************************/
    um=&gauge_field[thissite][0];
    s3 = tab->halfspinorBufferOffset(DECOMP_HVV_SCATTER,ix1,0);
    decomp_hvv_gamma0_plus(*sm, *um, *s3);
    
    /******************************* direction -1 *********************************/
    um=&gauge_field[thissite][1];
    s3 = tab->halfspinorBufferOffset(DECOMP_HVV_SCATTER,ix1,1);
    decomp_hvv_gamma1_plus(*sm, *um, *s3);

    /******************************* direction -2 *********************************/
    um=&gauge_field[thissite][2];
    s3 = tab->halfspinorBufferOffset(DECOMP_HVV_SCATTER,ix1,2);
    decomp_hvv_gamma2_plus(*sm, *um, *s3);

    /******************************* direction -3 *********************************/
    um=&gauge_field[thissite][3];
    s3 = tab->halfspinorBufferOffset(DECOMP_HVV_SCATTER,ix1,3);
    decomp_hvv_gamma3_plus(*sm, *um, *s3);
  }
}
/***************end of decomp_hvv****************/


/* the basic operations in this routine include loading the halfspinor 
 * from memory, multiplying it by the appropriate gauge field, doing the 
 * spin reconstruction, and summing over directions, and saving the partial 
 * sum over directions */

 void mvv_recons_plus(size_t lo,size_t hi, int id, const void *ptr)
{
  ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr;
  ShiftTable<HalfSpinor>* tab = (ShiftTable<HalfSpinor>*)a->s;
  int subgrid_vol_cb = tab->subgridVolCB();
 
  int cb = a->cb; 

  int ix1=0;
  const int low = cb*subgrid_vol_cb + lo; 
  const int high = cb*subgrid_vol_cb + hi;

  GaugeMatrix (*gauge_field)[4] = (GaugeMatrix(*)[4])a->u;
  FourSpinor *psi = (FourSpinor*)a->spinor;
  HalfSpinor *chi = (HalfSpinor*)a->half_spinor;
  GaugeMatrix *up ALIGN;
  HalfSpinor *s3 ALIGN, *s4 ALIGN;
  FourSpinor part_sum ALIGN, *result ALIGN;

  /* printf("ID: %d low=%d high=%d: MvvReconsPlus\n", id, low, high); */
 


  for (ix1=low;ix1<high;ix1++) {
    int thissite=tab->siteTable(ix1);
    result=&psi[thissite];

    up=&gauge_field[thissite][0];
    s3 = tab->halfspinorBufferOffset(RECONS_MVV_GATHER,ix1,0);
    mvv_recons_gamma0_plus(*s3, *up, part_sum);    

    up=&gauge_field[thissite][1];
    s3 = tab->halfspinorBufferOffset(RECONS_MVV_GATHER,ix1,1);
    mvv_recons_gamma1_plus_add(*s3, *up, part_sum);    

    up=&gauge_field[thissite][2];
    s3 = tab->halfspinorBufferOffset(RECONS_MVV_GATHER,ix1,2);
    mvv_recons_gamma2_plus_add(*s3, *up, part_sum);    

    up=&gauge_field[thissite][3];
    s3 = tab->halfspinorBufferOffset(RECONS_MVV_GATHER,ix1,3);
    mvv_recons_gamma3_plus_add_store(*s3, *up, part_sum, *result);    

  }
}


   
/* this routine takes the partial sum from mvv_recons() and loops 
 * over the output spin components, 2 at a time doing a sum over directions 
 * for each set, accumulating in xmm0-2 and loading the halfspinor 
 * temporaries into xmm3-5 */

 void recons_plus(size_t lo,size_t hi, int id, const void *ptr )	
{
  ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr; 
  ShiftTable<HalfSpinor>* tab = (ShiftTable<HalfSpinor>*)a->s;
  int subgrid_vol_cb = tab->subgridVolCB();

  int cb = a->cb; 
  int ix1=0;
  const int low = cb*subgrid_vol_cb + lo; 
  const int high = cb*subgrid_vol_cb + hi;

  FourSpinor *psi = (FourSpinor*)a->spinor;
  HalfSpinor *chi = (HalfSpinor*)a->half_spinor;

  HalfSpinor *hs0 ALIGN , *hs1 ALIGN, *hs2 ALIGN, *hs3 ALIGN;
  FourSpinor *rn ALIGN;

  /*  printf("ID: %d low=%d high=%d: ReconsPlus\n", id, low, high); */
  for (ix1=low;ix1<high;ix1++) 
  {
    int thissite = tab->siteTable(ix1);
    rn=&psi[thissite];

    /* first spin component of result */
    hs0 = tab->halfspinorBufferOffset(RECONS_GATHER,ix1,0);
    hs1 = tab->halfspinorBufferOffset(RECONS_GATHER,ix1,1);	  
    hs2 = tab->halfspinorBufferOffset(RECONS_GATHER,ix1,2);
    hs3 = tab->halfspinorBufferOffset(RECONS_GATHER,ix1,3);
    recons_4dir_plus(*hs0, *hs1, *hs2, *hs3, *rn);
  }
}


/*************** now for isign corresponding to -1  ****************************************/

 void decomp_minus(size_t lo,size_t hi, int id, const void *ptr ) /*need to fix decomp_minus */
{
  ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr; 
  ShiftTable<HalfSpinor>* tab = (ShiftTable<HalfSpinor>*)a->s;
  int subgrid_vol_cb = tab->subgridVolCB();

  int cb = a->cb; 
  int ix1=0;
  const int low = cb*subgrid_vol_cb + lo; 
  const int high = cb*subgrid_vol_cb + hi;

  FourSpinor *psi = (FourSpinor*)a->spinor;
  HalfSpinor *chi = (HalfSpinor*)a->half_spinor;

  HalfSpinor *s3 ALIGN; 
  FourSpinor *sp ALIGN;
 
  /*  printf("ID: %d low=%d high=%d: DecompMinus\n", id, low, high); */
  for (ix1=low;ix1<high;ix1++) {
    int thissite = tab->siteTable(ix1);

    sp=&psi[thissite];

    s3 = tab->halfspinorBufferOffset(DECOMP_SCATTER,ix1,0);
    decomp_gamma0_plus(*sp, *s3);

    s3 = tab->halfspinorBufferOffset(DECOMP_SCATTER,ix1,1);
    decomp_gamma1_plus(*sp, *s3);
    
    s3 = tab->halfspinorBufferOffset(DECOMP_SCATTER,ix1,2);
    decomp_gamma2_plus(*sp, *s3);
    
    s3 = tab->halfspinorBufferOffset(DECOMP_SCATTER,ix1,3);
    decomp_gamma3_plus(*sp, *s3); 
  
  }
}


/* need gauge fields on opposite cb */
 void decomp_hvv_minus(size_t lo,size_t hi, int id, const void *ptr )
{
  ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr; 
  ShiftTable<HalfSpinor>* tab = (ShiftTable<HalfSpinor>*)a->s;
  int subgrid_vol_cb = tab->subgridVolCB();

  int cb = a->cb; 
  int ix1 = 0;
  const int low = cb*subgrid_vol_cb + lo; 
  const int high = cb*subgrid_vol_cb + hi;

  GaugeMatrix (*gauge_field)[4] = (GaugeMatrix(*)[4])a->u;
  FourSpinor *psi = (FourSpinor*)a->spinor;
  HalfSpinor *chi = (HalfSpinor*)a->half_spinor;

  GaugeMatrix *um ALIGN;
  HalfSpinor *s3 ALIGN;
  FourSpinor  *sm ALIGN;

  /* printf("ID: %d low=%d high=%d: DecompHvvMinus\n", id, low, high); */

  for (ix1=low;ix1<high;ix1++) 
  {
    int thissite = tab->siteTable(ix1);
    sm=&psi[thissite];

    um=&gauge_field[thissite][0];
    s3 = tab->halfspinorBufferOffset(DECOMP_HVV_SCATTER,ix1,0);
    decomp_hvv_gamma0_minus(*sm, *um, *s3);	   

    um=&gauge_field[thissite][1];
    s3 = tab->halfspinorBufferOffset(DECOMP_HVV_SCATTER,ix1,1);
    decomp_hvv_gamma1_minus(*sm, *um, *s3);	   

    um=&gauge_field[thissite][2];
    s3 = tab->halfspinorBufferOffset(DECOMP_HVV_SCATTER,ix1,2);
    decomp_hvv_gamma2_minus(*sm, *um, *s3);	   
    
    um=&gauge_field[thissite][3];
    s3 = tab->halfspinorBufferOffset(DECOMP_HVV_SCATTER,ix1,3);
    decomp_hvv_gamma3_minus(*sm, *um, *s3);	   
  }
}


 void mvv_recons_minus(size_t lo,size_t hi, int id, const void *ptr )
{
  ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr; 

  ShiftTable<HalfSpinor>* tab = (ShiftTable<HalfSpinor>*)a->s;
  int subgrid_vol_cb = tab->subgridVolCB();
 
  int cb = a->cb; 
  int ix1;

  const int low = cb*subgrid_vol_cb + lo; 
  const int high = cb*subgrid_vol_cb + hi;

  GaugeMatrix (*gauge_field)[4] = (GaugeMatrix(*)[4])a->u;
  FourSpinor *psi = (FourSpinor*)a->spinor;
  HalfSpinor *chi = (HalfSpinor*)a->half_spinor;

  GaugeMatrix *up ALIGN;

  HalfSpinor *s3 ALIGN;

  FourSpinor rs ALIGN;
  FourSpinor *rn ALIGN;

  /*  printf("ID: %d low=%d high=%d: MvvReconsMinus\n", id, low, high); */

  for (ix1=low;ix1<high;ix1++) {
    int thissite = tab->siteTable(ix1);
    rn=&psi[thissite];

    up = &gauge_field[thissite][0];
    s3 = tab->halfspinorBufferOffset(RECONS_MVV_GATHER,ix1,0);
    mvv_recons_gamma0_minus(*s3, *up, rs);


    up = &gauge_field[thissite][1];
    s3 = tab->halfspinorBufferOffset(RECONS_MVV_GATHER,ix1,1);
    mvv_recons_gamma1_minus_add(*s3, *up, rs);


    up = &gauge_field[thissite][2];
    s3 = tab->halfspinorBufferOffset(RECONS_MVV_GATHER,ix1,2);
     mvv_recons_gamma2_minus_add(*s3, *up, rs);

    up = &gauge_field[thissite][3];
    s3 = tab->halfspinorBufferOffset(RECONS_MVV_GATHER,ix1,3);
    mvv_recons_gamma3_minus_add_store(*s3, *up, rs,*rn);

  }

}
/******************end of mvv_recons*************************/


 void recons_minus(size_t lo,size_t hi, int id, const void *ptr )	
{
  ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr; 
  ShiftTable<HalfSpinor>* tab = (ShiftTable<HalfSpinor>*)a->s;
  int subgrid_vol_cb = tab->subgridVolCB();

  int cb = a->cb; 
  int ix1=0;
  const int low = cb*subgrid_vol_cb + lo; 
  const int high = cb*subgrid_vol_cb + hi;

  FourSpinor *psi = (FourSpinor*)a->spinor;
  HalfSpinor *chi = (HalfSpinor*)a->half_spinor;

  HalfSpinor *hs0 ALIGN;
  HalfSpinor *hs1 ALIGN;
  HalfSpinor *hs2 ALIGN;
  HalfSpinor *hs3 ALIGN;   
  FourSpinor  *s1 ALIGN,  *rn ALIGN;

  for (ix1=low; ix1<high; ix1++)  {
    int thissite = tab->siteTable(ix1);
    rn=&psi[thissite];

    /* first spin component of result */
    hs0 = tab->halfspinorBufferOffset(RECONS_GATHER,ix1,0);
    hs1 = tab->halfspinorBufferOffset(RECONS_GATHER,ix1,1);
    hs2 = tab->halfspinorBufferOffset(RECONS_GATHER,ix1,2);
    hs3 = tab->halfspinorBufferOffset(RECONS_GATHER,ix1,3);
    recons_4dir_minus(*hs0, *hs1, *hs2, *hs3, *rn);
    /* end of loop */
  }
}
/*****************end of isign corresponding to -1 **********************/
  }
}


