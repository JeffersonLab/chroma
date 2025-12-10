#include "cpp_dslash_parscalar_utils_32bit.h"
#include "shift_table_parscalar.h"
#include "dispatch_parscalar.h"
#include "cpp_dslash_types.h"


namespace CPlusPlusWilsonDslash { 
  using namespace Dslash32BitTypes;

  namespace DslashParscalar32Bit { 
  
    
    void decomp_plus(size_t lo,size_t hi, int id, const void *ptr) /*need to fix decomp_minus */
    {
      int ix1, ix2, iz1;                           /* Site index - iz1 used at loop end */
    FourSpinor* sp1 ALIGN;                /* Spinor under consideration */
    FourSpinor* sp2 ALIGN;
    const ThreadWorkerArgs *a = (ThreadWorkerArgs *)ptr;          /* Cast the argument */
    ShiftTable<HalfSpinor>* s=(ShiftTable<HalfSpinor>*)a->s;
    int subgrid_vol_cb=s->subgridVolCB();

    HalfSpinor* chi =(HalfSpinor *)a->half_spinor; /* needs to be changed to HalfSpinor and be an array*/
    int cb = a->cb;
    int low = cb*subgrid_vol_cb + lo;
    int high = cb*subgrid_vol_cb + hi;


    HalfSpinor* s3 ALIGN;
    HalfSpinor* s4 ALIGN;
    HalfSpinor* s5 ALIGN;
    HalfSpinor* s6 ALIGN;

    
    FourSpinor* spinor_field=(FourSpinor *)a->spinor;

    int thissite;
    

    /************************ loop over all lattice sites *************************/
    thissite = s->siteTable(low);
    sp1=&spinor_field[thissite]; 
    PREFETCH(sp1, _MM_HINT_T0);
    
 
    s3 =  s->halfspinorBufferOffset(DECOMP_SCATTER,low,0);
    PREFETCH(s3, _MM_HINT_T0);

    s4 =  s->halfspinorBufferOffset(DECOMP_SCATTER,low,1);
    PREFETCH(s4, _MM_HINT_T0);

    s5 =  s->halfspinorBufferOffset(DECOMP_SCATTER,low,2);
    PREFETCH(s5, _MM_HINT_T0);


    s6 =  s->halfspinorBufferOffset(DECOMP_SCATTER,low,3);
    PREFETCH(s6, _MM_HINT_T0);

    for (ix1=low+1; ix1<high;ix1++) {
      thissite=s->siteTable(ix1); // Next site  
      sp2=&spinor_field[thissite]; // For prefetching
      PREFETCH(sp2, _MM_HINT_T0);
  
      decomp_gamma0_minus(sp1[0], *s3);
      s3 =  s->halfspinorBufferOffset(DECOMP_SCATTER,ix1,0);
      PREFETCH(s3, _MM_HINT_T0);

      decomp_gamma1_minus(sp1[0], *s4);
      s4 =  s->halfspinorBufferOffset(DECOMP_SCATTER,ix1,1);
      PREFETCH(s4, _MM_HINT_T0);

      decomp_gamma2_minus(sp1[0], *s5);
      s5 =  s->halfspinorBufferOffset(DECOMP_SCATTER,ix1,2);
      PREFETCH(s5, _MM_HINT_T0);




      decomp_gamma3_minus(sp1[0], *s6);
      s6 =  s->halfspinorBufferOffset(DECOMP_SCATTER,ix1,3);
      PREFETCH(s6, _MM_HINT_T0);
      
      sp1=sp2; // For prefetching
    }

    decomp_gamma0_minus(sp1[0], *s3);
    decomp_gamma1_minus(sp1[0], *s4);
    decomp_gamma2_minus(sp1[0], *s5);
    decomp_gamma3_minus(sp1[0], *s6);


  }


/* the basic operations in this routine include loading a spinor, doing 
 * the spin projection, and multiplying the halfspinor by the appropriate 
 * gauge field, and saving the resulting halfspinor to a lattice temporary */

/* need gauge fields on opposite cb */

void decomp_hvv_plus(size_t lo,size_t hi, int id, const void *ptr)
{

  int ix1, iz1;              /* Site addresses. ix1 = current. 
				iz1 is for next loop iteration to allow some loop peeling
			        with ix1 */

  GaugeMatrix* um1 ALIGN;    /* Gauge pointer for 1 site */
  GaugeMatrix* um2 ALIGN;    /* Gauge pointer for 2nd site */
  GaugeMatrix* um3 ALIGN;    /* Temporary gauge pointer for prefetching */
  GaugeMatrix* um4 ALIGN;

  FourSpinor* sm1 ALIGN;   /* spinor */
  FourSpinor* sm2 ALIGN;

  const ThreadWorkerArgs *a = (const ThreadWorkerArgs *)ptr;
  ShiftTable<HalfSpinor>* s=(ShiftTable<HalfSpinor>*)a->s;
  int subgrid_vol_cb=s->subgridVolCB();
  FourSpinor* spinor_field = (FourSpinor*)a->spinor;
  HalfSpinor* chi = (HalfSpinor *)a->half_spinor; /* a 1-d std::map of a 2-d array */
  GaugeMatrix  (*gauge_field)[4] = (GaugeMatrix (*)[4])a->u;

  HalfSpinor* s3 ALIGN;
  HalfSpinor* s4 ALIGN;
  HalfSpinor* s5 ALIGN;
  HalfSpinor* s6 ALIGN;

  int cb = a->cb;

  int low = cb*subgrid_vol_cb + lo;
  int high = cb*subgrid_vol_cb + hi;


  /************************ loop over all lattice sites *************************/
  int thissite = s->siteTable(low);
  s3 =  s->halfspinorBufferOffset(DECOMP_HVV_SCATTER,low,0);
  s4 =  s->halfspinorBufferOffset(DECOMP_HVV_SCATTER,low,1);    
  s5 =  s->halfspinorBufferOffset(DECOMP_HVV_SCATTER,low,2);    
  s6 =  s->halfspinorBufferOffset(DECOMP_HVV_SCATTER,low,3); 
  um1=&gauge_field[thissite][0];
  um2=&gauge_field[thissite][1];
  um3=&gauge_field[thissite][2];
  um4=&gauge_field[thissite][3];

  PREFETCH(s3, _MM_HINT_T0);
  PREFETCH(um1,_MM_HINT_T0);

  PREFETCH(s4, _MM_HINT_T0);
  PREFETCH(um2,_MM_HINT_T0);

  PREFETCH(s5, _MM_HINT_T0);
  PREFETCH(um3, _MM_HINT_T0);

  PREFETCH(s6, _MM_HINT_T0);
  PREFETCH(um4, _MM_HINT_T0);

  sm1=&spinor_field[thissite];

  for (ix1=low+1;ix1<high;ix1++) {
    thissite=s->siteTable(ix1); // Next site
    sm2=&spinor_field[thissite]; 


    /****************** direction +0 *********************************/
    decomp_hvv_gamma0_plus(*sm1,*um1,*s3);
    s3 =  s->halfspinorBufferOffset(DECOMP_HVV_SCATTER,ix1,0);
    um1=&gauge_field[thissite][0];
    PREFETCH(s3, _MM_HINT_T0);
    PREFETCH(um1,_MM_HINT_T0);

    decomp_hvv_gamma1_plus(*sm1,*um2,*s4);
    s4 =  s->halfspinorBufferOffset(DECOMP_HVV_SCATTER,ix1,1);
    um2=&gauge_field[thissite][1];
    PREFETCH(s4, _MM_HINT_T0);
    PREFETCH(um2, _MM_HINT_T0);

    PREFETCH(sm1,_MM_HINT_T0);

    decomp_hvv_gamma2_plus(*sm1,*um3,*s5);
    s5 =  s->halfspinorBufferOffset(DECOMP_HVV_SCATTER,ix1,2);
    um3=&gauge_field[thissite][2];
    PREFETCH(s5, _MM_HINT_T0);
    PREFETCH(um3, _MM_HINT_T0);


    decomp_hvv_gamma3_plus(*sm1,*um4,*s6);
    s6 =  s->halfspinorBufferOffset(DECOMP_HVV_SCATTER,ix1,3);
    um4=&gauge_field[thissite][3];
    PREFETCH(s6, _MM_HINT_T0);
    PREFETCH(um4, _MM_HINT_T0);
    
    sm1=sm2;
  }
  decomp_hvv_gamma0_plus(*sm1,*um1,*s3);
  decomp_hvv_gamma1_plus(*sm1,*um2,*s4);

  PREFETCH(sm1,_MM_HINT_T0);

  decomp_hvv_gamma2_plus(*sm1,*um3,*s5);
  decomp_hvv_gamma3_plus(*sm1,*um4,*s6);

}
/***************end of decomp_hvv****************/


/* the basic operations in this routine include loading the halfspinor 
 * from memory, multiplying it by the appropriate gauge field, doing the 
 * spin reconstruction, and summing over directions, and saving the partial 
 * sum over directions */

void mvv_recons_plus(size_t lo,size_t hi, int id, const void *ptr)
{

  int ix1, iz1;

  GaugeMatrix* up1 ALIGN;
  GaugeMatrix* up2 ALIGN;
  GaugeMatrix* up3 ALIGN;
  GaugeMatrix* up4 ALIGN;

  FourSpinor* sn1 ALIGN;
  FourSpinor* sn2 ALIGN;

  HalfSpinor r12_1 ALIGN, r34_1 ALIGN,r12_2 ALIGN,r34_2 ALIGN;


  const ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr;
  ShiftTable<HalfSpinor>* s=(ShiftTable<HalfSpinor>*)a->s;
  int subgrid_vol_cb=s->subgridVolCB();

  FourSpinor* spinor_field = (FourSpinor *)a->spinor;
  HalfSpinor* chi =(HalfSpinor *)a->half_spinor; /* a 1-d std::map of a 2-d array */
  GaugeMatrix (*gauge_field)[4] = (GaugeMatrix (*)[4])a->u;

  int cb = a->cb;


  HalfSpinor* s3 ALIGN;
  HalfSpinor* s4 ALIGN;
  HalfSpinor* s5 ALIGN;
  HalfSpinor* s6 ALIGN;

  int low = cb*subgrid_vol_cb + lo;
  int high = cb*subgrid_vol_cb + hi;
  
  /************************ loop over all lattice sites *************************/
  int thissite = s->siteTable( low );
  s3 =  s->halfspinorBufferOffset(RECONS_MVV_GATHER,low,0);
  s4 =  s->halfspinorBufferOffset(RECONS_MVV_GATHER,low,1);
  s5 =  s->halfspinorBufferOffset(RECONS_MVV_GATHER,low,2);
  s6 =  s->halfspinorBufferOffset(RECONS_MVV_GATHER,low,3);
  up1=&gauge_field[thissite][0]; 
  up2=&gauge_field[thissite][1];
  up3=&gauge_field[thissite][2];
  up4=&gauge_field[thissite][3];

  PREFETCH(s3, _MM_HINT_T0);
  PREFETCH(up1,_MM_HINT_T0);

  PREFETCH(s4, _MM_HINT_T0);
  PREFETCH(up2,_MM_HINT_T0);

  PREFETCH(s5, _MM_HINT_T0);
  PREFETCH(up3, _MM_HINT_T0);

  PREFETCH(s6, _MM_HINT_T0);
  PREFETCH(up4, _MM_HINT_T0);

  sn1=&spinor_field[thissite];    


  for (ix1=low+1;ix1<high;ix1++) {

    thissite=s->siteTable(ix1);
    sn2=&spinor_field[thissite];    

    mvv_recons_gamma0_plus(*s3, *up1, r12_1, r34_1);
    s3 =  s->halfspinorBufferOffset(RECONS_MVV_GATHER,ix1,0);
    up1=&gauge_field[thissite][0]; 
    PREFETCH(s3, _MM_HINT_T0);
    PREFETCH(up1,_MM_HINT_T0);

    mvv_recons_gamma1_plus_add(*s4, *up2, r12_1, r34_1);
    s4 =  s->halfspinorBufferOffset(RECONS_MVV_GATHER,ix1,1);
    up2=&gauge_field[thissite][1]; 
    PREFETCH(s4, _MM_HINT_T0);
    PREFETCH(up2, _MM_HINT_T0);

    
    PREFETCH(sn1,_MM_HINT_T0);

    mvv_recons_gamma2_plus_add(*s5, *up3, r12_1, r34_1);
    s5 =  s->halfspinorBufferOffset(RECONS_MVV_GATHER,ix1,2);
    up3=&gauge_field[thissite][2];
    PREFETCH(s5, _MM_HINT_T0);
    PREFETCH(up3, _MM_HINT_T0);

    mvv_recons_gamma3_plus_add_store(*s6, *up4, r12_1, r34_1,*sn1);
    s6 =  s->halfspinorBufferOffset(RECONS_MVV_GATHER,ix1,3);
    up4=&gauge_field[thissite][3];
    PREFETCH(s6, _MM_HINT_T0);
    PREFETCH(up4, _MM_HINT_T0);
  
    sn1=sn2;
  }

  mvv_recons_gamma0_plus(*s3, *up1, r12_1, r34_1);
  mvv_recons_gamma1_plus_add(*s4, *up2, r12_1, r34_1);

  PREFETCH(sn1,_MM_HINT_T0);

  mvv_recons_gamma2_plus_add(*s5, *up3, r12_1, r34_1);
  mvv_recons_gamma3_plus_add_store(*s6, *up4, r12_1, r34_1,*sn1);

}


   
/* this routine takes the partial sum from mvv_recons() and loops 
 * over the output spin components, 2 at a time doing a sum over directions 
 * for each set, accumulating in xmm0-2 and loading the halfspinor 
 * temporaries into xmm3-5 */

void recons_plus(size_t lo,size_t hi, int id, const void *ptr )	
{
  int ix1;
  FourSpinor* sn1 ALIGN;
  FourSpinor* sn2 ALIGN;
  

  const ThreadWorkerArgs *a = (ThreadWorkerArgs *)ptr;
  ShiftTable<HalfSpinor>* s=(ShiftTable<HalfSpinor>*)a->s;
  int subgrid_vol_cb=s->subgridVolCB();

  FourSpinor* spinor_field = (FourSpinor *)a->spinor;
  HalfSpinor* chi = (HalfSpinor *)a->half_spinor;
  int cb = a->cb;

  HalfSpinor *hs0 ALIGN;
  HalfSpinor *hs1 ALIGN;
  HalfSpinor *hs2 ALIGN;
  HalfSpinor *hs3 ALIGN;

  int low = cb*subgrid_vol_cb + lo;
  int high = cb*subgrid_vol_cb + hi;
  
  const int PREFDIST=4;

  /************************ loop over all lattice sites *************************/
  int thissite = s->siteTable(low);



  hs0 =  s->halfspinorBufferOffset(RECONS_GATHER,low,0); 
  PREFETCH(hs0, _MM_HINT_NTA);

  hs1 =  s->halfspinorBufferOffset(RECONS_GATHER,low,1); 
  PREFETCH(hs1, _MM_HINT_NTA);

  hs2 =  s->halfspinorBufferOffset(RECONS_GATHER,low,2); 
  PREFETCH(hs2, _MM_HINT_NTA);

  hs3 =  s->halfspinorBufferOffset(RECONS_GATHER,low,3);
  PREFETCH(hs3, _MM_HINT_NTA);

  sn1=&spinor_field[thissite];   
  PREFETCH(sn1, _MM_HINT_NTA);

  for (ix1=low+1;ix1<high;ix1++) {
    thissite = s->siteTable(ix1);
    sn2 = &spinor_field[thissite];   
    PREFETCH(sn2, _MM_HINT_NTA);
   
    recons_4dir_plus(*hs0, *hs1, *hs2, *hs3, *sn1);

    hs0 =  s->halfspinorBufferOffset(RECONS_GATHER,ix1,0); 
    PREFETCH(hs0, _MM_HINT_NTA);

    hs1 =  s->halfspinorBufferOffset(RECONS_GATHER,ix1,1); 
    PREFETCH(hs1, _MM_HINT_NTA);

    hs2 =  s->halfspinorBufferOffset(RECONS_GATHER,ix1,2); 
    PREFETCH(hs2, _MM_HINT_NTA);

    hs3 =  s->halfspinorBufferOffset(RECONS_GATHER,ix1,3); 
    PREFETCH(hs3, _MM_HINT_NTA);

    sn1=sn2;
    /*************************end of loop ****************************/
  }
  recons_4dir_plus(*hs0, *hs1, *hs2, *hs3, *sn1);


}


/*************** now for isign corresponding to -1  ****************************************/

void decomp_minus(size_t lo,size_t hi, int id, const void *ptr ) /*need to fix decomp_minus */
{

  int ix1,iz1;

  FourSpinor* sp1 ALIGN;
  FourSpinor* sp2 ALIGN;

  const ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr;
  ShiftTable<HalfSpinor>* s=(ShiftTable<HalfSpinor>*)a->s;
  int subgrid_vol_cb=s->subgridVolCB();

  HalfSpinor* chi = (HalfSpinor *)a->half_spinor; /* needs to be changed to HalfSpinor and be an array*/

  HalfSpinor* s3 ALIGN;
  HalfSpinor* s4 ALIGN;
  HalfSpinor* s5 ALIGN;
  HalfSpinor* s6 ALIGN;

  int cb = a->cb;

  int low = cb*subgrid_vol_cb + lo;
  int high = cb*subgrid_vol_cb + hi;

  FourSpinor* spinor_field=(FourSpinor *) a->spinor;

  int thissite;

  /************************ loop over all lattice sites *************************/

  thissite = s->siteTable(low);
  sp1=&spinor_field[thissite]; 
  PREFETCH(sp1, _MM_HINT_T0);
  
  s3 =  s->halfspinorBufferOffset(DECOMP_SCATTER,low,0);
  PREFETCH(s3, _MM_HINT_T0);
  
  s4 =  s->halfspinorBufferOffset(DECOMP_SCATTER,low,1);
  PREFETCH(s4, _MM_HINT_T0);
  
  s5 =  s->halfspinorBufferOffset(DECOMP_SCATTER,low,2);
  PREFETCH(s5, _MM_HINT_T0);
  
  s6 =  s->halfspinorBufferOffset(DECOMP_SCATTER,low,3);
  PREFETCH(s6, _MM_HINT_T0);

  
  for (ix1=low+1;ix1<high;ix1++) {
    thissite=s->siteTable(ix1); // Next site  
    sp2=&spinor_field[thissite]; // For prefetching
    PREFETCH(sp2, _MM_HINT_T0);
  
    
    decomp_gamma0_plus(sp1[0], *s3);
    s3 =  s->halfspinorBufferOffset(DECOMP_SCATTER,ix1,0);
    PREFETCH(s3, _MM_HINT_T0);

    decomp_gamma1_plus(sp1[0], *s4);
    s4 =  s->halfspinorBufferOffset(DECOMP_SCATTER,ix1,1);
    PREFETCH(s4, _MM_HINT_T0);
      
    decomp_gamma2_plus(sp1[0], *s5);
    s5 =  s->halfspinorBufferOffset(DECOMP_SCATTER,ix1,2);
    PREFETCH(s5, _MM_HINT_T0);

    decomp_gamma3_plus(sp1[0], *s6);    
    s6 =  s->halfspinorBufferOffset(DECOMP_SCATTER,ix1,3);
    PREFETCH(s6, _MM_HINT_T0);
      
    sp1=sp2; // For prefetching
  }
  decomp_gamma0_plus(sp1[0], *s3);
  decomp_gamma1_plus(sp1[0], *s4);
  decomp_gamma2_plus(sp1[0], *s5);
  decomp_gamma3_plus(sp1[0], *s6);
  
}


/* need gauge fields on opposite cb */

void decomp_hvv_minus(size_t lo,size_t hi, int id, const void *ptr )
{

  int ix1,iz1;
  GaugeMatrix* um1 ALIGN;
  GaugeMatrix* um2 ALIGN;
  GaugeMatrix* um3 ALIGN;
  GaugeMatrix* um4 ALIGN;

  FourSpinor* sm1 ALIGN;
  FourSpinor* sm2 ALIGN;



  const ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr;
  ShiftTable<HalfSpinor>* s=(ShiftTable<HalfSpinor>*)a->s;
  int subgrid_vol_cb=s->subgridVolCB();

  FourSpinor* spinor_field = (FourSpinor *)a->spinor;
  HalfSpinor* chi = (HalfSpinor *)a->half_spinor; /* a 1-d std::map of a 2-d array */
  GaugeMatrix (*gauge_field)[4] = (GaugeMatrix(*)[4])a->u;

  int cb = a->cb;

  HalfSpinor* s3 ALIGN;
  HalfSpinor* s4 ALIGN;
  HalfSpinor* s5 ALIGN;
  HalfSpinor* s6 ALIGN;

  int low = cb*subgrid_vol_cb + lo;
  int high = cb*subgrid_vol_cb + hi;

  /************************ loop over all lattice sites *************************/
  int thissite = s->siteTable(low);
  s3 =  s->halfspinorBufferOffset(DECOMP_HVV_SCATTER,low,0);
  s4 =  s->halfspinorBufferOffset(DECOMP_HVV_SCATTER,low,1);    
  s5 =  s->halfspinorBufferOffset(DECOMP_HVV_SCATTER,low,2);    
  s6 =  s->halfspinorBufferOffset(DECOMP_HVV_SCATTER,low,3); 
  um1=&gauge_field[thissite][0];
  um2=&gauge_field[thissite][1];
  um3=&gauge_field[thissite][2];
  um4=&gauge_field[thissite][3];

  PREFETCH(s3, _MM_HINT_T0);
  PREFETCH(um1,_MM_HINT_T0);

  PREFETCH(s4, _MM_HINT_T0);
  PREFETCH(um2,_MM_HINT_T0);

  PREFETCH(s5, _MM_HINT_T0);
  PREFETCH(um3, _MM_HINT_T0);

  PREFETCH(s6, _MM_HINT_T0);
  PREFETCH(um4, _MM_HINT_T0);


  sm1 = &spinor_field[thissite];

  for (ix1=low+1;ix1<high;ix1++) {
    thissite = s->siteTable(ix1);
    sm2=&spinor_field[thissite]; 

    /***************** direction +0 *********************************/
    decomp_hvv_gamma0_minus(*sm1, *um1, *s3);
    s3 =  s->halfspinorBufferOffset(DECOMP_HVV_SCATTER,ix1,0);
    um1=&gauge_field[thissite][0];
    PREFETCH(s3, _MM_HINT_T0);
    PREFETCH(um1,_MM_HINT_T0);

    decomp_hvv_gamma1_minus(*sm1, *um2, *s4);
    s4 =  s->halfspinorBufferOffset(DECOMP_HVV_SCATTER,ix1,1);
    um2=&gauge_field[thissite][1];
    PREFETCH(s4, _MM_HINT_T0);
    PREFETCH(um2, _MM_HINT_T0);


    decomp_hvv_gamma2_minus(*sm1, *um3, *s5);
    s5 =  s->halfspinorBufferOffset(DECOMP_HVV_SCATTER,ix1,2);
    um3=&gauge_field[thissite][2];
    PREFETCH(s5, _MM_HINT_T0);
    PREFETCH(um3, _MM_HINT_T0);

    PREFETCH(sm1,_MM_HINT_T0);
    
    decomp_hvv_gamma3_minus(*sm1, *um4, *s6);
    s6 =  s->halfspinorBufferOffset(DECOMP_HVV_SCATTER,ix1,3);
    um4=&gauge_field[thissite][3];

    PREFETCH(s6, _MM_HINT_T0);
    PREFETCH(um4, _MM_HINT_T0);
  
    sm1=sm2;

  }
  decomp_hvv_gamma0_minus(*sm1, *um1, *s3);
  decomp_hvv_gamma1_minus(*sm1, *um2, *s4);

  PREFETCH(sm1,_MM_HINT_T0);

  decomp_hvv_gamma2_minus(*sm1, *um3, *s5);
  decomp_hvv_gamma3_minus(*sm1, *um4, *s6);

}


void mvv_recons_minus(size_t lo,size_t hi, int id, const void *ptr )
{
  int ix1, iz1;
  FourSpinor* sn1 ALIGN;  /* The spinor to store to */
  FourSpinor* sn2 ALIGN; 
  GaugeMatrix* um1 ALIGN;
  GaugeMatrix* um2 ALIGN;
  GaugeMatrix* um3 ALIGN;
  GaugeMatrix* um4 ALIGN;

  /* Temporaries for the top and bottom parts of spinors. */
  HalfSpinor r12_1 ALIGN, r34_1 ALIGN, r12_2 ALIGN,r34_2 ALIGN;

  /* if going to support unpacked gauge fields, need to treat site ix1 and site ix1+1 separately */
  /* to support unpacked gauge fields the prefetches will need to be changed */
  const ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr;
 
  ShiftTable<HalfSpinor>* s=(ShiftTable<HalfSpinor>*)a->s;
  int subgrid_vol_cb=s->subgridVolCB();

  FourSpinor* spinor_field = (FourSpinor *)a->spinor;
  HalfSpinor* chi = (HalfSpinor *)a->half_spinor; /* a 1-d std::map of a 2-d array */
  GaugeMatrix (*gauge_field)[4] =(GaugeMatrix(*)[4]) a->u;
  int cb = a->cb;

  HalfSpinor* s3 ALIGN;
  HalfSpinor* s4 ALIGN;
  HalfSpinor* s5 ALIGN;
  HalfSpinor* s6 ALIGN;
  

  int low = cb*subgrid_vol_cb + lo;
  int high = cb*subgrid_vol_cb + hi;

  int thissite = s->siteTable(low);
  s3 =  s->halfspinorBufferOffset(RECONS_MVV_GATHER,low,0);
  s4 =  s->halfspinorBufferOffset(RECONS_MVV_GATHER,low,1);
  s5 =  s->halfspinorBufferOffset(RECONS_MVV_GATHER,low,2);
  s6 =  s->halfspinorBufferOffset(RECONS_MVV_GATHER,low,3);
  um1=&gauge_field[thissite][0]; 
  um2=&gauge_field[thissite][1];
  um3=&gauge_field[thissite][2];
  um4=&gauge_field[thissite][3];

  PREFETCH(s3, _MM_HINT_T0);
  PREFETCH(um1,_MM_HINT_T0);

  PREFETCH(s4, _MM_HINT_T0);
  PREFETCH(um2,_MM_HINT_T0);

  PREFETCH(s5, _MM_HINT_T0);
  PREFETCH(um3, _MM_HINT_T0);

  PREFETCH(s6, _MM_HINT_T0);
  PREFETCH(um4, _MM_HINT_T0);

  sn1=&spinor_field[thissite];    
/************************ loop over all lattice sites *************************/
  for (ix1=low+1;ix1<high;ix1++) {
    thissite = s->siteTable(ix1);
    sn2 = &spinor_field[thissite];   
 
    mvv_recons_gamma0_minus(*s3, *um1, r12_1, r34_1);
    s3 =  s->halfspinorBufferOffset(RECONS_MVV_GATHER,ix1,0);
    um1=&gauge_field[thissite][0];  
    PREFETCH(s3, _MM_HINT_T0);
    PREFETCH(um1,_MM_HINT_T0);

    mvv_recons_gamma1_minus_add(*s4, *um2, r12_1, r34_1);
    s4 =  s->halfspinorBufferOffset(RECONS_MVV_GATHER,ix1,1); 
    um2= &gauge_field[thissite][1];
    PREFETCH(s4, _MM_HINT_T0);
    PREFETCH(um2, _MM_HINT_T0);


    mvv_recons_gamma2_minus_add(*s5, *um3, r12_1, r34_1);
    s5 =  s->halfspinorBufferOffset(RECONS_MVV_GATHER,ix1,2);
    um3=&gauge_field[thissite][2];  
    PREFETCH(s5, _MM_HINT_T0);
    PREFETCH(um3, _MM_HINT_T0);

    PREFETCH(sn1,_MM_HINT_T0);

    mvv_recons_gamma3_minus_add_store(*s6, *um4, r12_1, r34_1,*sn1);

    s6 =  s->halfspinorBufferOffset(RECONS_MVV_GATHER,ix1,3);
    um4=&gauge_field[thissite][3];  
    PREFETCH(s6, _MM_HINT_T0);
    PREFETCH(um4, _MM_HINT_T0);
  
    sn1=sn2;
    /******************************** end of loop *********************************/
  }
  mvv_recons_gamma0_minus(*s3, *um1, r12_1, r34_1);
  mvv_recons_gamma1_minus_add(*s4, *um2, r12_1, r34_1);
  mvv_recons_gamma2_minus_add(*s5, *um3, r12_1, r34_1);
  PREFETCH(sn1,_MM_HINT_T0);
  mvv_recons_gamma3_minus_add_store(*s6, *um4, r12_1, r34_1,*sn1);
}
/******************end of mvv_recons*************************/


void recons_minus(size_t lo,size_t hi, int id, const void *ptr )	
{
  int ix1;
  FourSpinor* sn1 ALIGN;
  FourSpinor* sn2 ALIGN;


  const ThreadWorkerArgs *a = (ThreadWorkerArgs *)ptr;
  ShiftTable<HalfSpinor>* s=(ShiftTable<HalfSpinor>*)a->s;
  int subgrid_vol_cb= s->subgridVolCB();

  FourSpinor* spinor_field = (FourSpinor*)a->spinor;
  HalfSpinor* chi = (HalfSpinor*)a->half_spinor; /* a 1-d std::map of a 2-d array */
  int cb = a->cb;

  HalfSpinor *hs0 ALIGN;
  HalfSpinor *hs1 ALIGN;
  HalfSpinor *hs2 ALIGN;
  HalfSpinor *hs3 ALIGN;

  int low = cb*subgrid_vol_cb + lo;
  int high = cb*subgrid_vol_cb + hi;
  int thissite = s->siteTable(low);  
  hs0 =  s->halfspinorBufferOffset(RECONS_GATHER,low,0); 
  PREFETCH(hs0, _MM_HINT_NTA);

  hs1 =  s->halfspinorBufferOffset(RECONS_GATHER,low,1); 
  PREFETCH(hs1, _MM_HINT_NTA);

  hs2 =  s->halfspinorBufferOffset(RECONS_GATHER,low,2); 
  PREFETCH(hs2, _MM_HINT_NTA);

  hs3 =  s->halfspinorBufferOffset(RECONS_GATHER,low,3);
  PREFETCH(hs3, _MM_HINT_NTA);

  sn1=&spinor_field[thissite];   
  PREFETCH(sn1, _MM_HINT_NTA);

  for (ix1=low+1;ix1<high;ix1++) {
    thissite = s->siteTable(ix1);
     sn2 = &spinor_field[thissite];   
    PREFETCH(sn2, _MM_HINT_NTA);

    recons_4dir_minus(*hs0, *hs1, *hs2, *hs3, *sn1);

    hs0 =  s->halfspinorBufferOffset(RECONS_GATHER,ix1,0); 
    PREFETCH(hs0, _MM_HINT_NTA);
    
    hs1 =  s->halfspinorBufferOffset(RECONS_GATHER,ix1,1); 
    PREFETCH(hs1, _MM_HINT_NTA);
    
    hs2 =  s->halfspinorBufferOffset(RECONS_GATHER,ix1,2); 
    PREFETCH(hs2, _MM_HINT_NTA);
    
    hs3 =  s->halfspinorBufferOffset(RECONS_GATHER,ix1,3); 
    PREFETCH(hs3, _MM_HINT_NTA);
   
    sn1=sn2;

  }
  recons_4dir_minus(*hs0, *hs1, *hs2, *hs3, *sn1);
}
/*****************end of isign corresponding to -1 **********************/



  }
}

