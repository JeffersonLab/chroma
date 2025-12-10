#include <cpp_dslash_parscalar.h>

#include <cstdlib>

#include <cache.h>


#include <cpp_dslash_parscalar_utils_32bit.h>

#include <tables_parscalar.h>
#include <dispatch_parscalar.h>

#include <cpp_dslash_types.h>
using namespace CPlusPlusWilsonDslash::Dslash32BitTypes;

#include <shift_table_3d_parscalar.h>


#define QMP_COMMS 

namespace CPlusPlusWilsonDslash {

ShiftTable3D<Dslash<float>::HalfSpinor>* Dslash3D<float>::s_tab = nullptr;
DslashTables<Dslash<float>::HalfSpinor, 3>* Dslash3D<float>::tab = nullptr;

  namespace DslashParscalar32Bit {

    void decomp_plus_3d(size_t lo,size_t hi, int id, const void *ptr) /*need to fix decomp_minus */
    {
      int ix1;                           /* Site index - iz1 used at loop end */
      FourSpinor* sp1 ALIGN;                /* Spinor under consideration */
      
      const ThreadWorkerArgs *a = (ThreadWorkerArgs *)ptr;          /* Cast the argument */
      ShiftTable3D<HalfSpinor>* s = (ShiftTable3D<HalfSpinor>*)a->s;
      int subgrid_vol_cb = s->subgridVolCB();

      HalfSpinor* chi =(HalfSpinor *)a->half_spinor; /* needs to be changed to HalfSpinor and be an array*/
      FourSpinor* spinor_field = (FourSpinor *)a->spinor;
      int cb = a->cb;
      
      
      HalfSpinor* s3;
      
      int low;
      int high;
      
      low = cb*subgrid_vol_cb + lo;
      high = cb*subgrid_vol_cb + hi;
      
      
      /************************ loop over all lattice sites *************************/
      for (ix1 = low; ix1 < high ; ix1++) {
	
	int thissite = s->siteTable(ix1);
	
	sp1=&spinor_field[ thissite ];
	
	/******************************* direction +0 *********************************/
	/* first of two sites */
	s3 = s->halfspinorBufferOffset(DECOMP_SCATTER, ix1, 0);
	
	decomp_gamma0_minus(sp1[0], *s3);
	
	/******************************* direction +1 *********************************/
	s3 = s->halfspinorBufferOffset(DECOMP_SCATTER, ix1, 1);
	
	decomp_gamma1_minus(sp1[0], *s3);
	
	/******************************* direction +2 *********************************/
	s3 = s->halfspinorBufferOffset(DECOMP_SCATTER, ix1, 2);
	
	decomp_gamma2_minus(sp1[0], *s3);
	
	//      s3 = s->halfspinorBufferOffset(DECOMP_SCATTER, ix1, 3);
	// decomp_gamma3_minus(sp1[0], *s3);
	
      }
    }
    
    
    /* the basic operations in this routine include loading a spinor, doing 
     * the spin projection, and multiplying the halfspinor by the appropriate 
     * gauge field, and saving the resulting halfspinor to a lattice temporary */
    
    /* need gauge fields on opposite cb */
    void decomp_hvv_plus_3d(size_t lo,size_t hi, int id, const void *ptr)
    {
      
      int ix1;                   /* Site addresses. ix1 = current. 
				    iz1 is for next loop iteration to allow some loop peeling
				    with ix1 */
      
      GaugeMatrix* um1 ALIGN;    /* Gauge pointer for 1 site */
      
      FourSpinor* sm1 ALIGN;   /* spinor */
      
      
      ThreadWorkerArgs *a = (ThreadWorkerArgs *)ptr;
      ShiftTable3D<HalfSpinor>* s = (ShiftTable3D<HalfSpinor>*)a->s;
      int subgrid_vol_cb = s->subgridVolCB();

      FourSpinor* spinor_field =(FourSpinor*)a->spinor;
      HalfSpinor* chi = (HalfSpinor*)a->half_spinor; /* a 1-d std::map of a 2-d array */
      GaugeMatrix (*gauge_field)[4] =(GaugeMatrix (*)[4]) a->u;
      
      int cb = a->cb;
      
      HalfSpinor* s3;
      
      int low;
      int high;
      
      
      low = cb*subgrid_vol_cb + lo;
      high = cb*subgrid_vol_cb + hi;
      
      
      
      /************************ loop over all lattice sites *************************/
      for (ix1 =low; ix1 <high ;ix1++) {
	
	int thissite = s->siteTable( ix1 );
	
	sm1=&spinor_field[ thissite ];
	um1 = &gauge_field[ thissite ][0];
	
	s3 = s->halfspinorBufferOffset(DECOMP_HVV_SCATTER, ix1, 0);    
	decomp_hvv_gamma0_plus(*sm1,*um1,*s3);
	
	
	/******************************* direction +1 *********************************/
	um1 = &gauge_field[ thissite ][1];
	s3 = s->halfspinorBufferOffset(DECOMP_HVV_SCATTER, ix1, 1);
	decomp_hvv_gamma1_plus(*sm1,*um1,*s3);
	
	/******************************* direction +2 *********************************/
	um1 = &gauge_field[ thissite ][2];
	s3 = s->halfspinorBufferOffset(DECOMP_HVV_SCATTER, ix1, 2);
	decomp_hvv_gamma2_plus(*sm1,*um1,*s3);
	
	/******************************* direction +3 *********************************/
	// um1 = &gauge_field[ thissite ][3];
	//      s3 = s->halfspinorBufferOffset(DECOMP_HVV_SCATTER, ix1, 3);
	// decomp_hvv_gamma3_plus(*sm1,*um1,*s3);
	
      /******************************** end of loop *********************************/
      }
    }
    /***************end of decomp_hvv****************/
    
    
    /* the basic operations in this routine include loading the halfspinor 
     * from memory, multiplying it by the appropriate gauge field, doing the 
     * spin reconstruction, and summing over directions, and saving the partial 
     * sum over directions */
    
    void mvv_recons_plus_3d(size_t lo,size_t hi, int id, const void *ptr)
    {
      
      int ix1;
      int low;
      int high;
      
      
      GaugeMatrix* up1 ALIGN;
      FourSpinor* sn1 ALIGN;
      HalfSpinor r12_1 ALIGN, r34_1 ALIGN;
      
      const ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr;
      ShiftTable3D<HalfSpinor>* s = (ShiftTable3D<HalfSpinor>*)a->s;
      int subgrid_vol_cb = s->subgridVolCB();

      FourSpinor* spinor_field = (FourSpinor *)a->spinor;
      HalfSpinor* chi = (HalfSpinor*)a->half_spinor; /* a 1-d std::map of a 2-d array */
      
      GaugeMatrix (*gauge_field)[4] = (GaugeMatrix (*)[4])a->u;
      int cb = a->cb;
      
      HalfSpinor* s3;
      
      
      low = cb*subgrid_vol_cb + lo;
      high = cb*subgrid_vol_cb + hi;
      
      /************************ loop over all lattice sites *************************/
      for (ix1 =low; ix1 <high ;ix1++) {
	
	//    int thissite = lookup_site(cb,ix1);
	int thissite = s->siteTable( ix1 );
	
	up1=&gauge_field[thissite][0];
	s3 = s->halfspinorBufferOffset(RECONS_MVV_GATHER, ix1, 0);
	
	mvv_recons_gamma0_plus(*s3, *up1, r12_1, r34_1);
	
	
	/***************************** direction +1 ***********************************/
	up1=&gauge_field[thissite][1];
	s3 = s->halfspinorBufferOffset(RECONS_MVV_GATHER, ix1, 1);
	mvv_recons_gamma1_plus_add(*s3, *up1, r12_1, r34_1);
	
	
	/******************************* direction +2 *********************************/
	up1=&gauge_field[thissite][2];
	s3 = s->halfspinorBufferOffset(RECONS_MVV_GATHER, ix1, 2); 
	sn1=&spinor_field[thissite];
	mvv_recons_gamma2_plus_add_store(*s3, *up1, r12_1, r34_1, *sn1);
	
	/******************************* direction +2 *********************************/
	//    up1=&gauge_field[thissite][3];
	// s3 = s->halfspinorBufferOffset(RECONS_MVV_GATHER, ix1, 3); 
	// sn1=&spinor_field[thissite];
	// mvv_recons_gamma3_plus_add_store(*s3, *up1, r12_1, r34_1, *sn1);
	
	
      }
    }


   
    /* this routine takes the partial sum from mvv_recons() and loops 
     * over the output spin components, 2 at a time doing a sum over directions 
     * for each set, accumulating in xmm0-2 and loading the halfspinor 
     * temporaries into xmm3-5 */
    
    void recons_plus_3d(size_t lo,size_t hi, int id, const void *ptr )	
    {
      int ix1;
      FourSpinor* sn1 ALIGN;
      
      
      const ThreadWorkerArgs *a = (ThreadWorkerArgs *)ptr;
      ShiftTable3D<HalfSpinor>* s = (ShiftTable3D<HalfSpinor>*)a->s;
      int subgrid_vol_cb = s->subgridVolCB();


      FourSpinor* spinor_field = (FourSpinor *)a->spinor;
      HalfSpinor* chi = (HalfSpinor*)a->half_spinor;
      int cb = a->cb;
      
      
      
      HalfSpinor *hs0, *hs1, *hs2;  
      int low;
      int high;
      
      low = cb*subgrid_vol_cb + lo;
      high = cb*subgrid_vol_cb + hi;
      
      
      /************************ loop over all lattice sites *************************/
      for (ix1 =low; ix1 <high ;ix1++) {
	
	//    int thissite=lookup_site(cb,ix1);
	int thissite = s->siteTable( ix1 );
	
	hs0 = s->halfspinorBufferOffset(RECONS_GATHER,ix1,0); 
	hs1 = s->halfspinorBufferOffset(RECONS_GATHER,ix1,1);
	hs2 = s->halfspinorBufferOffset(RECONS_GATHER,ix1,2);
	//    hs3 = s->halfspinorBufferOffset(RECONS_GATHER,ix1,3);
	sn1=&spinor_field[thissite];   
	recons_3dir_plus(*hs0, *hs1, *hs2, *sn1);
	
	/*************************end of loop ****************************/
      }
    }
    /*****************end of recons**************/
    
    
    
    
    
    /*************** now for isign corresponding to -1  ****************************************/
    
    void decomp_minus_3d(size_t lo,size_t hi, int id, const void *ptr ) /*need to fix decomp_minus */
    {
      
      int ix1;
      FourSpinor* sp1 ALIGN;
      
      const ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr;
      ShiftTable3D<HalfSpinor>* s = (ShiftTable3D<HalfSpinor>*)a->s;
      int subgrid_vol_cb = s->subgridVolCB();

      HalfSpinor* chi = (HalfSpinor*)a->half_spinor; /* needs to be changed to HalfSpinor and be an array*/
      
      HalfSpinor* s3;
      
      int cb = a->cb;
      FourSpinor* spinor_field= (FourSpinor *)a->spinor;
      
      /************************ loop over all lattice sites *************************/
      int low = cb*subgrid_vol_cb + lo;
      int high = cb*subgrid_vol_cb + hi;
      
      
      /************************ loop over all lattice sites *************************/
      for (ix1 =low; ix1 <high ;ix1++) {
	
	int thissite = s->siteTable( ix1 );       
	sp1=&spinor_field[thissite];
	
	/******************************* direction +0 *********************************/
	/* ...(1-gamma(0))... */
	s3 = s->halfspinorBufferOffset(DECOMP_SCATTER,ix1,0);
	decomp_gamma0_plus(sp1[0], *s3);
	
	/******************************* direction +1 *********************************/
	s3 = s->halfspinorBufferOffset(DECOMP_SCATTER,ix1,1);    
	decomp_gamma1_plus(sp1[0], *s3);
	
	/******************************* direction +2 *********************************/
	s3 = s->halfspinorBufferOffset(DECOMP_SCATTER,ix1,2);
	decomp_gamma2_plus(sp1[0], *s3);
	
	/******************************* direction +2 *********************************/
	//    s3 = s->halfspinorBufferOffset(DECOMP_SCATTER,ix1,3);
	//    decomp_gamma3_plus(sp1[0], *s3);
	
	
      }
    }
    
    
    /* need gauge fields on opposite cb */
    void decomp_hvv_minus_3d(size_t lo,size_t hi, int id, const void *ptr )
    {
      
      int ix1;
      GaugeMatrix* um1 ALIGN;
      FourSpinor* sm1 ALIGN;
      
      
      const ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr;
      ShiftTable3D<HalfSpinor>* s = (ShiftTable3D<HalfSpinor>*)a->s;
      int subgrid_vol_cb = s->subgridVolCB();
      
      FourSpinor* spinor_field = (FourSpinor *)a->spinor;
      
      HalfSpinor* chi = (HalfSpinor*)a->half_spinor; /* a 1-d std::map of a 2-d array */
      
      GaugeMatrix (*gauge_field)[4] = (GaugeMatrix (*)[4])a->u;
      
      int cb = a->cb;
      
      HalfSpinor* s3;
      
      
      /************************ loop over all lattice sites *************************/
      int low = cb*subgrid_vol_cb + lo;
      int high = cb*subgrid_vol_cb + hi;
      
      
      /************************ loop over all lattice sites *************************/
      for (ix1 =low; ix1 <high ;ix1++) {
	
	int thissite = s->siteTable( ix1 );       
	
	um1=&gauge_field[thissite][0]; 
	sm1=&spinor_field[thissite];
	
	s3 = s->halfspinorBufferOffset(DECOMP_HVV_SCATTER,ix1,0);
	decomp_hvv_gamma0_minus(*sm1, *um1, *s3);
	
	
	/******************************* direction -1 *********************************/
	um1=&gauge_field[thissite][1]; 
	s3 = s->halfspinorBufferOffset(DECOMP_HVV_SCATTER,ix1,1);
	decomp_hvv_gamma1_minus(*sm1, *um1, *s3);
	
	
	/******************************* direction -2 *********************************/    
	um1=&gauge_field[thissite][2]; 
	s3 = s->halfspinorBufferOffset(DECOMP_HVV_SCATTER,ix1,2);
	decomp_hvv_gamma2_minus(*sm1, *um1, *s3);
	
	/******************************* direction -2 *********************************/    
	// um1=&gauge_field[thissite][3]; 
	// s3 = s->halfspinorBufferOffset(DECOMP_HVV_SCATTER,ix1,3);
	// decomp_hvv_gamma3_minus(*sm1, *um1, *s3);
	
      }
    }
    
    
    void mvv_recons_minus_3d(size_t lo,size_t hi, int id, const void *ptr )
    {
      int ix1;
      GaugeMatrix* up1 ALIGN;   /* Gauge pointers for site x and x+1 */
      FourSpinor* sn1 ALIGN;  /* The spinor to store to */
      
      
      /* Temporaries for the top and bottom parts of spinors. */
      HalfSpinor r12_1 ALIGN,r34_1 ALIGN;
      
      /* if going to support unpacked gauge fields, need to treat site ix1 and site ix1+1 separately */
      /* to support unpacked gauge fields the prefetches will need to be changed */
      const ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr;
      ShiftTable3D<HalfSpinor>* s = (ShiftTable3D<HalfSpinor>*)a->s;
      int subgrid_vol_cb = s->subgridVolCB();

      FourSpinor* spinor_field = (FourSpinor *)a->spinor;
      HalfSpinor* chi = (HalfSpinor*)a->half_spinor; /* a 1-d std::map of a 2-d array */
      GaugeMatrix (*gauge_field)[4] = (GaugeMatrix (*)[4])a->u;
      int cb = a->cb;
      
      HalfSpinor* s3;
      
      
      /************************ loop over all lattice sites *************************/
      int low = cb*subgrid_vol_cb + lo;
      int high = cb*subgrid_vol_cb + hi;
      
      
      /************************ loop over all lattice sites *************************/
      for (ix1 =low; ix1 <high ;ix1++) {
	
	int thissite = s->siteTable( ix1 );       
	
	up1=&gauge_field[thissite][0];
	/******************************* direction +0 *********************************/
	s3 = s->halfspinorBufferOffset(RECONS_MVV_GATHER,ix1,0);
	mvv_recons_gamma0_minus(*s3, *up1, r12_1, r34_1);
	
	
	/******************************* direction +1 *********************************/
	up1=&gauge_field[thissite][1];
	s3 = s->halfspinorBufferOffset(RECONS_MVV_GATHER,ix1,1);
	mvv_recons_gamma1_minus_add(*s3, *up1, r12_1, r34_1);
	
	/******************************* direction +2 *********************************/
	up1=&gauge_field[thissite][2];
	s3 = s->halfspinorBufferOffset(RECONS_MVV_GATHER,ix1,2);
	sn1=&spinor_field[thissite];
	mvv_recons_gamma2_minus_add_store(*s3, *up1, r12_1, r34_1, *sn1);
	
	//    up1=&gauge_field[thissite][3];
	// s3 = s->halfspinorBufferOffset(RECONS_MVV_GATHER,ix1,3);
	// mvv_recons_gamma3_minus_add_store(*s3, *up1, r12_1, r34_1, *sn1);
	
	/******************************** end of loop *********************************/
      }
    }
    /******************end of mvv_recons*************************/
    
    
    void recons_minus_3d(size_t lo,size_t hi, int id, const void *ptr )	
    {
      int ix1;
      FourSpinor* sn1 ALIGN;
      
      
      const ThreadWorkerArgs *a = (ThreadWorkerArgs *)ptr;
      ShiftTable3D<HalfSpinor>* s = (ShiftTable3D<HalfSpinor>*)a->s;
      int subgrid_vol_cb = s->subgridVolCB();

      FourSpinor* spinor_field = (FourSpinor *)a->spinor;
      HalfSpinor* chi = (HalfSpinor*)a->half_spinor; /* a 1-d std::map of a 2-d array */
      int cb = a->cb;
      HalfSpinor* hs0, *hs1, *hs2;
      
      /************************ loop over all lattice sites *************************/
      int low = cb*subgrid_vol_cb + lo;
      int high = cb*subgrid_vol_cb + hi;
      
      
      /************************ loop over all lattice sites *************************/
      for (ix1 =low; ix1 <high ;ix1++) {
	
	int thissite = s->siteTable( ix1 );       
	
	sn1=&spinor_field[thissite];   
	
	hs0 = s->halfspinorBufferOffset(RECONS_GATHER,ix1,0); 
	hs1 = s->halfspinorBufferOffset(RECONS_GATHER,ix1,1);
	hs2 = s->halfspinorBufferOffset(RECONS_GATHER,ix1,2);

	
	recons_3dir_minus(*hs0, *hs1, *hs2, *sn1);
      }
      
    }
  }
  



  // Your actual operator
  void Dslash3D<float>::operator()(float* res, 
				   float* psi, 
				   float* u, 
				   int isign,
				   int cb)
  {
    HalfSpinor* chi1 = tab->getChi1();
    HalfSpinor* chi2 = tab->getChi2();
    int subgrid_vol_cb = s_tab->subgridVolCB();
    
    
    if(isign==1) {
      
#ifndef SSEDSLASH_4D_NOCOMMS
      tab->startReceives();
#endif
      
#ifndef SSEDSLASH_4D_NOCOMPUTE
      dispatchToThreads(DslashParscalar32Bit::decomp_plus_3d,
			(void *)psi,
			(void *)chi1,
			(void *)u,
			(void *)s_tab,
			cb,
			subgrid_vol_cb);
#endif
      
#ifndef SSEDSLASH_4D_NOCOMMS
      tab->startSendForward(); 
#endif
      
#ifndef SSEDSLASH_4D_NOCOMPUTE
      dispatchToThreads(DslashParscalar32Bit::decomp_hvv_plus_3d,
			(void*)psi,
			(void*)chi2,
			(void*)u,
			(void*)s_tab,
			cb,
			(int)subgrid_vol_cb);
      
#endif	
      
#ifndef SSEDSLASH_4D_NOCOMMS
      tab->finishSendForward();
      tab->finishReceiveFromBack();
      tab->startSendBack();
#endif   // NOCOMMS
      
#ifndef SSEDSLASH_4D_NOCOMPUTE
      dispatchToThreads(DslashParscalar32Bit::mvv_recons_plus_3d,
			(void*)res,
			(void *)chi1,
			(void *)u,
			(void *)s_tab,
			1-cb,
			subgrid_vol_cb);
#endif
      
#ifndef SSEDSLASH_4D_NOCOMMS
      tab->finishSendBack();
      tab->finishReceiveFromForward();    
#endif  // NOCOMMS
      
      
#ifndef SSEDSLASH_4D_NOCOMPUTE
      dispatchToThreads(DslashParscalar32Bit::recons_plus_3d,
			(void*)res, 
			(void*)chi2,
			(void*)u,	
			(void*)s_tab,
			1-cb,
			subgrid_vol_cb);
#endif
      
    }		
    
    if(isign==-1) {
	
  
#ifndef SSEDSLASH_4D_NOCOMMS
      tab->startReceives();
#endif // NOCOMMS
      
      
#ifndef SSEDSLASH_4D_NOCOMPUTE
      dispatchToThreads(DslashParscalar32Bit::decomp_minus_3d,
			(void*)psi,
			(void *)chi1,
			(void *)u,
			(void *)s_tab,
			cb,
			subgrid_vol_cb);
#endif
      
#ifndef SSEDSLASH_4D_NOCOMMS
      tab->startSendForward();
#endif
      
#ifndef SSEDSLASH_4D_NOCOMPUTE
      dispatchToThreads(DslashParscalar32Bit::decomp_hvv_minus_3d,
			(void*)psi,
			(void*)chi2,
			(void*)u,
			(void*)s_tab,
			cb,
			subgrid_vol_cb);
#endif
      
#ifndef SSEDSLASH_4D_NOCOMMS
      tab->finishSendForward();
      tab->finishReceiveFromBack();
      tab->startSendBack();
#endif
      
#ifndef SSEDSLASH_4D_NOCOMPUTE
      dispatchToThreads(DslashParscalar32Bit::mvv_recons_minus_3d,
			(void*)res,
			(void *)chi1,
			(void *)u,
			(void *)s_tab,
			1-cb,
			subgrid_vol_cb);
#endif
      
#ifndef SSEDSLASH_4D_NOCOMMS
      tab->finishSendBack();
      tab->finishReceiveFromForward();
#endif // #ifndef NOCOMMS
      
#ifndef SSEDSLASH_4D_NOCOMPUTE
      dispatchToThreads(DslashParscalar32Bit::recons_minus_3d,
			(void*)res, 
			(void *)chi2,
			(void *)u,	
			(void *)s_tab,
			1-cb,
			subgrid_vol_cb);
#endif
    }	
  }


  /* INITIALIZE ROUTINE */
  /* Constructor */
  Dslash3D<float>::Dslash3D(const int latt_size[],      
		 void (*getSiteCoords)(int coord[], int node, int linearsite),
		 int (*getLinearSiteIndex)(const int coord[]),
		 int (*getNodeNumber)(const int coord[])
		 ) 
  {
    /* Get the dimensions of the machine */
    const int *machine_size = QMP_get_logical_dimensions();
    
    /* Check we are in 4D */
    if (QMP_get_logical_number_of_dimensions() != 4) {
      QMP_error("init_sse_su3dslash: number of logical dimensions does not match problem");
      QMP_abort(1);
    }

    /* Check problem size in 4D */
    for(int mu=0; mu < 4; mu++)  {
      if ( latt_size[mu] % 2 != 0 ) {
	fprintf(stderr,"This is a Dslash with checkerboarding in 4 dimensions. Each GLOBAL dimension must be even. In addition LOCAL dimension 0 (x) has to be even ,  Your lattice does not meet the GLOBAL requirement latt_size[%d]=%d\n", 
		mu, latt_size[mu]);
	
	exit(1);
      }
    }
    
    /* Check x-checkerboarding */
    if ( (latt_size[0] / machine_size[0]) % 2 != 0 ) {
      fprintf(stderr,"This is a Dslash with checkerboarding in 4 dimensions. Each GLOBAL dimension must be even. In addition LOCAL dimension 0 (x) has to be even ,  Your lattice does not meet the LOCAL requirement\n");
      QMP_abort(1);
    }

    /* Check requisite number of even local sizes */
    int num_even=0;
    int sx = latt_size[0]/machine_size[0]; if( sx%2 == 0 ) num_even++;
    int sy = latt_size[1]/machine_size[1]; if( sy%2 == 0 ) num_even++;
    int sz = latt_size[2]/machine_size[2]; if( sz%2 == 0 ) num_even++;
    int st = latt_size[3]/machine_size[3]; if( st%2 == 0 ) num_even++;

    if( num_even < 2 ) { 
      fprintf(stderr, "Need at least 2 subdimensions to be even");
      QMP_abort(1);
    }

    
    int subgrid[4];
    subgrid[0]=sx;
    subgrid[1]=sy;
    subgrid[2]=sz;
    subgrid[3]=st;

    tab = new DslashTables<HalfSpinor,3>(subgrid);
    s_tab = new ShiftTable3D<HalfSpinor>(subgrid, 
				       tab->getChi1(), 
				       tab->getChi2(), 
				       (HalfSpinor*(*)[3])(tab->getRecvBufptr()), 
				       (HalfSpinor*(*)[3])(tab->getSendBufptr()), 
				       getSiteCoords,
				       getLinearSiteIndex,
				       getNodeNumber);

  }

  /* Destructor */
  Dslash3D<float>::~Dslash3D() 
  {
    delete s_tab;
    delete tab;
  }


} // Namespace
