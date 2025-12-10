#include <cpp_dslash_parscalar.h>

#include <cstdlib>

#include <cache.h>

#include <cpp_dslash_parscalar_utils_64bit.h>

#include <tables_parscalar.h>
#include <dispatch_parscalar.h>

#include <cpp_dslash_types.h>
using namespace CPlusPlusWilsonDslash::Dslash64BitTypes;

#include <shift_table_3d_parscalar.h>


namespace CPlusPlusWilsonDslash {
ShiftTable3D<Dslash<double>::HalfSpinor>* Dslash3D<double>::s_tab = nullptr;
DslashTables<Dslash<double>::HalfSpinor, 3>* Dslash3D<double>::tab = nullptr;
  namespace DslashParscalar64Bit {
    /* this routine is similar to wnxtsu3dslash, except instead of handling the second site's worth in the same loop, the second spin component's worth must be handled seperately */
    void decomp_plus_3d(size_t lo, size_t hi, int id, const void *ptr)
    {
      const ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr; 
      ShiftTable3D<HalfSpinor>* stab = (ShiftTable3D<HalfSpinor>*)a->s;
      int subgrid_vol_cb = stab->subgridVolCB();
      const int cb = a->cb; 
      FourSpinor *psi ALIGN = (FourSpinor*)a->spinor;
      HalfSpinor *chi ALIGN = (HalfSpinor*)a->half_spinor;
      
      int ix1;
      HalfSpinor *s3 ALIGN ;
      FourSpinor* sp ALIGN;
      const int low = cb*subgrid_vol_cb + lo; 
      const int high = cb*subgrid_vol_cb + hi;

      
      for (ix1=low;ix1<high;++ix1) {
	int thissite = stab->siteTable( ix1 );
	sp=&psi[thissite];
	s3 = stab->halfspinorBufferOffset(DECOMP_SCATTER,ix1,0);
	decomp_gamma0_minus(*sp, *s3);
	/******************************* direction +1 *********************************/
	s3 = stab->halfspinorBufferOffset(DECOMP_SCATTER,ix1,1);
	decomp_gamma1_minus(*sp, *s3);
	
	/******************************* direction +2 *********************************/
	s3 = stab->halfspinorBufferOffset(DECOMP_SCATTER,ix1,2);
	decomp_gamma2_minus(*sp, *s3);
	
      }
    }
    
    
    void decomp_hvv_plus_3d(size_t lo, size_t hi, int id, const void *ptr)
    {
      const ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr; 
      ShiftTable3D<HalfSpinor>* stab = (ShiftTable3D<HalfSpinor>*)a->s;
      int subgrid_vol_cb = stab->subgridVolCB();

      const int cb = a->cb; 
      GaugeMatrix (*gauge_field)[4] = (GaugeMatrix (*)[4])a->u;
      FourSpinor *psi = (FourSpinor*)a->spinor;
      HalfSpinor *chi = (HalfSpinor*)a->half_spinor;
      
      int ix1=0;
      
      GaugeMatrix *um ALIGN;

      HalfSpinor *s3 ALIGN;
      
      FourSpinor *sm ALIGN; 
      
      
      
      const int low = cb*subgrid_vol_cb + lo; 
      const int high = cb*subgrid_vol_cb + hi;
      
      
      for (ix1=low;ix1<high;++ix1) {
	int thissite = stab->siteTable( ix1 );
	
	sm=&psi[thissite];
	
	/******************************* direction -0 *********************************/
	um=&gauge_field[thissite][0];
	s3 = stab->halfspinorBufferOffset(DECOMP_HVV_SCATTER,ix1,0);
	decomp_hvv_gamma0_plus(*sm, *um, *s3);
	
	/******************************* direction -1 *********************************/
	um=&gauge_field[thissite][1];
	s3 = stab->halfspinorBufferOffset(DECOMP_HVV_SCATTER,ix1,1);
	decomp_hvv_gamma1_plus(*sm, *um, *s3);
	
	/******************************* direction -2 *********************************/
	um=&gauge_field[thissite][2];
	s3 = stab->halfspinorBufferOffset(DECOMP_HVV_SCATTER,ix1,2);
	decomp_hvv_gamma2_plus(*sm, *um, *s3);
      }
    }
    
    
    void mvv_recons_plus_3d(size_t lo, size_t hi, int id, const void *ptr)
    {
      const ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr; 
      ShiftTable3D<HalfSpinor>* stab = (ShiftTable3D<HalfSpinor>*)a->s;
      int subgrid_vol_cb = stab->subgridVolCB();

      const int cb = a->cb; 
      int ix1=0;
      
      GaugeMatrix (*gauge_field)[4] = (GaugeMatrix (*)[4])a->u;
      FourSpinor *psi = (FourSpinor*)a->spinor;
      HalfSpinor *chi = (HalfSpinor*)a->half_spinor;
      GaugeMatrix *up ALIGN;
      HalfSpinor *s3 ALIGN, *s4 ALIGN;
      FourSpinor part_sum ALIGN, *result ALIGN;
      
      

      const int low = cb*subgrid_vol_cb + lo; 
      const int high = cb*subgrid_vol_cb + hi;
      
      
      for (ix1=low;ix1<high;++ix1) {
	int thissite = stab->siteTable( ix1 );
	
	
	/******************************* direction +0 *********************************/	
	up=&gauge_field[thissite][0];    
	s3 = stab->halfspinorBufferOffset(RECONS_MVV_GATHER,ix1,0);
	mvv_recons_gamma0_plus(*s3, *up, part_sum);
	
	
	/******************************* direction +1 *********************************/
	up=&gauge_field[thissite][1];
	s3 = stab->halfspinorBufferOffset(RECONS_MVV_GATHER,ix1,1);
	mvv_recons_gamma1_plus_add(*s3, *up, part_sum);
	
	
	/******************************* direction +2 *********************************/
	up=&gauge_field[thissite][2];
	s3 = stab->halfspinorBufferOffset(RECONS_MVV_GATHER,ix1,2);
	mvv_recons_gamma2_plus_add_store(*s3, *up, part_sum, psi[thissite]);
      }
    }
    

    
    /*optimized for SZIN spin basis */
    void recons_plus_3d(size_t lo, size_t hi, int id, const void *ptr)
    {
      
      
      const ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr; 
      ShiftTable3D<HalfSpinor>* stab = (ShiftTable3D<HalfSpinor>*)a->s;
      int subgrid_vol_cb = stab->subgridVolCB();

      const int cb = a->cb; 
      int ix1=0;
      
      FourSpinor *psi = (FourSpinor*)a->spinor;
      HalfSpinor *chi = (HalfSpinor*)a->half_spinor;
      
      HalfSpinor *hs0 ALIGN , *hs1 ALIGN, *hs2 ALIGN;
      FourSpinor *s1 ALIGN, *rn ALIGN;
      
      int low = cb*subgrid_vol_cb + lo;
      int high = cb*subgrid_vol_cb + hi;
      
      
      /************************ loop over all lattice sites *************************/
      for (ix1 =low; ix1 <high ;ix1++) {
	
	int thissite = stab->siteTable( ix1 );
	
	/* first spin component of result */
	hs0 = stab->halfspinorBufferOffset(RECONS_GATHER,ix1,0);
	hs1 = stab->halfspinorBufferOffset(RECONS_GATHER,ix1,1);	  
	hs2 = stab->halfspinorBufferOffset(RECONS_GATHER,ix1,2);
	recons_3dir_plus(*hs0, *hs1, *hs2, psi[thissite]);
      }
      
    }
    
    
    /************now for isign = -1  **********************/
    
    
    void decomp_minus_3d(size_t lo, size_t hi, int id, const void *ptr)
    {
      
      
      const ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr; 
      ShiftTable3D<HalfSpinor>* stab = (ShiftTable3D<HalfSpinor>*)a->s;
      int subgrid_vol_cb = stab->subgridVolCB();

      const int cb = a->cb; 
      int ix1=0;
      
      FourSpinor *psi = (FourSpinor*)a->spinor;
      HalfSpinor *chi = (HalfSpinor*)a->half_spinor;
      
      HalfSpinor *s3 ALIGN; 
      FourSpinor *sp ALIGN;
      
      int low = cb*subgrid_vol_cb + lo;
      int high = cb*subgrid_vol_cb + hi;
      
      
      /************************ loop over all lattice sites *************************/
      for (ix1 =low; ix1 <high ;ix1++) {
	
	int thissite = stab->siteTable( ix1 );       
	sp=&psi[thissite];
	
	s3 = stab->halfspinorBufferOffset(DECOMP_SCATTER,ix1,0);
	decomp_gamma0_plus(*sp, *s3);
	
	s3 = stab->halfspinorBufferOffset(DECOMP_SCATTER,ix1,1);
	decomp_gamma1_plus(*sp, *s3);
	
	s3 = stab->halfspinorBufferOffset(DECOMP_SCATTER,ix1,2);
	decomp_gamma2_plus(*sp, *s3);
	
      }
    }
    
    
    void decomp_hvv_minus_3d(size_t lo, size_t hi, int id, const void *ptr)
    {
      
      const ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr; 

      ShiftTable3D<HalfSpinor>* stab = (ShiftTable3D<HalfSpinor>*)a->s;
      int subgrid_vol_cb = stab->subgridVolCB();

      const int cb = a->cb; 
      int ix1=0;
      
      
      GaugeMatrix (*gauge_field)[4] = (GaugeMatrix (*)[4])a->u;
      
      FourSpinor *psi = (FourSpinor*)a->spinor;
      HalfSpinor *chi = (HalfSpinor*)a->half_spinor;
      
      GaugeMatrix *um ALIGN;
      HalfSpinor *s3 ALIGN;
      FourSpinor *s1 ALIGN, *sm ALIGN;
      
      
      int low = cb*subgrid_vol_cb + lo;
      int high = cb*subgrid_vol_cb + hi;
      
      for (ix1 =low; ix1 <high ;ix1++) {
	
	int thissite = stab->siteTable( ix1 );    
	
	sm=&psi[thissite];
	um=&gauge_field[thissite][0];
	
	s3 = stab->halfspinorBufferOffset(DECOMP_HVV_SCATTER,ix1,0);
	decomp_hvv_gamma0_minus(*sm, *um, *s3);	   
	
	um=&gauge_field[thissite][1];
	s3 = stab->halfspinorBufferOffset(DECOMP_HVV_SCATTER,ix1,1);
	decomp_hvv_gamma1_minus(*sm, *um, *s3);	   
	
	um=&gauge_field[thissite][2];
	s3 = stab->halfspinorBufferOffset(DECOMP_HVV_SCATTER,ix1,2);
	decomp_hvv_gamma2_minus(*sm, *um, *s3);	   
	
      }
    }
    
    
    void mvv_recons_minus_3d(size_t lo, size_t hi, int id, const void *ptr)
    {
      
      
      const ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr; 
      ShiftTable3D<HalfSpinor>* stab = (ShiftTable3D<HalfSpinor>*)a->s;
      int subgrid_vol_cb = stab->subgridVolCB();

      const int cb = a->cb; 
      int ix1=0;
      
      GaugeMatrix (*gauge_field)[4] = (GaugeMatrix (*)[4])a->u;
      FourSpinor *psi = (FourSpinor*)a->spinor;
      HalfSpinor *chi = (HalfSpinor*)a->half_spinor;
      
      GaugeMatrix *up ALIGN;
      
      HalfSpinor *s3 ALIGN, *s4 ALIGN;
      
      FourSpinor rs ALIGN,*rn ALIGN;
      
      
      int low = cb*subgrid_vol_cb + lo;
      int high = cb*subgrid_vol_cb + hi;
      
      for (ix1 =low; ix1 <high ;ix1++) {
	int thissite = stab->siteTable( ix1 );       
	
	up=&gauge_field[thissite][0];
	s3 = stab->halfspinorBufferOffset(RECONS_MVV_GATHER,ix1,0);
	mvv_recons_gamma0_minus(*s3, *up, rs);
	
	up = &gauge_field[thissite][1];
	s3 = stab->halfspinorBufferOffset(RECONS_MVV_GATHER,ix1,1);
	mvv_recons_gamma1_minus_add(*s3, *up, rs);
	
	up = &gauge_field[thissite][2];
	s3 = stab->halfspinorBufferOffset(RECONS_MVV_GATHER,ix1,2);
	mvv_recons_gamma2_minus_add_store(*s3, *up, rs, psi[thissite]);
      }
      
    }
    
    
    
    void recons_minus_3d(size_t lo, size_t hi, int id, const void *ptr)
    {
      
      const ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr; 
      ShiftTable3D<HalfSpinor>* stab = (ShiftTable3D<HalfSpinor>*)a->s;
      int subgrid_vol_cb = stab->subgridVolCB();

      const int cb = a->cb; 
      int ix1=0;
      
      
      FourSpinor *psi = (FourSpinor*)a->spinor;
      HalfSpinor *chi = (HalfSpinor*)a->half_spinor;
      
      HalfSpinor *hs0,*hs1,*hs2,*hs3;   
      FourSpinor  *s1 ALIGN,  *rn ALIGN;
      
      int low = cb*subgrid_vol_cb + lo;
      int high = cb*subgrid_vol_cb + hi;
      
      
      /************************ loop over all lattice sites *************************/
      for (ix1 =low; ix1 <high ;ix1++) {
	
	int thissite = stab->siteTable( ix1 );     
	
	rn=&psi[thissite];
	/* first spin component of result */
	hs0 = stab->halfspinorBufferOffset(RECONS_GATHER,ix1,0);
	hs1 = stab->halfspinorBufferOffset(RECONS_GATHER,ix1,1);
	hs2 = stab->halfspinorBufferOffset(RECONS_GATHER,ix1,2);
	recons_3dir_minus(*hs0, *hs1, *hs2, psi[thissite]);
	/* end of loop */
      }
    }
    
  }




// Your actual operator
void Dslash3D<double>::operator()(double* res, 
			double* psi, 
			double* u, 
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
    dispatchToThreads(DslashParscalar64Bit::decomp_plus_3d,
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
    dispatchToThreads(DslashParscalar64Bit::decomp_hvv_plus_3d,
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
    dispatchToThreads(DslashParscalar64Bit::mvv_recons_plus_3d,
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
    dispatchToThreads(DslashParscalar64Bit::recons_plus_3d,
		      (void*)res, 
		      (void*)chi2,
		      (void*)u,	
		      (void*)s_tab,
		      1-cb,
		      subgrid_vol_cb);
#endif

  }		

  if(isign==-1) 
  {

  
#ifndef SSEDSLASH_4D_NOCOMMS
    tab->startReceives();
#endif // NOCOMMS


#ifndef SSEDSLASH_4D_NOCOMPUTE
    dispatchToThreads(DslashParscalar64Bit::decomp_minus_3d,
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
    dispatchToThreads(DslashParscalar64Bit::decomp_hvv_minus_3d,
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
    dispatchToThreads(DslashParscalar64Bit::mvv_recons_minus_3d,
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
    dispatchToThreads(DslashParscalar64Bit::recons_minus_3d,
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
  Dslash3D<double>::Dslash3D(const int latt_size[],      
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
	fprintf(stderr,"This is a Dslash3D with checkerboarding in 4 dimensions. Each GLOBAL dimension must be even. In addition LOCAL dimension 0 (x) has to be even ,  Your lattice does not meet the GLOBAL requirement latt_size[%d]=%d\n", 
		mu, latt_size[mu]);
	
	exit(1);
      }
    }
    
    /* Check x-checkerboarding */
    if ( (latt_size[0] / machine_size[0]) % 2 != 0 ) {
      fprintf(stderr,"This is a Dslash3D with checkerboarding in 4 dimensions. Each GLOBAL dimension must be even. In addition LOCAL dimension 0 (x) has to be even ,  Your lattice does not meet the LOCAL requirement\n");
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
  Dslash3D<double>::~Dslash3D() 
  {
    delete s_tab;
    delete tab;
  }


} // Namespace
