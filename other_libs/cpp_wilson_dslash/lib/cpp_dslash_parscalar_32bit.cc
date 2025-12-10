#include <cpp_dslash_parscalar.h>

#include <cstdlib>

#include <cache.h>
#include <cpp_dslash_parscalar_utils_32bit.h>
#include <tables_parscalar.h>
#include <dispatch_parscalar.h>
#include <cpp_dslash_types.h>
#include <shift_table_parscalar.h>

using namespace CPlusPlusWilsonDslash::Dslash32BitTypes;

#define QMP_COMMS 

namespace CPlusPlusWilsonDslash { 

// Temporary space...
ShiftTable<Dslash<float>::HalfSpinor>* Dslash<float>::s_tab = nullptr;
DslashTables<Dslash<float>::HalfSpinor, 4>* Dslash<float>::tab = nullptr;

// Your actual operator
void Dslash<float>::operator()(float* res, 
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
    dispatchToThreads(DslashParscalar32Bit::decomp_plus,
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
    dispatchToThreads(DslashParscalar32Bit::decomp_hvv_plus,
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
    dispatchToThreads(DslashParscalar32Bit::mvv_recons_plus,
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
    dispatchToThreads(DslashParscalar32Bit::recons_plus,
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
    dispatchToThreads(DslashParscalar32Bit::decomp_minus,
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
    dispatchToThreads(DslashParscalar32Bit::decomp_hvv_minus,
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
    dispatchToThreads(DslashParscalar32Bit::mvv_recons_minus,
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
    dispatchToThreads(DslashParscalar32Bit::recons_minus,
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
  Dslash<float>::Dslash(const int latt_size[],      
		 void (*getSiteCoords)(int coord[], int node, int linearsite),
		 int (*getLinearSiteIndex)(const int coord[]),
		 int (*getNodeNumber)(const int coord[])
		 ) 
  {
	  if ( tab == nullptr && s_tab == nullptr ) {
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


    tab = new DslashTables<HalfSpinor,4>(subgrid);


    s_tab = new ShiftTable<HalfSpinor>(subgrid, 
				       tab->getChi1(), 
				       tab->getChi2(), 
				       (HalfSpinor*(*)[4])(tab->getRecvBufptr()), 
				       (HalfSpinor*(*)[4])(tab->getSendBufptr()), 
				       getSiteCoords,
				       getLinearSiteIndex,
				       getNodeNumber);
    }
  }

  /* Destructor */
  Dslash<float>::~Dslash() 
  {
	  // Never free
#if 0
    delete s_tab;
    delete tab;
#endif
  }


} // Namespace
