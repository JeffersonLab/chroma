#include <cpp_clover_parscalar.h>

#include <cstdlib>

#include <cache.h>

#include <cpp_dslash_parscalar_utils_64bit.h>
#include <cpp_clover_odd_odd_parscalar_64bit.h>
#include <tables_parscalar.h>
#include <dispatch_parscalar.h>

#include <cpp_dslash_types.h>
#include <cpp_clover_types.h>
#include "allocate.h"
#include <shift_table_parscalar.h>

#include <qmp.h>
#define QMP_COMMS 

using namespace CPlusPlusWilsonDslash::Dslash64BitTypes;
using namespace CPlusPlusWilsonDslash::DslashParscalar64Bit;

namespace CPlusPlusClover {


  using namespace Clover64BitTypes;
  using namespace CPlusPlusClover64Bit;

  namespace CloverParscalar64Bit {
 



void recons_plus_clov_inv(size_t lo,size_t hi, int id, const void *ptr )	
{
  CloverThreadWorkerArgs *a =(CloverThreadWorkerArgs *)ptr; 
  ShiftTable<HalfSpinor>* tab = (ShiftTable<HalfSpinor>*)a->s;
  int subgrid_vol_cb = tab->subgridVolCB();

  
  int ix1=0;
  // Always on checkerboard 0
  const int low = lo;
  const int high = hi;

  FourSpinor *psi = (FourSpinor*)a->spinor; // Accumulated sum from mvv_recons
  HalfSpinor *chi = (HalfSpinor*)a->half_spinor; // stuff to recons here
  FourSpinor* res =(FourSpinor *) a->spinor2; // result
  CloverTerm* invclov = (CloverTerm *)a->clov; // clover term 

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
    cloverSiteApply(res[thissite], invclov[thissite], *rn);
  }
}

void recons_plus_final(size_t lo,size_t hi, int id, const void *ptr )	
{
  CloverThreadWorkerArgs *a =(CloverThreadWorkerArgs *)ptr; 
  ShiftTable<HalfSpinor>* tab = (ShiftTable<HalfSpinor>*)a->s;
  int subgrid_vol_cb = tab->subgridVolCB();

  // cb is always 1
  int ix1=0;
  const int low = subgrid_vol_cb + lo; 
  const int high = subgrid_vol_cb + hi;

  FourSpinor *psi = (FourSpinor*)a->spinor;
  HalfSpinor *chi = (HalfSpinor*)a->half_spinor;
  FourSpinor* clov_oo_psi = (FourSpinor *)a->spinor2; // A_oo in here 

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
    (*rn)[0][0][0] = clov_oo_psi[thissite][0][0][0] - (*rn)[0][0][0];
    (*rn)[0][0][1] = clov_oo_psi[thissite][0][0][1] - (*rn)[0][0][1];
    (*rn)[0][1][0] = clov_oo_psi[thissite][0][1][0] - (*rn)[0][1][0];
    (*rn)[0][1][1] = clov_oo_psi[thissite][0][1][1] - (*rn)[0][1][1];
    (*rn)[0][2][0] = clov_oo_psi[thissite][0][2][0] - (*rn)[0][2][0];
    (*rn)[0][2][1] = clov_oo_psi[thissite][0][2][1] - (*rn)[0][2][1];

    (*rn)[1][0][0] = clov_oo_psi[thissite][1][0][0] - (*rn)[1][0][0];
    (*rn)[1][0][1] = clov_oo_psi[thissite][1][0][1] - (*rn)[1][0][1];
    (*rn)[1][1][0] = clov_oo_psi[thissite][1][1][0] - (*rn)[1][1][0];
    (*rn)[1][1][1] = clov_oo_psi[thissite][1][1][1] - (*rn)[1][1][1];
    (*rn)[1][2][0] = clov_oo_psi[thissite][1][2][0] - (*rn)[1][2][0];
    (*rn)[1][2][1] = clov_oo_psi[thissite][1][2][1] - (*rn)[1][2][1];

    (*rn)[2][0][0] = clov_oo_psi[thissite][2][0][0] - (*rn)[2][0][0];
    (*rn)[2][0][1] = clov_oo_psi[thissite][2][0][1] - (*rn)[2][0][1];
    (*rn)[2][1][0] = clov_oo_psi[thissite][2][1][0] - (*rn)[2][1][0];
    (*rn)[2][1][1] = clov_oo_psi[thissite][2][1][1] - (*rn)[2][1][1];
    (*rn)[2][2][0] = clov_oo_psi[thissite][2][2][0] - (*rn)[2][2][0];
    (*rn)[2][2][1] = clov_oo_psi[thissite][2][2][1] - (*rn)[2][2][1];

    (*rn)[3][0][0] = clov_oo_psi[thissite][3][0][0] - (*rn)[3][0][0];
    (*rn)[3][0][1] = clov_oo_psi[thissite][3][0][1] - (*rn)[3][0][1];
    (*rn)[3][1][0] = clov_oo_psi[thissite][3][1][0] - (*rn)[3][1][0];
    (*rn)[3][1][1] = clov_oo_psi[thissite][3][1][1] - (*rn)[3][1][1];
    (*rn)[3][2][0] = clov_oo_psi[thissite][3][2][0] - (*rn)[3][2][0];
    (*rn)[3][2][1] = clov_oo_psi[thissite][3][2][1] - (*rn)[3][2][1];
 
  }
}



void recons_minus_clov_inv(size_t lo,size_t hi, int id, const void *ptr )	
{
  CloverThreadWorkerArgs *a =(CloverThreadWorkerArgs *)ptr; 
  ShiftTable<HalfSpinor>* tab = (ShiftTable<HalfSpinor>*)a->s;
  int subgrid_vol_cb = tab->subgridVolCB();

  // cp = 0 always 
  int ix1=0;
  const int low = lo; 
  const int high = hi;

  FourSpinor *psi = (FourSpinor*)a->spinor;
  HalfSpinor *chi = (HalfSpinor*)a->half_spinor;
  FourSpinor* res =(FourSpinor *) a->spinor2; // result
  CloverTerm* invclov = (CloverTerm *)a->clov; // clover term 

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
    cloverSiteApply(res[thissite], invclov[thissite], *rn);
    /* end of loop */
  }
}

void recons_minus_final(size_t lo,size_t hi, int id, const void *ptr )	
{
  CloverThreadWorkerArgs *a =(CloverThreadWorkerArgs *)ptr; 
  ShiftTable<HalfSpinor>* tab = (ShiftTable<HalfSpinor>*)a->s;
  int subgrid_vol_cb = tab->subgridVolCB();

  // CB = 1 always
  int ix1=0;
  const int low = subgrid_vol_cb + lo; 
  const int high = subgrid_vol_cb + hi;

  FourSpinor *psi = (FourSpinor*)a->spinor;
  HalfSpinor *chi = (HalfSpinor*)a->half_spinor;
  FourSpinor* clov_oo_psi = (FourSpinor *)a->spinor2;
 
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
    (*rn)[0][0][0] = clov_oo_psi[thissite][0][0][0] - (*rn)[0][0][0];
    (*rn)[0][0][1] = clov_oo_psi[thissite][0][0][1] - (*rn)[0][0][1];
    (*rn)[0][1][0] = clov_oo_psi[thissite][0][1][0] - (*rn)[0][1][0];
    (*rn)[0][1][1] = clov_oo_psi[thissite][0][1][1] - (*rn)[0][1][1];
    (*rn)[0][2][0] = clov_oo_psi[thissite][0][2][0] - (*rn)[0][2][0];
    (*rn)[0][2][1] = clov_oo_psi[thissite][0][2][1] - (*rn)[0][2][1];

    (*rn)[1][0][0] = clov_oo_psi[thissite][1][0][0] - (*rn)[1][0][0];
    (*rn)[1][0][1] = clov_oo_psi[thissite][1][0][1] - (*rn)[1][0][1];
    (*rn)[1][1][0] = clov_oo_psi[thissite][1][1][0] - (*rn)[1][1][0];
    (*rn)[1][1][1] = clov_oo_psi[thissite][1][1][1] - (*rn)[1][1][1];
    (*rn)[1][2][0] = clov_oo_psi[thissite][1][2][0] - (*rn)[1][2][0];
    (*rn)[1][2][1] = clov_oo_psi[thissite][1][2][1] - (*rn)[1][2][1];

    (*rn)[2][0][0] = clov_oo_psi[thissite][2][0][0] - (*rn)[2][0][0];
    (*rn)[2][0][1] = clov_oo_psi[thissite][2][0][1] - (*rn)[2][0][1];
    (*rn)[2][1][0] = clov_oo_psi[thissite][2][1][0] - (*rn)[2][1][0];
    (*rn)[2][1][1] = clov_oo_psi[thissite][2][1][1] - (*rn)[2][1][1];
    (*rn)[2][2][0] = clov_oo_psi[thissite][2][2][0] - (*rn)[2][2][0];
    (*rn)[2][2][1] = clov_oo_psi[thissite][2][2][1] - (*rn)[2][2][1];

    (*rn)[3][0][0] = clov_oo_psi[thissite][3][0][0] - (*rn)[3][0][0];
    (*rn)[3][0][1] = clov_oo_psi[thissite][3][0][1] - (*rn)[3][0][1];
    (*rn)[3][1][0] = clov_oo_psi[thissite][3][1][0] - (*rn)[3][1][0];
    (*rn)[3][1][1] = clov_oo_psi[thissite][3][1][1] - (*rn)[3][1][1];
    (*rn)[3][2][0] = clov_oo_psi[thissite][3][2][0] - (*rn)[3][2][0];
    (*rn)[3][2][1] = clov_oo_psi[thissite][3][2][1] - (*rn)[3][2][1];

    /* end of loop */
  }
}

/*****************end of isign corresponding to -1 **********************/

  } // Namespace 



// Your actual operator
void CloverSchur4D<double>::operator()(double* res, 
				      const double* psi, 
				      const double* u, 
				      const double* clov_oo,
				      const double* invclov_ee,
				      int isign)
{
  HalfSpinor* chi1 = tab->getChi1();
  HalfSpinor* chi2 = tab->getChi2();
  int subgrid_vol_cb = s_tab->subgridVolCB();


  if(isign==1) {

    // Apply Dslash (source cb = 1)
    tab->startReceives();

    CPlusPlusWilsonDslash::dispatchToThreads(DslashParscalar64Bit::decomp_plus,
		      (void *)psi,
		      (void *)chi1,
		      (void *)u,
		      (void *)s_tab,
		      1,
		      subgrid_vol_cb);

    tab->startSendForward(); 

    CPlusPlusWilsonDslash::dispatchToThreads(DslashParscalar64Bit::decomp_hvv_plus,
		      (void*)psi,
		      (void*)chi2,
		      (void*)u,
		      (void*)s_tab,
		      1,
		      (int)subgrid_vol_cb);

    tab->startSendBack();

    // All comms started do the clover odd odd bit...
    CPlusPlusClover::dispatchToThreads(CPlusPlusClover::CPlusPlusClover64Bit::clovOddOddApply,
				       (void *)psi,   
				       (void *)0, // No Half Spinor
				       (void *)0, // No Gauge Field
				       (void *)clov_oo,
				       (void *)t_spinor2, // Only on odd part
				       (void *)s_tab,
				       0, // First half
				       (int)subgrid_vol_cb/2);
				     
				       
    

    tab->finishSendForward();
    tab->finishReceiveFromBack();


    CPlusPlusWilsonDslash::dispatchToThreads(DslashParscalar64Bit::mvv_recons_plus,
					     (void*)t_spinor2,  // Even part
					     (void *)chi1,
					     (void *)u,
					     (void *)s_tab,
					     0,
					     subgrid_vol_cb);

    tab->finishSendBack();
    tab->finishReceiveFromForward();    



    CPlusPlusClover::dispatchToThreads(CloverParscalar64Bit::recons_plus_clov_inv,
				       (void*)t_spinor2,
				       (void*)chi2,
				       (void*)u,	
				       (void *)invclov_ee,
				       (void *)t_spinor, // Result A^{-1} D in even part
				       (void*)s_tab, // 
				       0, // Half is irrelevant
				       subgrid_vol_cb);




    // Apply Dslash (source cb = 0)
    tab->startReceives();

    CPlusPlusWilsonDslash::dispatchToThreads(DslashParscalar64Bit::decomp_plus,
					     (void *)t_spinor,         // Decompose on even part
					     (void *)chi1,
					     (void *)u,
					     (void *)s_tab,
					     0,
					     subgrid_vol_cb);

    tab->startSendForward(); 

    CPlusPlusWilsonDslash::dispatchToThreads(DslashParscalar64Bit::decomp_hvv_plus,
					     (void*)t_spinor,         // Decompose on even part
					     (void*)chi2,
					     (void*)u,
					     (void*)s_tab,
					     0,
					     (int)subgrid_vol_cb);

    tab->startSendBack();


    // All comms started 
    // All comms started do the clover odd odd bit...
    CPlusPlusClover::dispatchToThreads(CPlusPlusClover::CPlusPlusClover64Bit::clovOddOddApply,
				       (void *)psi,   
				       (void *)0, // No Half Spinor
				       (void *)0, // No Gauge Field
				       (void *)clov_oo, 
				       (void *)t_spinor2, // Even part  (t_spinor odd now holds A_oo psi
				       (void *)s_tab,
				       1, // Second  half
				       (int)subgrid_vol_cb/2);

				       
				     


    tab->finishSendForward();
    tab->finishReceiveFromBack();


    CPlusPlusWilsonDslash::dispatchToThreads(DslashParscalar64Bit::mvv_recons_plus,
					     (void*)res,
					     (void *)chi1,
					     (void *)u,
					     (void *)s_tab,
					     1,
					     subgrid_vol_cb);

    tab->finishSendBack();
    tab->finishReceiveFromForward();    

    CPlusPlusClover::dispatchToThreads(CloverParscalar64Bit::recons_plus_final,
				       (void*)res, 
				       (void*)chi2,
				       (void*)0,  // No Gauge Field	
				       (void *)0, // No clover piece
				       (void *)t_spinor2, // No spinor 2
				       (void*)s_tab, // 
				       0, // Half is irrelevant
				       subgrid_vol_cb);



  }		

  if(isign==-1) 
  {

    // Apply Dslash (source cb = 1)
    tab->startReceives();

    CPlusPlusWilsonDslash::dispatchToThreads(DslashParscalar64Bit::decomp_minus,
		      (void *)psi,
		      (void *)chi1,
		      (void *)u,
		      (void *)s_tab,
		      1,
		      subgrid_vol_cb);

    tab->startSendForward(); 

    CPlusPlusWilsonDslash::dispatchToThreads(DslashParscalar64Bit::decomp_hvv_minus,
		      (void*)psi,
		      (void*)chi2,
		      (void*)u,
		      (void*)s_tab,
		      1,
		      (int)subgrid_vol_cb);

    tab->startSendBack();


    // All comms started do the clover odd odd bit...
    CPlusPlusClover::dispatchToThreads(CPlusPlusClover::CPlusPlusClover64Bit::clovOddOddApply,
				       (void *)psi,   
				       (void *)0, // No Half Spinor
				       (void *)0, // No Gauge Field
				       (void *)clov_oo, 
				       (void *)t_spinor2,   // odd odd first half
				       (void *)s_tab,
				       0, // First half
				       (int)subgrid_vol_cb/2);

				       
				     
				       
    

    tab->finishSendForward();
    tab->finishReceiveFromBack();


    CPlusPlusWilsonDslash::dispatchToThreads(DslashParscalar64Bit::mvv_recons_minus,
					     (void*)t_spinor2, // even half
					     (void *)chi1,
					     (void *)u,
					     (void *)s_tab,
					     0,
					     subgrid_vol_cb);

    tab->finishSendBack();
    tab->finishReceiveFromForward();    

    CPlusPlusClover::dispatchToThreads(CloverParscalar64Bit::recons_minus_clov_inv,
				       (void*)t_spinor2, // even half 
				       (void*)chi2,
				       (void*)u,	
				       (void *)invclov_ee,
				       (void *)t_spinor, // No spinor 2
				       (void*)s_tab, // 
				       0, // Half is irrelevant
				       subgrid_vol_cb);


    // Apply Dslash (source cb = 0)
    tab->startReceives();

    CPlusPlusWilsonDslash::dispatchToThreads(DslashParscalar64Bit::decomp_minus,
		      (void *)t_spinor,
		      (void *)chi1,
		      (void *)u,
		      (void *)s_tab,
		      0,
		      subgrid_vol_cb);

    tab->startSendForward(); 

    CPlusPlusWilsonDslash::dispatchToThreads(DslashParscalar64Bit::decomp_hvv_minus,
		      (void*)t_spinor,
		      (void*)chi2,
		      (void*)u,
		      (void*)s_tab,
		      0,
		      (int)subgrid_vol_cb);

    tab->startSendBack();


    // All comms started 
    // All comms started do the clover odd odd bit...
    CPlusPlusClover::dispatchToThreads(CPlusPlusClover::CPlusPlusClover64Bit::clovOddOddApply,
				       (void *)psi,   
				       (void *)0, // No Half Spinor
				       (void *)0, // No Gauge Field
				       (void *)clov_oo, 
				       (void *)t_spinor2,    // Odd part 
				       (void *)s_tab,
				       1, // Second  half
				       (int)subgrid_vol_cb/2);

				       
				     


    tab->finishSendForward();
    tab->finishReceiveFromBack();


    CPlusPlusWilsonDslash::dispatchToThreads(DslashParscalar64Bit::mvv_recons_minus,
					     (void*)res,
					     (void *)chi1,
					     (void *)u,
					     (void *)s_tab,
					     1,
					     subgrid_vol_cb);

    tab->finishSendBack();
    tab->finishReceiveFromForward();    

    CPlusPlusClover::dispatchToThreads(CloverParscalar64Bit::recons_minus_final,
		      (void*)res, 
		      (void*)chi2,
				       (void*)0,  // No Gauge Field	
				       (void *)0, // No clover piece
		      (void *)t_spinor2, // No spinor 2
		      (void*)s_tab, // 
		      0, // Half is irrelevant
		      subgrid_vol_cb);

  }		
}







  /* INITIALIZE ROUTINE */
  /* Constructor */
CloverSchur4D<double>::CloverSchur4D(const int latt_size[],      
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

    tab = new DslashTables<HalfSpinor,4>(subgrid);
    s_tab = new ShiftTable<HalfSpinor>(subgrid, 
				       tab->getChi1(), 
				       tab->getChi2(), 
				       (HalfSpinor*(*)[4])(tab->getRecvBufptr()), 
				       (HalfSpinor*(*)[4])(tab->getSendBufptr()), 
				       getSiteCoords,
				       getLinearSiteIndex,
				       getNodeNumber);

    xt_spinor = (FourSpinor *)CPlusPlusWilsonDslash::alloc(2*s_tab->subgridVolCB()*sizeof(FourSpinor)
				     +Cache::CacheLineSize);
    if( xt_spinor == (FourSpinor *)NULL ) { 
      QMP_error("Couldnt allocate temporary xt_spinor");
      exit(1);
    }
    unsigned long pad = 0;
    if ( (unsigned long)xt_spinor % Cache::CacheLineSize != 0 ) { 
      pad = Cache::CacheLineSize - (unsigned long)xt_spinor % Cache::CacheLineSize;
    }
    t_spinor = (FourSpinor *)((unsigned char *)xt_spinor + pad);

    xt_spinor2 = (FourSpinor *)CPlusPlusWilsonDslash::alloc(2*s_tab->subgridVolCB()*sizeof(FourSpinor)
				     +Cache::CacheLineSize);
    if( xt_spinor2 == (FourSpinor *)NULL ) { 
      QMP_error("Couldnt allocate temporary xt_spinor");
      exit(1);
    }
    
    pad = 0;
    if ( (unsigned long)xt_spinor2 % Cache::CacheLineSize != 0 ) { 
      pad = Cache::CacheLineSize - (unsigned long)xt_spinor2 % Cache::CacheLineSize;
    }
    t_spinor2 = (FourSpinor *)((unsigned char *)xt_spinor2 + pad);

  }

  /* Destructor */
  CloverSchur4D<double>::~CloverSchur4D() 
  {
    delete s_tab;
    delete tab;
    CPlusPlusWilsonDslash::dealloc(xt_spinor);
    CPlusPlusWilsonDslash::dealloc(xt_spinor2);

  }


} // Namespace
