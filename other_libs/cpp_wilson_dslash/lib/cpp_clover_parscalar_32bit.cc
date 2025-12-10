#include <cpp_clover_parscalar.h>

#include <cstdlib>

#include <cache.h>

#include <cpp_dslash_parscalar_utils_32bit.h>
#include <cpp_clover_odd_odd_parscalar_32bit.h>
#include <tables_parscalar.h>
#include <dispatch_parscalar.h>

#include <cpp_dslash_types.h>
#include <cpp_clover_types.h>

#include <shift_table_parscalar.h>
#include "allocate.h"
#include <qmp.h>
#define QMP_COMMS 

using namespace CPlusPlusWilsonDslash::Dslash32BitTypes;
using namespace CPlusPlusWilsonDslash::DslashParscalar32Bit;

namespace CPlusPlusClover {


  using namespace Clover32BitTypes;
  using namespace CPlusPlusClover32Bit;

  namespace CloverParscalar32Bit {
 



void recons_plus_clov_inv(size_t lo,size_t hi, int id, const void *ptr )	
{
  int ix1;
  FourSpinor* sn1 ALIGN;
  FourSpinor* sn2 ALIGN;
  
  FourSpinor tmp ALIGN;

  const CloverThreadWorkerArgs *a = (CloverThreadWorkerArgs *)ptr;
  ShiftTable<HalfSpinor>* s=(ShiftTable<HalfSpinor>*)a->s;
  int subgrid_vol_cb=s->subgridVolCB();

  // This has a half accumulated sum in it... 
  FourSpinor* spinor_field = (FourSpinor *)a->spinor;
  FourSpinor* res =(FourSpinor *) a->spinor2;

  HalfSpinor* chi = (HalfSpinor *)a->half_spinor;
  CloverTerm* invclov = (CloverTerm *)a->clov;

  // Even even is always cb = 0;

  HalfSpinor *hs0 ALIGN;
  HalfSpinor *hs1 ALIGN;
  HalfSpinor *hs2 ALIGN;
  HalfSpinor *hs3 ALIGN;

  int low = lo;
  int high = hi;
  
  const int PREFDIST=4;

  /************************ loop over all lattice sites *************************/
  int thissite = s->siteTable(low);
  int clovsite = thissite;

  hs0 =  s->halfspinorBufferOffset(RECONS_GATHER,low,0); 
  _mm_prefetch(hs0, _MM_HINT_NTA);

  hs1 =  s->halfspinorBufferOffset(RECONS_GATHER,low,1); 
  _mm_prefetch(hs1, _MM_HINT_NTA);

  hs2 =  s->halfspinorBufferOffset(RECONS_GATHER,low,2); 
  _mm_prefetch(hs2, _MM_HINT_NTA);

  hs3 =  s->halfspinorBufferOffset(RECONS_GATHER,low,3);
  _mm_prefetch(hs3, _MM_HINT_NTA);

  sn1=&spinor_field[thissite];   
  _mm_prefetch(sn1, _MM_HINT_NTA);

  for (ix1=low+1;ix1<high;ix1++) {
    thissite = s->siteTable(ix1);
    sn2 = &spinor_field[thissite];   
    _mm_prefetch(sn2, _MM_HINT_NTA);
   
    // Reconstruct into temporary...
    recons_4dir_plus(*hs0, *hs1, *hs2, *hs3, *sn1);
    cloverSiteApply(res[clovsite], invclov[clovsite], *sn1);


    hs0 =  s->halfspinorBufferOffset(RECONS_GATHER,ix1,0); 
    _mm_prefetch(hs0, _MM_HINT_NTA);

    hs1 =  s->halfspinorBufferOffset(RECONS_GATHER,ix1,1); 
    _mm_prefetch(hs1, _MM_HINT_NTA);

    hs2 =  s->halfspinorBufferOffset(RECONS_GATHER,ix1,2); 
    _mm_prefetch(hs2, _MM_HINT_NTA);

    hs3 =  s->halfspinorBufferOffset(RECONS_GATHER,ix1,3); 
    _mm_prefetch(hs3, _MM_HINT_NTA);

    sn1=sn2;
    clovsite = thissite;

    /*************************end of loop ****************************/
  }
  recons_4dir_plus(*hs0, *hs1, *hs2, *hs3, *sn1);
  cloverSiteApply(res[clovsite], invclov[clovsite], *sn1);

}

void recons_plus_final(size_t lo,size_t hi, int id, const void *ptr )	
{
  int ix1;
  FourSpinor* sn1 ALIGN;
  FourSpinor* sn2 ALIGN;
  
  FourSpinor tmp ALIGN;

  const CloverThreadWorkerArgs *a = (CloverThreadWorkerArgs *)ptr;
  ShiftTable<HalfSpinor>* s=(ShiftTable<HalfSpinor>*)a->s;
  int subgrid_vol_cb=s->subgridVolCB();

  FourSpinor* spinor_field = (FourSpinor *)a->spinor;
  HalfSpinor* chi = (HalfSpinor *)a->half_spinor;
  FourSpinor* clov_oo_psi = (FourSpinor *)a->spinor2;


  // Even even is always cb = 1;

  HalfSpinor *hs0 ALIGN;
  HalfSpinor *hs1 ALIGN;
  HalfSpinor *hs2 ALIGN;
  HalfSpinor *hs3 ALIGN;

  int low = lo+subgrid_vol_cb;
  int high = hi+subgrid_vol_cb;
  
  const int PREFDIST=4;

  /************************ loop over all lattice sites *************************/
  int thissite = s->siteTable(low);
  hs0 =  s->halfspinorBufferOffset(RECONS_GATHER,low,0); 
  _mm_prefetch(hs0, _MM_HINT_NTA);

  hs1 =  s->halfspinorBufferOffset(RECONS_GATHER,low,1); 
  _mm_prefetch(hs1, _MM_HINT_NTA);

  hs2 =  s->halfspinorBufferOffset(RECONS_GATHER,low,2); 
  _mm_prefetch(hs2, _MM_HINT_NTA);

  hs3 =  s->halfspinorBufferOffset(RECONS_GATHER,low,3);
  _mm_prefetch(hs3, _MM_HINT_NTA);

  sn1=&spinor_field[thissite];   
  _mm_prefetch(sn1, _MM_HINT_NTA);

  for (ix1=low+1;ix1<high;ix1++) {
    int nextsite = s->siteTable(ix1);
    sn2 = &spinor_field[nextsite];   
    _mm_prefetch(sn2, _MM_HINT_NTA);
   
    // Reconstruct into temporary...
    recons_4dir_plus(*hs0, *hs1, *hs2, *hs3, *sn1);

    (*sn1)[0][0][0] = clov_oo_psi[thissite][0][0][0] - (*sn1)[0][0][0];
    (*sn1)[0][0][1] = clov_oo_psi[thissite][0][0][1] - (*sn1)[0][0][1];
    (*sn1)[0][1][0] = clov_oo_psi[thissite][0][1][0] - (*sn1)[0][1][0];
    (*sn1)[0][1][1] = clov_oo_psi[thissite][0][1][1] - (*sn1)[0][1][1];
    (*sn1)[0][2][0] = clov_oo_psi[thissite][0][2][0] - (*sn1)[0][2][0];
    (*sn1)[0][2][1] = clov_oo_psi[thissite][0][2][1] - (*sn1)[0][2][1];

    (*sn1)[1][0][0] = clov_oo_psi[thissite][1][0][0] - (*sn1)[1][0][0];
    (*sn1)[1][0][1] = clov_oo_psi[thissite][1][0][1] - (*sn1)[1][0][1];
    (*sn1)[1][1][0] = clov_oo_psi[thissite][1][1][0] - (*sn1)[1][1][0];
    (*sn1)[1][1][1] = clov_oo_psi[thissite][1][1][1] - (*sn1)[1][1][1];
    (*sn1)[1][2][0] = clov_oo_psi[thissite][1][2][0] - (*sn1)[1][2][0];
    (*sn1)[1][2][1] = clov_oo_psi[thissite][1][2][1] - (*sn1)[1][2][1];

    (*sn1)[2][0][0] = clov_oo_psi[thissite][2][0][0] - (*sn1)[2][0][0];
    (*sn1)[2][0][1] = clov_oo_psi[thissite][2][0][1] - (*sn1)[2][0][1];
    (*sn1)[2][1][0] = clov_oo_psi[thissite][2][1][0] - (*sn1)[2][1][0];
    (*sn1)[2][1][1] = clov_oo_psi[thissite][2][1][1] - (*sn1)[2][1][1];
    (*sn1)[2][2][0] = clov_oo_psi[thissite][2][2][0] - (*sn1)[2][2][0];
    (*sn1)[2][2][1] = clov_oo_psi[thissite][2][2][1] - (*sn1)[2][2][1];

    (*sn1)[3][0][0] = clov_oo_psi[thissite][3][0][0] - (*sn1)[3][0][0];
    (*sn1)[3][0][1] = clov_oo_psi[thissite][3][0][1] - (*sn1)[3][0][1];
    (*sn1)[3][1][0] = clov_oo_psi[thissite][3][1][0] - (*sn1)[3][1][0];
    (*sn1)[3][1][1] = clov_oo_psi[thissite][3][1][1] - (*sn1)[3][1][1];
    (*sn1)[3][2][0] = clov_oo_psi[thissite][3][2][0] - (*sn1)[3][2][0];
    (*sn1)[3][2][1] = clov_oo_psi[thissite][3][2][1] - (*sn1)[3][2][1];
 


    hs0 =  s->halfspinorBufferOffset(RECONS_GATHER,ix1,0); 
    _mm_prefetch(hs0, _MM_HINT_NTA);

    hs1 =  s->halfspinorBufferOffset(RECONS_GATHER,ix1,1); 
    _mm_prefetch(hs1, _MM_HINT_NTA);

    hs2 =  s->halfspinorBufferOffset(RECONS_GATHER,ix1,2); 
    _mm_prefetch(hs2, _MM_HINT_NTA);

    hs3 =  s->halfspinorBufferOffset(RECONS_GATHER,ix1,3); 
    _mm_prefetch(hs3, _MM_HINT_NTA);

    sn1=sn2;
    thissite = nextsite;

    /*************************end of loop ****************************/
  }
  recons_4dir_plus(*hs0, *hs1, *hs2, *hs3, (*sn1));
  (*sn1)[0][0][0] = clov_oo_psi[thissite][0][0][0] - (*sn1)[0][0][0];
  (*sn1)[0][0][1] = clov_oo_psi[thissite][0][0][1] - (*sn1)[0][0][1];
  (*sn1)[0][1][0] = clov_oo_psi[thissite][0][1][0] - (*sn1)[0][1][0];
  (*sn1)[0][1][1] = clov_oo_psi[thissite][0][1][1] - (*sn1)[0][1][1];
  (*sn1)[0][2][0] = clov_oo_psi[thissite][0][2][0] - (*sn1)[0][2][0];
  (*sn1)[0][2][1] = clov_oo_psi[thissite][0][2][1] - (*sn1)[0][2][1];
  
  (*sn1)[1][0][0] = clov_oo_psi[thissite][1][0][0] - (*sn1)[1][0][0];
  (*sn1)[1][0][1] = clov_oo_psi[thissite][1][0][1] - (*sn1)[1][0][1];
  (*sn1)[1][1][0] = clov_oo_psi[thissite][1][1][0] - (*sn1)[1][1][0];
  (*sn1)[1][1][1] = clov_oo_psi[thissite][1][1][1] - (*sn1)[1][1][1];
  (*sn1)[1][2][0] = clov_oo_psi[thissite][1][2][0] - (*sn1)[1][2][0];
  (*sn1)[1][2][1] = clov_oo_psi[thissite][1][2][1] - (*sn1)[1][2][1];
  
  (*sn1)[2][0][0] = clov_oo_psi[thissite][2][0][0] - (*sn1)[2][0][0];
  (*sn1)[2][0][1] = clov_oo_psi[thissite][2][0][1] - (*sn1)[2][0][1];
  (*sn1)[2][1][0] = clov_oo_psi[thissite][2][1][0] - (*sn1)[2][1][0];
  (*sn1)[2][1][1] = clov_oo_psi[thissite][2][1][1] - (*sn1)[2][1][1];
  (*sn1)[2][2][0] = clov_oo_psi[thissite][2][2][0] - (*sn1)[2][2][0];
  (*sn1)[2][2][1] = clov_oo_psi[thissite][2][2][1] - (*sn1)[2][2][1];
  
  (*sn1)[3][0][0] = clov_oo_psi[thissite][3][0][0] - (*sn1)[3][0][0];
  (*sn1)[3][0][1] = clov_oo_psi[thissite][3][0][1] - (*sn1)[3][0][1];
  (*sn1)[3][1][0] = clov_oo_psi[thissite][3][1][0] - (*sn1)[3][1][0];
  (*sn1)[3][1][1] = clov_oo_psi[thissite][3][1][1] - (*sn1)[3][1][1];
  (*sn1)[3][2][0] = clov_oo_psi[thissite][3][2][0] - (*sn1)[3][2][0];
  (*sn1)[3][2][1] = clov_oo_psi[thissite][3][2][1] - (*sn1)[3][2][1];
}



void recons_minus_clov_inv(size_t lo,size_t hi, int id, const void *ptr )	
{
  int ix1;
  FourSpinor* sn1 ALIGN;
  FourSpinor* sn2 ALIGN;
  FourSpinor tmp ALIGN;


  const CloverThreadWorkerArgs *a = (CloverThreadWorkerArgs *)ptr;
  ShiftTable<HalfSpinor>* s=(ShiftTable<HalfSpinor>*)a->s;
  int subgrid_vol_cb= s->subgridVolCB();

  FourSpinor* spinor_field = (FourSpinor*)a->spinor;
  FourSpinor* res = (FourSpinor *)a->spinor2;

  CloverTerm* invclov = (CloverTerm *)a->clov;

  HalfSpinor* chi = (HalfSpinor*)a->half_spinor; /* a 1-d std::map of a 2-d array */

  // cb = 0 always here

  HalfSpinor *hs0 ALIGN;
  HalfSpinor *hs1 ALIGN;
  HalfSpinor *hs2 ALIGN;
  HalfSpinor *hs3 ALIGN;

  // Always cb=0
  int low = lo;
  int high = hi;

  int thissite = s->siteTable(low);  
  int clovsite = thissite;

  hs0 =  s->halfspinorBufferOffset(RECONS_GATHER,low,0); 
  _mm_prefetch(hs0, _MM_HINT_NTA);

  hs1 =  s->halfspinorBufferOffset(RECONS_GATHER,low,1); 
  _mm_prefetch(hs1, _MM_HINT_NTA);

  hs2 =  s->halfspinorBufferOffset(RECONS_GATHER,low,2); 
  _mm_prefetch(hs2, _MM_HINT_NTA);

  hs3 =  s->halfspinorBufferOffset(RECONS_GATHER,low,3);
  _mm_prefetch(hs3, _MM_HINT_NTA);

  sn1=&spinor_field[thissite];   
  _mm_prefetch(sn1, _MM_HINT_NTA);

  for (ix1=low+1;ix1<high;ix1++) {
    thissite = s->siteTable(ix1);
     sn2 = &spinor_field[thissite];   
    _mm_prefetch(sn2, _MM_HINT_NTA);

    recons_4dir_minus(*hs0, *hs1, *hs2, *hs3, *sn1);
    cloverSiteApply(res[clovsite], invclov[clovsite], *sn1);

    hs0 =  s->halfspinorBufferOffset(RECONS_GATHER,ix1,0); 
    _mm_prefetch(hs0, _MM_HINT_NTA);
    
    hs1 =  s->halfspinorBufferOffset(RECONS_GATHER,ix1,1); 
    _mm_prefetch(hs1, _MM_HINT_NTA);
    
    hs2 =  s->halfspinorBufferOffset(RECONS_GATHER,ix1,2); 
    _mm_prefetch(hs2, _MM_HINT_NTA);
    
    hs3 =  s->halfspinorBufferOffset(RECONS_GATHER,ix1,3); 
    _mm_prefetch(hs3, _MM_HINT_NTA);
   
    sn1=sn2;
    clovsite = thissite;
  }
  recons_4dir_minus(*hs0, *hs1, *hs2, *hs3, *sn1);
  cloverSiteApply(res[clovsite], invclov[clovsite], *sn1);

}

void recons_minus_final(size_t lo,size_t hi, int id, const void *ptr )	
{
  int ix1;
  FourSpinor* sn1 ALIGN;
  FourSpinor* sn2 ALIGN;
  FourSpinor tmp ALIGN;

  const CloverThreadWorkerArgs *a = (CloverThreadWorkerArgs *)ptr;
  ShiftTable<HalfSpinor>* s=(ShiftTable<HalfSpinor>*)a->s;
  int subgrid_vol_cb= s->subgridVolCB();

  FourSpinor* spinor_field = (FourSpinor*)a->spinor;
  FourSpinor* clov_oo_psi = (FourSpinor*)a->spinor2;
  HalfSpinor* chi = (HalfSpinor*)a->half_spinor; /* a 1-d std::map of a 2-d array */

  // cb = 0 always here
  HalfSpinor *hs0 ALIGN;
  HalfSpinor *hs1 ALIGN;
  HalfSpinor *hs2 ALIGN;
  HalfSpinor *hs3 ALIGN;

  // Always cb=1
  int low = lo+subgrid_vol_cb;
  int high = hi+subgrid_vol_cb;

  int thissite = s->siteTable(low);  


  hs0 =  s->halfspinorBufferOffset(RECONS_GATHER,low,0); 
  _mm_prefetch(hs0, _MM_HINT_NTA);

  hs1 =  s->halfspinorBufferOffset(RECONS_GATHER,low,1); 
  _mm_prefetch(hs1, _MM_HINT_NTA);

  hs2 =  s->halfspinorBufferOffset(RECONS_GATHER,low,2); 
  _mm_prefetch(hs2, _MM_HINT_NTA);

  hs3 =  s->halfspinorBufferOffset(RECONS_GATHER,low,3);
  _mm_prefetch(hs3, _MM_HINT_NTA);

  sn1=&spinor_field[thissite];   
  _mm_prefetch(sn1, _MM_HINT_NTA);

  for (ix1=low+1;ix1<high;ix1++) {
    int nextsite = s->siteTable(ix1);
     sn2 = &spinor_field[nextsite];   
    _mm_prefetch(sn2, _MM_HINT_NTA);

    recons_4dir_minus(*hs0, *hs1, *hs2, *hs3, (*sn1));
    (*sn1)[0][0][0] = clov_oo_psi[thissite][0][0][0] - (*sn1)[0][0][0];
    (*sn1)[0][0][1] = clov_oo_psi[thissite][0][0][1] - (*sn1)[0][0][1];
    (*sn1)[0][1][0] = clov_oo_psi[thissite][0][1][0] - (*sn1)[0][1][0];
    (*sn1)[0][1][1] = clov_oo_psi[thissite][0][1][1] - (*sn1)[0][1][1];
    (*sn1)[0][2][0] = clov_oo_psi[thissite][0][2][0] - (*sn1)[0][2][0];
    (*sn1)[0][2][1] = clov_oo_psi[thissite][0][2][1] - (*sn1)[0][2][1];
  
    (*sn1)[1][0][0] = clov_oo_psi[thissite][1][0][0] - (*sn1)[1][0][0];
    (*sn1)[1][0][1] = clov_oo_psi[thissite][1][0][1] - (*sn1)[1][0][1];
    (*sn1)[1][1][0] = clov_oo_psi[thissite][1][1][0] - (*sn1)[1][1][0];
    (*sn1)[1][1][1] = clov_oo_psi[thissite][1][1][1] - (*sn1)[1][1][1];
    (*sn1)[1][2][0] = clov_oo_psi[thissite][1][2][0] - (*sn1)[1][2][0];
    (*sn1)[1][2][1] = clov_oo_psi[thissite][1][2][1] - (*sn1)[1][2][1];
    
    (*sn1)[2][0][0] = clov_oo_psi[thissite][2][0][0] - (*sn1)[2][0][0];
    (*sn1)[2][0][1] = clov_oo_psi[thissite][2][0][1] - (*sn1)[2][0][1];
    (*sn1)[2][1][0] = clov_oo_psi[thissite][2][1][0] - (*sn1)[2][1][0];
    (*sn1)[2][1][1] = clov_oo_psi[thissite][2][1][1] - (*sn1)[2][1][1];
    (*sn1)[2][2][0] = clov_oo_psi[thissite][2][2][0] - (*sn1)[2][2][0];
    (*sn1)[2][2][1] = clov_oo_psi[thissite][2][2][1] - (*sn1)[2][2][1];
    
    (*sn1)[3][0][0] = clov_oo_psi[thissite][3][0][0] - (*sn1)[3][0][0];
    (*sn1)[3][0][1] = clov_oo_psi[thissite][3][0][1] - (*sn1)[3][0][1];
    (*sn1)[3][1][0] = clov_oo_psi[thissite][3][1][0] - (*sn1)[3][1][0];
    (*sn1)[3][1][1] = clov_oo_psi[thissite][3][1][1] - (*sn1)[3][1][1];
    (*sn1)[3][2][0] = clov_oo_psi[thissite][3][2][0] - (*sn1)[3][2][0];
    (*sn1)[3][2][1] = clov_oo_psi[thissite][3][2][1] - (*sn1)[3][2][1];
    
    
    hs0 =  s->halfspinorBufferOffset(RECONS_GATHER,ix1,0); 
    _mm_prefetch(hs0, _MM_HINT_NTA);
    
    hs1 =  s->halfspinorBufferOffset(RECONS_GATHER,ix1,1); 
    _mm_prefetch(hs1, _MM_HINT_NTA);
    
    hs2 =  s->halfspinorBufferOffset(RECONS_GATHER,ix1,2); 
    _mm_prefetch(hs2, _MM_HINT_NTA);
    
    hs3 =  s->halfspinorBufferOffset(RECONS_GATHER,ix1,3); 
    _mm_prefetch(hs3, _MM_HINT_NTA);
   
    sn1=sn2;
    thissite = nextsite;
  }
  recons_4dir_minus(*hs0, *hs1, *hs2, *hs3, (*sn1));

  (*sn1)[0][0][0] = clov_oo_psi[thissite][0][0][0] - (*sn1)[0][0][0];
  (*sn1)[0][0][1] = clov_oo_psi[thissite][0][0][1] - (*sn1)[0][0][1];
  (*sn1)[0][1][0] = clov_oo_psi[thissite][0][1][0] - (*sn1)[0][1][0];
  (*sn1)[0][1][1] = clov_oo_psi[thissite][0][1][1] - (*sn1)[0][1][1];
  (*sn1)[0][2][0] = clov_oo_psi[thissite][0][2][0] - (*sn1)[0][2][0];
  (*sn1)[0][2][1] = clov_oo_psi[thissite][0][2][1] - (*sn1)[0][2][1];
  
  (*sn1)[1][0][0] = clov_oo_psi[thissite][1][0][0] - (*sn1)[1][0][0];
  (*sn1)[1][0][1] = clov_oo_psi[thissite][1][0][1] - (*sn1)[1][0][1];
  (*sn1)[1][1][0] = clov_oo_psi[thissite][1][1][0] - (*sn1)[1][1][0];
  (*sn1)[1][1][1] = clov_oo_psi[thissite][1][1][1] - (*sn1)[1][1][1];
  (*sn1)[1][2][0] = clov_oo_psi[thissite][1][2][0] - (*sn1)[1][2][0];
  (*sn1)[1][2][1] = clov_oo_psi[thissite][1][2][1] - (*sn1)[1][2][1];
  
  (*sn1)[2][0][0] = clov_oo_psi[thissite][2][0][0] - (*sn1)[2][0][0];
  (*sn1)[2][0][1] = clov_oo_psi[thissite][2][0][1] - (*sn1)[2][0][1];
  (*sn1)[2][1][0] = clov_oo_psi[thissite][2][1][0] - (*sn1)[2][1][0];
  (*sn1)[2][1][1] = clov_oo_psi[thissite][2][1][1] - (*sn1)[2][1][1];
  (*sn1)[2][2][0] = clov_oo_psi[thissite][2][2][0] - (*sn1)[2][2][0];
  (*sn1)[2][2][1] = clov_oo_psi[thissite][2][2][1] - (*sn1)[2][2][1];
  
  (*sn1)[3][0][0] = clov_oo_psi[thissite][3][0][0] - (*sn1)[3][0][0];
  (*sn1)[3][0][1] = clov_oo_psi[thissite][3][0][1] - (*sn1)[3][0][1];
  (*sn1)[3][1][0] = clov_oo_psi[thissite][3][1][0] - (*sn1)[3][1][0];
  (*sn1)[3][1][1] = clov_oo_psi[thissite][3][1][1] - (*sn1)[3][1][1];
  (*sn1)[3][2][0] = clov_oo_psi[thissite][3][2][0] - (*sn1)[3][2][0];
  (*sn1)[3][2][1] = clov_oo_psi[thissite][3][2][1] - (*sn1)[3][2][1];
  
}









/*****************end of isign corresponding to -1 **********************/

  }



// Your actual operator
void CloverSchur4D<float>::operator()(float* res, 
				      const float* psi, 
				      const float* u, 
				      const float* clov_oo,
				      const float* invclov_ee,
				      int isign)
{
  HalfSpinor* chi1 = tab->getChi1();
  HalfSpinor* chi2 = tab->getChi2();
  int subgrid_vol_cb = s_tab->subgridVolCB();


  if(isign==1) {

    // Apply Dslash (source cb = 1)
    tab->startReceives();

    CPlusPlusWilsonDslash::dispatchToThreads(DslashParscalar32Bit::decomp_plus,
		      (void *)psi,
		      (void *)chi1,
		      (void *)u,
		      (void *)s_tab,
		      1,
		      subgrid_vol_cb);

    tab->startSendForward(); 

    CPlusPlusWilsonDslash::dispatchToThreads(DslashParscalar32Bit::decomp_hvv_plus,
		      (void*)psi,
		      (void*)chi2,
		      (void*)u,
		      (void*)s_tab,
		      1,
		      (int)subgrid_vol_cb);

    tab->startSendBack();

    // All comms started do the clover odd odd bit...
    CPlusPlusClover::dispatchToThreads(CPlusPlusClover::CPlusPlusClover32Bit::clovOddOddApply,
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


    CPlusPlusWilsonDslash::dispatchToThreads(DslashParscalar32Bit::mvv_recons_plus,
					     (void*)t_spinor2,  // Even part
					     (void *)chi1,
					     (void *)u,
					     (void *)s_tab,
					     0,
					     subgrid_vol_cb);

    tab->finishSendBack();
    tab->finishReceiveFromForward();    



    CPlusPlusClover::dispatchToThreads(CloverParscalar32Bit::recons_plus_clov_inv,
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

    CPlusPlusWilsonDslash::dispatchToThreads(DslashParscalar32Bit::decomp_plus,
					     (void *)t_spinor,         // Decompose on even part
					     (void *)chi1,
					     (void *)u,
					     (void *)s_tab,
					     0,
					     subgrid_vol_cb);

    tab->startSendForward(); 

    CPlusPlusWilsonDslash::dispatchToThreads(DslashParscalar32Bit::decomp_hvv_plus,
					     (void*)t_spinor,         // Decompose on even part
					     (void*)chi2,
					     (void*)u,
					     (void*)s_tab,
					     0,
					     (int)subgrid_vol_cb);

    tab->startSendBack();


    // All comms started 
    // All comms started do the clover odd odd bit...
    CPlusPlusClover::dispatchToThreads(CPlusPlusClover::CPlusPlusClover32Bit::clovOddOddApply,
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


    CPlusPlusWilsonDslash::dispatchToThreads(DslashParscalar32Bit::mvv_recons_plus,
					     (void*)res,
					     (void *)chi1,
					     (void *)u,
					     (void *)s_tab,
					     1,
					     subgrid_vol_cb);

    tab->finishSendBack();
    tab->finishReceiveFromForward();    

    CPlusPlusClover::dispatchToThreads(CloverParscalar32Bit::recons_plus_final,
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

    CPlusPlusWilsonDslash::dispatchToThreads(DslashParscalar32Bit::decomp_minus,
		      (void *)psi,
		      (void *)chi1,
		      (void *)u,
		      (void *)s_tab,
		      1,
		      subgrid_vol_cb);

    tab->startSendForward(); 

    CPlusPlusWilsonDslash::dispatchToThreads(DslashParscalar32Bit::decomp_hvv_minus,
		      (void*)psi,
		      (void*)chi2,
		      (void*)u,
		      (void*)s_tab,
		      1,
		      (int)subgrid_vol_cb);

    tab->startSendBack();


    // All comms started do the clover odd odd bit...
    CPlusPlusClover::dispatchToThreads(CPlusPlusClover::CPlusPlusClover32Bit::clovOddOddApply,
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


    CPlusPlusWilsonDslash::dispatchToThreads(DslashParscalar32Bit::mvv_recons_minus,
					     (void*)t_spinor2, // even half
					     (void *)chi1,
					     (void *)u,
					     (void *)s_tab,
					     0,
					     subgrid_vol_cb);

    tab->finishSendBack();
    tab->finishReceiveFromForward();    

    CPlusPlusClover::dispatchToThreads(CloverParscalar32Bit::recons_minus_clov_inv,
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

    CPlusPlusWilsonDslash::dispatchToThreads(DslashParscalar32Bit::decomp_minus,
		      (void *)t_spinor,
		      (void *)chi1,
		      (void *)u,
		      (void *)s_tab,
		      0,
		      subgrid_vol_cb);

    tab->startSendForward(); 

    CPlusPlusWilsonDslash::dispatchToThreads(DslashParscalar32Bit::decomp_hvv_minus,
		      (void*)t_spinor,
		      (void*)chi2,
		      (void*)u,
		      (void*)s_tab,
		      0,
		      (int)subgrid_vol_cb);

    tab->startSendBack();


    // All comms started 
    // All comms started do the clover odd odd bit...
    CPlusPlusClover::dispatchToThreads(CPlusPlusClover::CPlusPlusClover32Bit::clovOddOddApply,
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


    CPlusPlusWilsonDslash::dispatchToThreads(DslashParscalar32Bit::mvv_recons_minus,
					     (void*)res,
					     (void *)chi1,
					     (void *)u,
					     (void *)s_tab,
					     1,
					     subgrid_vol_cb);

    tab->finishSendBack();
    tab->finishReceiveFromForward();    

    CPlusPlusClover::dispatchToThreads(CloverParscalar32Bit::recons_minus_final,
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
CloverSchur4D<float>::CloverSchur4D(const int latt_size[],      
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
  CloverSchur4D<float>::~CloverSchur4D() 
  {
    delete s_tab;
    delete tab;
    CPlusPlusWilsonDslash::dealloc(xt_spinor);
    CPlusPlusWilsonDslash::dealloc(xt_spinor2);

  }


} // Namespace
