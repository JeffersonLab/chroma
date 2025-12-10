#ifndef CPP_CLOVER_ODD_ODD_PARSCALAR_32BIT_H
#define CPP_CLOVER_ODD_ODD_PARSCALAR_32BIT_H

#include "cpp_dslash_types.h"
#include "cpp_clover_types.h"
#include "shift_table_parscalar.h"
#include "dispatch_parscalar.h"
#include "cpp_clover_site_apply_32bit.h"

namespace CPlusPlusClover { 
  
  
  namespace CPlusPlusClover32Bit { 
    using namespace Dslash32BitTypes;
    using namespace Clover32BitTypes;

    inline void clovOddOddApply(size_t lo, size_t hi, int id, const void *args)
    {
      const CloverThreadWorkerArgs *a=(const CloverThreadWorkerArgs *)args;
      const CloverTerm *clover_term = (const CloverTerm *)a->clov;
      FourSpinor *dst_spinor = (FourSpinor *)a->spinor2;
      const FourSpinor *src_spinor = (FourSpinor *)a->spinor;
      ShiftTable<HalfSpinor>* s_tab = (ShiftTable<HalfSpinor> *)a->s;
      const int subgrid_vol_cb = s_tab->subgridVolCB();
      int subgrid_vol_cb_by_two = subgrid_vol_cb/2;
      int half = a->half;
      
      // Always odd subset so low = lo+subgrid_vol_cb
      //                      high = hi+subgrid_vol_cb
      
      // We can also add on half a subgrid volume depending on
      // which half we are doing
      // half = 0 => first half
      // half = 1 => second half
      //
      
      int low = lo + subgrid_vol_cb+subgrid_vol_cb_by_two*half;
      int high = hi + subgrid_vol_cb+subgrid_vol_cb_by_two*half;
      
      for(int i = low; i < high; i++) {
	int thissite = s_tab->siteTable(i);
	cloverSiteApply(dst_spinor[thissite], clover_term[thissite], src_spinor[thissite]);
      }
    }
  }
}

#endif
