#ifndef SHIFT_TABLE_TYPES_PARSCALAR_H
#define SHIFT_TABLE_TYPES_PARSCALAR_H

namespace CPlusPlusWilsonDslash {

    enum HalfSpinorOffsetType {
      DECOMP_SCATTER=0,
      DECOMP_HVV_SCATTER,
      RECONS_MVV_GATHER,
      RECONS_GATHER
    } ;

    struct InvTab { 
      int cb;
      int linearcb;
    } ;

}

#endif
