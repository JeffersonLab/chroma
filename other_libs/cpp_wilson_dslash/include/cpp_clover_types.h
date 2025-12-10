#ifndef CPP_CLOVER_TYPES_H
#define CPP_CLOVER_TYPES_H

namespace CPlusPlusClover {

  namespace Clover32BitTypes { 

    struct CloverStruct { 
      float diag[8]; // Use only 6
      float off_diag[16][2]; // Use only 15
    };
    typedef CloverStruct CloverTerm[2];

    // typedef float ClovFourSpinor[12][2];
  }


  namespace Clover64BitTypes {

    struct CloverStruct { 
      double diag[8]; // Use only 6
      double off_diag[16][2]; // Use only 15
    };

    typedef CloverStruct CloverTerm[2];
    // typedef double ClovFourSpinor[12][2];
  }



}


#endif
