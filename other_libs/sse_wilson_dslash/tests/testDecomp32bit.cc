#include "qdp.h"
#include "testDecomp.h"
#include "unittest.h"

using namespace QDP;
using namespace Assertions;
#include <sse32.h>
#include <sse_align.h>
#include <types32.h>
#include <cmath>
#include <decomp.h>

/* Spin Matrices */
/* gamma 0 */
#define _sse_42_gamma0_minus()     _sse_vector_xch_i_sub()
#define _sse_42_gamma0_plus()      _sse_vector_xch_i_add()


/* gamma 1 */
#define _sse_42_gamma1_minus()     _sse_vector_xch();\
				   _sse_vector_addsub()

#define _sse_42_gamma1_plus()      _sse_vector_xch();\
				   _sse_vector_subadd()


/* gamma 2 */
#define _sse_42_gamma2_minus()     _sse_vector_i_subadd()
#define _sse_42_gamma2_plus()      _sse_vector_i_addsub()

/* gamma 3 */
#define _sse_42_gamma3_minus()     _sse_vector_sub()
#define _sse_42_gamma3_plus()      _sse_vector_add()

void
testDecomp0Minus::run()
{
  spinor_array spinor;
  halfspinor_array hspinor1, hspinor2;

  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	hspinor1[col][spin2][reim] = 0;
	hspinor2[col][spin2][reim] = 0;
      }
    }

    for(int spin4=0; spin4 < 4; spin4++) {
      for(int reim=0; reim < 2; reim++) { 

	spinor[spin4][col][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;

      }
    }
  } // Color

  spinor_array* sp1 = &spinor;
  halfspinor_array *s3 = &hspinor1;

  // Copy from code 
  _sse_pair_load((*sp1)[0],(*sp1)[1]);
      
  //  s3 = chi+ halfspinor_buffer_offset(DECOMP_SCATTER,ix1,0);
  
  _sse_pair_load_up((*sp1)[2],(*sp1)[3]);
  
  _sse_42_gamma0_minus();
  
  _sse_vector_store(*s3);


  // My version
  decomp_gamma0_minus(spinor, hspinor2);

  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	float diff = hspinor1[col][spin2][reim] - hspinor2[col][spin2][reim];
	diff /= (float)(3*2*2); // per number
	//	QDPIO::cout << "   col="<<col<<" s="<<spin2<<" reim=" <<reim<< "  diff = " << diff << std::endl;
	assertion( fabs(diff) < 1.0e-9 );
      }
    }
  }


}

void
testDecomp1Minus::run()
{
  spinor_array spinor;
  halfspinor_array hspinor1, hspinor2;

  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	hspinor1[col][spin2][reim] = 0;
	hspinor2[col][spin2][reim] = 0;
      }
    }

    for(int spin4=0; spin4 < 4; spin4++) {
      for(int reim=0; reim < 2; reim++) { 

	spinor[spin4][col][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;

      }
    }
  } // Color

  spinor_array* sp1 = &spinor;
  halfspinor_array *s3 = &hspinor1;

  // Copy from code 
  _sse_pair_load((*sp1)[0],(*sp1)[1]);
  _sse_pair_load_up((*sp1)[2],(*sp1)[3]);
  
  _sse_42_gamma1_minus();
  
  _sse_vector_store(*s3);


  // My version
  decomp_gamma1_minus(spinor, hspinor2);

  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	float diff = hspinor1[col][spin2][reim] - hspinor2[col][spin2][reim];
	diff /= (float)(3*2*2); // per number
	//	QDPIO::cout << "   col="<<col<<" s="<<spin2<<" reim=" <<reim<< "  diff = " << diff << std::endl;
	assertion( fabs(diff) < 1.0e-9 );
      }
    }
  }


}

void
testDecomp2Minus::run()
{
  spinor_array spinor;
  halfspinor_array hspinor1, hspinor2;

  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	hspinor1[col][spin2][reim] = 0;
	hspinor2[col][spin2][reim] = 0;
      }
    }

    for(int spin4=0; spin4 < 4; spin4++) {
      for(int reim=0; reim < 2; reim++) { 

	spinor[spin4][col][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;

      }
    }
  } // Color

  spinor_array* sp1 = &spinor;
  halfspinor_array *s3 = &hspinor1;

  // Copy from code 
  _sse_pair_load((*sp1)[0],(*sp1)[1]);
      
  //  s3 = chi+ halfspinor_buffer_offset(DECOMP_SCATTER,ix1,0);
  
  _sse_pair_load_up((*sp1)[2],(*sp1)[3]);
  
  _sse_42_gamma2_minus();
  
  _sse_vector_store(*s3);


  // My version
  decomp_gamma2_minus(spinor, hspinor2);

  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	float diff = hspinor1[col][spin2][reim] - hspinor2[col][spin2][reim];
	diff /= (float)(3*2*2); // per number
	//	QDPIO::cout << "   col="<<col<<" s="<<spin2<<" reim=" <<reim<< "  diff = " << diff << std::endl;
	assertion( fabs(diff) < 1.0e-9 );
      }
    }
  }


}

void
testDecomp3Minus::run()
{
  spinor_array spinor;
  halfspinor_array hspinor1, hspinor2;

  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	hspinor1[col][spin2][reim] = 0;
	hspinor2[col][spin2][reim] = 0;
      }
    }

    for(int spin4=0; spin4 < 4; spin4++) {
      for(int reim=0; reim < 2; reim++) { 

	spinor[spin4][col][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;

      }
    }
  } // Color

  spinor_array* sp1 = &spinor;
  halfspinor_array *s3 = &hspinor1;

  // Copy from code 
  _sse_pair_load((*sp1)[0],(*sp1)[1]);
      
  //  s3 = chi+ halfspinor_buffer_offset(DECOMP_SCATTER,ix1,0);
  
  _sse_pair_load_up((*sp1)[2],(*sp1)[3]);
  
  _sse_42_gamma3_minus();
  
  _sse_vector_store(*s3);


  // My version
  decomp_gamma3_minus(spinor, hspinor2);

  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	float diff = hspinor1[col][spin2][reim] - hspinor2[col][spin2][reim];
	diff /= (float)(3*2*2); // per number
	//	QDPIO::cout << "   col="<<col<<" s="<<spin2<<" reim=" <<reim<< "  diff = " << diff << std::endl;
	assertion( fabs(diff) < 1.0e-9 );
      }
    }
  }


}

void
testDecomp0Plus::run()
{
  spinor_array spinor;
  halfspinor_array hspinor1, hspinor2;

  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	hspinor1[col][spin2][reim] = 0;
	hspinor2[col][spin2][reim] = 0;
      }
    }

    for(int spin4=0; spin4 < 4; spin4++) {
      for(int reim=0; reim < 2; reim++) { 

	spinor[spin4][col][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;

      }
    }
  } // Color

  spinor_array* sp1 = &spinor;
  halfspinor_array *s3 = &hspinor1;

  // Copy from code 
  _sse_pair_load((*sp1)[0],(*sp1)[1]);
  _sse_pair_load_up((*sp1)[2],(*sp1)[3]);
  _sse_42_gamma0_plus();
  _sse_vector_store(*s3);


  // My version
  decomp_gamma0_plus(spinor, hspinor2);

  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	float diff = hspinor1[col][spin2][reim] - hspinor2[col][spin2][reim];
	diff /= (float)(3*2*2); // per number
	//	QDPIO::cout << "   col="<<col<<" s="<<spin2<<" reim=" <<reim<< "  diff = " << diff << std::endl;
	assertion( fabs(diff) < 1.0e-9 );
      }
    }
  }


}

void
testDecomp1Plus::run()
{
  spinor_array spinor;
  halfspinor_array hspinor1, hspinor2;

  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	hspinor1[col][spin2][reim] = 0;
	hspinor2[col][spin2][reim] = 0;
      }
    }

    for(int spin4=0; spin4 < 4; spin4++) {
      for(int reim=0; reim < 2; reim++) { 

	spinor[spin4][col][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;

      }
    }
  } // Color

  spinor_array* sp1 = &spinor;
  halfspinor_array *s3 = &hspinor1;

  // Copy from code 
  _sse_pair_load((*sp1)[0],(*sp1)[1]);
  _sse_pair_load_up((*sp1)[2],(*sp1)[3]);
  _sse_42_gamma1_plus();
  _sse_vector_store(*s3);


  // My version
  decomp_gamma1_plus(spinor, hspinor2);

  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	float diff = hspinor1[col][spin2][reim] - hspinor2[col][spin2][reim];
	diff /= (float)(3*2*2); // per number
	//	QDPIO::cout << "   col="<<col<<" s="<<spin2<<" reim=" <<reim<< "  diff = " << diff << std::endl;
	assertion( fabs(diff) < 1.0e-9 );
      }
    }
  }


}

void
testDecomp2Plus::run()
{
  spinor_array spinor;
  halfspinor_array hspinor1, hspinor2;

  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	hspinor1[col][spin2][reim] = 0;
	hspinor2[col][spin2][reim] = 0;
      }
    }

    for(int spin4=0; spin4 < 4; spin4++) {
      for(int reim=0; reim < 2; reim++) { 

	spinor[spin4][col][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;

      }
    }
  } // Color

  spinor_array* sp1 = &spinor;
  halfspinor_array *s3 = &hspinor1;

  // Copy from code 
  _sse_pair_load((*sp1)[0],(*sp1)[1]);
  _sse_pair_load_up((*sp1)[2],(*sp1)[3]);
  _sse_42_gamma2_plus();
  _sse_vector_store(*s3);


  // My version
  decomp_gamma2_plus(spinor, hspinor2);

  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	float diff = hspinor1[col][spin2][reim] - hspinor2[col][spin2][reim];
	diff /= (float)(3*2*2); // per number
	//	QDPIO::cout << "   col="<<col<<" s="<<spin2<<" reim=" <<reim<< "  diff = " << diff << std::endl;
	assertion( fabs(diff) < 1.0e-9 );
      }
    }
  }


}

void
testDecomp3Plus::run()
{
  spinor_array spinor;
  halfspinor_array hspinor1, hspinor2;

  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	hspinor1[col][spin2][reim] = 0;
	hspinor2[col][spin2][reim] = 0;
      }
    }

    for(int spin4=0; spin4 < 4; spin4++) {
      for(int reim=0; reim < 2; reim++) { 

	spinor[spin4][col][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;

      }
    }
  } // Color

  spinor_array* sp1 = &spinor;
  halfspinor_array *s3 = &hspinor1;

  // Copy from code 
  _sse_pair_load((*sp1)[0],(*sp1)[1]);
  _sse_pair_load_up((*sp1)[2],(*sp1)[3]);
  _sse_42_gamma3_plus();
  _sse_vector_store(*s3);


  // My version
  decomp_gamma3_plus(spinor, hspinor2);

  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	float diff = hspinor1[col][spin2][reim] - hspinor2[col][spin2][reim];
	diff /= (float)(3*2*2); // per number
	//	QDPIO::cout << "   col="<<col<<" s="<<spin2<<" reim=" <<reim<< "  diff = " << diff << std::endl;
	assertion( fabs(diff) < 1.0e-9 );
      }
    }
  }


}
