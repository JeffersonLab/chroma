#include "qdp.h"
#include "testDecomp.h"
#include "unittest.h"

using namespace QDP;
using namespace Assertions;
#include <sse64.h>
#include <sse_align.h>
#include <types64.h>
#include <cmath>
#include <decomp.h>


#define _sse_42_1_gamma0_minus(sp) \
      _sse_load((sp)[0]); \
      _sse_load_up((sp)[3]);\
      _sse_vector_i_mul();\
      _sse_vector_sub()
	  
#define _sse_42_2_gamma0_minus(sp) \
      _sse_load((sp)[1]);\
      _sse_load_up((sp)[2]);\
      _sse_vector_i_mul();\
      _sse_vector_sub()


#define _sse_42_1_gamma0_plus(sm) \
      _sse_load((sm)[0]);\
      _sse_load_up((sm)[3]);\
      _sse_vector_i_mul();\
      _sse_vector_add()

#define _sse_42_2_gamma0_plus(sm) \
	  _sse_load((sm)[1]);\
      _sse_load_up((sm)[2]);\
      _sse_vector_i_mul();\
      _sse_vector_add()

/* gamma 1 */


#define _sse_42_1_gamma1_minus(sp) \
      _sse_load((sp)[0]);\
      _sse_load_up((sp)[3]);\
      _sse_vector_add()
	  
#define _sse_42_2_gamma1_minus(sp) \
      _sse_load((sp)[1]);\
      _sse_load_up((sp)[2]);\
      _sse_vector_sub()


#define _sse_42_1_gamma1_plus(sm) \
      _sse_load((sm)[0]);\
      _sse_load_up((sm)[3]);\
      _sse_vector_sub()


#define _sse_42_2_gamma1_plus(sm) \
	  _sse_load((sm)[1]);\
      _sse_load_up((sm)[2]);\
      _sse_vector_add()


/* gamma 2 */


#define _sse_42_1_gamma2_minus(sp) \
      _sse_load((sp)[0]);\
      _sse_load_up((sp)[2]);\
      _sse_vector_i_mul();\
      _sse_vector_sub()
	  
#define _sse_42_2_gamma2_minus(sp) \
      _sse_load((sp)[1]);\
      _sse_load_up((sp)[3]);\
      _sse_vector_i_mul();\
      _sse_vector_add()

#define _sse_42_1_gamma2_plus(sm) \
      _sse_load((sm)[0]);\
      _sse_load_up((sm)[2]);\
      _sse_vector_i_mul();\
      _sse_vector_add()
 
#define _sse_42_2_gamma2_plus(sm) \
	  _sse_load((sm)[1]);\
      _sse_load_up((sm)[3]);\
      _sse_vector_i_mul();\
      _sse_vector_sub()


/* gamma 3 */
#define _sse_42_1_gamma3_minus(sp) \
	  _sse_load((sp)[0]); \
	  _sse_load_up((sp)[2]); \
      _sse_vector_sub()

	  
#define _sse_42_2_gamma3_minus(sp) \
      _sse_load((sp)[1]);\
      _sse_load_up((sp)[3]);\
      _sse_vector_sub()


#define _sse_42_1_gamma3_plus(sm) \
      _sse_load((sm)[0]);\
      _sse_load_up((sm)[2]);\
      _sse_vector_add()


#define _sse_42_2_gamma3_plus(sm) \
	  _sse_load((sm)[1]);\
      _sse_load_up((sm)[3]);\
      _sse_vector_add()



void
testDecomp0Minus::run()
{
  spinor_array spinor;
  halfspinor_array hspinor1, hspinor2;

  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	hspinor1[spin2][col][reim] = 0;
	hspinor2[spin2][col][reim] = 0;
      }
    }

    for(int spin4=0; spin4 < 4; spin4++) {
      for(int reim=0; reim < 2; reim++) { 

	spinor[spin4][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;

      }
    }
  } // Color

  spinor_array* sp = &spinor;
  halfspinor_array *s3 = &hspinor1;

 _sse_42_1_gamma0_minus(*sp);
 _sse_store((*s3)[0]);
 /*spin decomposition of first component of halfspinor */
 _sse_42_2_gamma0_minus(*sp);
 _sse_store((*s3)[1]);

  // My version
  decomp_gamma0_minus(spinor, hspinor2);


  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	double diff = hspinor1[spin2][col][reim] - hspinor2[spin2][col][reim];
	diff /= (double)(3*2*2);
	assertion( fabs(diff)< 1.0e-17 );
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
	hspinor1[spin2][col][reim] = 0;
	hspinor2[spin2][col][reim] = 0;
      }
    }

    for(int spin4=0; spin4 < 4; spin4++) {
      for(int reim=0; reim < 2; reim++) { 

	spinor[spin4][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;

      }
    }
  } // Color

  spinor_array* sp = &spinor;
  halfspinor_array *s3 = &hspinor1;

  _sse_42_1_gamma1_minus(*sp);
  _sse_store((*s3)[0]);
  
  _sse_42_2_gamma1_minus(*sp);
  _sse_store((*s3)[1]);

  // My version
  decomp_gamma1_minus(spinor, hspinor2);

  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	double diff = hspinor1[spin2][col][reim] - hspinor2[spin2][col][reim];
	diff /= (double)(3*2*2);
	//	QDPIO::cout << "   col="<<col<<" s="<<spin2<<" reim=" <<reim<< "  diff = " << diff << std::endl;
	assertion( fabs(diff)< 1.0e-17 );
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
	hspinor1[spin2][col][reim] = 0;
	hspinor2[spin2][col][reim] = 0;
      }
    }

    for(int spin4=0; spin4 < 4; spin4++) {
      for(int reim=0; reim < 2; reim++) { 

	spinor[spin4][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;

      }
    }
  } // Color

  spinor_array* sp = &spinor;
  halfspinor_array *s3 = &hspinor1;

  _sse_42_1_gamma2_minus(*sp);
  _sse_store((*s3)[0]);
  
  _sse_42_2_gamma2_minus(*sp);
  _sse_store((*s3)[1]);

  // My version
  decomp_gamma2_minus(spinor, hspinor2);

  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	double diff = hspinor1[spin2][col][reim] - hspinor2[spin2][col][reim];
	diff /= (double)(3*2*2);
	//QDPIO::cout << "   col="<<col<<" s="<<spin2<<" reim=" <<reim<< "  diff = " << diff << std::endl;
	assertion( fabs(diff)< 1.0e-17 );
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
	hspinor1[spin2][col][reim] = 0;
	hspinor2[spin2][col][reim] = 0;
      }
    }

    for(int spin4=0; spin4 < 4; spin4++) {
      for(int reim=0; reim < 2; reim++) { 

	spinor[spin4][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;

      }
    }
  } // Color

  spinor_array* sp = &spinor;
  halfspinor_array *s3 = &hspinor1;

  // Copy from code 
    _sse_42_1_gamma3_minus(*sp);
    _sse_store((*s3)[0]);

    _sse_42_2_gamma3_minus(*sp);
    _sse_store((*s3)[1]);

  // My version
  decomp_gamma3_minus(spinor, hspinor2);

  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	double diff = hspinor1[spin2][col][reim] - hspinor2[spin2][col][reim];
	diff /= (double)(3*2*2);
	//	QDPIO::cout << "   col="<<col<<" s="<<spin2<<" reim=" <<reim<< "  diff = " << diff << std::endl;
	assertion( fabs(diff)< 1.0e-17 );
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
	hspinor1[spin2][col][reim] = 0;
	hspinor2[spin2][col][reim] = 0;
      }
    }

    for(int spin4=0; spin4 < 4; spin4++) {
      for(int reim=0; reim < 2; reim++) { 

	spinor[spin4][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;

      }
    }
  } // Color

  spinor_array* sp = &spinor;
  halfspinor_array *s3 = &hspinor1;

  // Copy from code 
    _sse_42_1_gamma0_plus(*sp);
    _sse_store((*s3)[0]);

    _sse_42_2_gamma0_plus(*sp);
    _sse_store((*s3)[1]);


  // My version
  decomp_gamma0_plus(spinor, hspinor2);

  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	double diff = hspinor1[spin2][col][reim] - hspinor2[spin2][col][reim];
	diff /= (double)(3*2*2);
	//	QDPIO::cout << "   col="<<col<<" s="<<spin2<<" reim=" <<reim<< "  diff = " << diff << std::endl;
	assertion( fabs(diff)< 1.0e-17 );
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
	hspinor1[spin2][col][reim] = 0;
	hspinor2[spin2][col][reim] = 0;
      }
    }

    for(int spin4=0; spin4 < 4; spin4++) {
      for(int reim=0; reim < 2; reim++) { 

	spinor[spin4][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;

      }
    }
  } // Color

  spinor_array* sp = &spinor;
  halfspinor_array *s3 = &hspinor1;

  // Copy from code 
    _sse_42_1_gamma1_plus(*sp);
    _sse_store((*s3)[0]);

    _sse_42_2_gamma1_plus(*sp);
    _sse_store((*s3)[1]);


  // My version
  decomp_gamma1_plus(spinor, hspinor2);

  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	double diff = hspinor1[spin2][col][reim] - hspinor2[spin2][col][reim];
	diff /= (double)(3*2*2);
	//	QDPIO::cout << "   col="<<col<<" s="<<spin2<<" reim=" <<reim<< "  diff = " << diff << std::endl;
	assertion( fabs(diff)< 1.0e-17 );
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
	hspinor1[spin2][col][reim] = 0;
	hspinor2[spin2][col][reim] = 0;
      }
    }

    for(int spin4=0; spin4 < 4; spin4++) {
      for(int reim=0; reim < 2; reim++) { 

	spinor[spin4][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;

      }
    }
  } // Color

  spinor_array* sp = &spinor;
  halfspinor_array *s3 = &hspinor1;

  // Copy from code 

    _sse_42_1_gamma2_plus(*sp);
    _sse_store((*s3)[0]);

    _sse_42_2_gamma2_plus(*sp);
    _sse_store((*s3)[1]);


  // My version
  decomp_gamma2_plus(spinor, hspinor2);

  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	double diff = hspinor1[spin2][col][reim] - hspinor2[spin2][col][reim];
	diff /= (double)(3*2*2);
	//	QDPIO::cout << "   col="<<col<<" s="<<spin2<<" reim=" <<reim<< "  diff = " << diff << std::endl;
	assertion( fabs(diff)< 1.0e-17 );
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
	hspinor1[spin2][col][reim] = 0;
	hspinor2[spin2][col][reim] = 0;
      }
    }

    for(int spin4=0; spin4 < 4; spin4++) {
      for(int reim=0; reim < 2; reim++) { 

	spinor[spin4][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;

      }
    }
  } // Color

  spinor_array* sp = &spinor;
  halfspinor_array *s3 = &hspinor1;

    _sse_42_1_gamma3_plus(*sp);
    _sse_store((*s3)[0]);

    _sse_42_2_gamma3_plus(*sp);
    _sse_store((*s3)[1]);

  // My version
  decomp_gamma3_plus(spinor, hspinor2);

  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	double diff = hspinor1[spin2][col][reim] - hspinor2[spin2][col][reim];
	diff /= (double)(3*2*2);
	//	QDPIO::cout << "   col="<<col<<" s="<<spin2<<" reim=" <<reim<< "  diff = " << diff << std::endl;
	assertion( fabs(diff)< 1.0e-17 );
      }
    }
  }


}
