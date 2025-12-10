#include "testDecompHvv.h"
#include "sse64.h"
#include "types64.h"
#include "decomp_hvv.h"
#include "unittest.h"

using namespace QDP;
using namespace Assertions;

/* gamma 0 */

#define _sse_42_1_gamma0_minus(sp) \
      _sse_load((sp)[0]); \
      _sse_load_up((sp)[3]);\
      _sse_vector_i_mul();\
      _sse_vector_sub()

#define _sse_24_1_gamma0_minus_set() \
      _sse_store_up(rs[0]);\
      _sse_vector_i_mul_up();\
      _sse_store_up(rs[3]) 

#define _sse_24_1_gamma0_minus_add() \
      _sse_load(rs[0]);\
      _sse_vector_add();\
      _sse_store(rs[0]);\
      _sse_load(rs[3]);\
      _sse_vector_i_mul();\
      _sse_vector_add();\
      _sse_store(rs[3]) 
	  
#define _sse_42_2_gamma0_minus(sp) \
      _sse_load((sp)[1]);\
      _sse_load_up((sp)[2]);\
      _sse_vector_i_mul();\
      _sse_vector_sub()

#define _sse_24_2_gamma0_minus_set() \
	  _sse_store_up(rs[1]);\
      _sse_vector_i_mul_up();\
      _sse_store_up(rs[2])

#define _sse_24_2_gamma0_minus_add() \
	  _sse_load(rs[1]);\
      _sse_vector_add();\
      _sse_store(rs[1]);\
      _sse_load(rs[2]);\
      _sse_vector_i_mul();\
      _sse_vector_add();\
      _sse_store(rs[2])

#define _sse_42_1_gamma0_plus(sm) \
      _sse_load((sm)[0]);\
      _sse_load_up((sm)[3]);\
      _sse_vector_i_mul();\
      _sse_vector_add()

#define _sse_24_1_gamma0_plus_set() \
	  _sse_store_up(rs[0]);\
      _sse_vector_i_mul_neg_up();\
      _sse_store_up(rs[3])

#define _sse_24_1_gamma0_plus_add() \
	  _sse_load(rs[0]);\
      _sse_vector_add();\
      _sse_store(rs[0]);\
      _sse_load(rs[3]);\
      _sse_vector_i_mul();\
      _sse_vector_sub();\
      _sse_store(rs[3])

#define _sse_42_2_gamma0_plus(sm) \
	  _sse_load((sm)[1]);\
      _sse_load_up((sm)[2]);\
      _sse_vector_i_mul();\
      _sse_vector_add()

#define _sse_24_2_gamma0_plus_set() \
      _sse_store_up(rs[1]);\
      _sse_vector_i_mul_neg_up();  \
      _sse_store_up(rs[2])

#define _sse_24_2_gamma0_plus_add() \
       _sse_load(rs[1]);\
      _sse_vector_add();\
      _sse_store(rs[1]);\
      _sse_load(rs[2]);\
      _sse_vector_i_mul();  \
      _sse_vector_sub();\
      _sse_store(rs[2])



/* gamma 1 */


#define _sse_42_1_gamma1_minus(sp) \
      _sse_load((sp)[0]);\
      _sse_load_up((sp)[3]);\
      _sse_vector_add()

#define _sse_24_1_gamma1_minus() \
      _sse_load(rs[0]);\
      _sse_vector_add();\
      _sse_store(rs[0]);\
      _sse_load(rs[3]);\
      _sse_vector_add();\
      _sse_store(rs[3])
	  
#define _sse_42_2_gamma1_minus(sp) \
      _sse_load((sp)[1]);\
      _sse_load_up((sp)[2]);\
      _sse_vector_sub()

#define _sse_24_2_gamma1_minus() \
	  _sse_load(rs[1]);\
      _sse_vector_add();\
      _sse_store(rs[1]);\
      _sse_load(rs[2]);\
      _sse_vector_sub();\
      _sse_store(rs[2])

#define _sse_42_1_gamma1_plus(sm) \
      _sse_load((sm)[0]);\
      _sse_load_up((sm)[3]);\
      _sse_vector_sub()

#define _sse_24_1_gamma1_plus() \
      _sse_load(rs[0]);\
      _sse_vector_add();\
      _sse_store(rs[0]);\
      _sse_load(rs[3]);\
      _sse_vector_sub();\
      _sse_store(rs[3])

#define _sse_42_2_gamma1_plus(sm) \
	  _sse_load((sm)[1]);\
      _sse_load_up((sm)[2]);\
      _sse_vector_add()

#define _sse_24_2_gamma1_plus() \
       _sse_load(rs[1]);\
      _sse_vector_add();\
      _sse_store(rs[1]);\
      _sse_load(rs[2]);\
      _sse_vector_add();\
      _sse_store(rs[2])





/* gamma 2 */


#define _sse_42_1_gamma2_minus(sp) \
      _sse_load((sp)[0]);\
      _sse_load_up((sp)[2]);\
      _sse_vector_i_mul();\
      _sse_vector_sub()

#define _sse_24_1_gamma2_minus() \
      _sse_load(rs[0]);\
      _sse_vector_add();\
      _sse_store(rs[0]);\
      _sse_load(rs[2]);\
      _sse_vector_i_mul();   \
      _sse_vector_add();\
      _sse_store(rs[2])
	  
#define _sse_42_2_gamma2_minus(sp) \
      _sse_load((sp)[1]);\
      _sse_load_up((sp)[3]);\
      _sse_vector_i_mul();\
      _sse_vector_add()

#define _sse_24_2_gamma2_minus() \
	   _sse_load(rs[1]);\
      _sse_vector_add();\
      _sse_store(rs[1]);\
      _sse_load(rs[3]);\
      _sse_vector_i_mul();   \
      _sse_vector_sub();\
      _sse_store(rs[3])

#define _sse_42_1_gamma2_plus(sm) \
      _sse_load((sm)[0]);\
      _sse_load_up((sm)[2]);\
      _sse_vector_i_mul();\
      _sse_vector_add()

#define _sse_24_1_gamma2_plus() \
      _sse_load(rs[0]);\
      _sse_vector_add();\
      _sse_store(rs[0]);\
      _sse_load(rs[2]);\
      _sse_vector_i_mul();   \
      _sse_vector_sub();\
      _sse_store(rs[2]);

#define _sse_42_2_gamma2_plus(sm) \
	  _sse_load((sm)[1]);\
      _sse_load_up((sm)[3]);\
      _sse_vector_i_mul();\
      _sse_vector_sub()

#define _sse_24_2_gamma2_plus() \
      _sse_load(rs[1]);\
      _sse_vector_add();\
      _sse_store(rs[1]);\
      _sse_load(rs[3]);\
      _sse_vector_i_mul();     \
      _sse_vector_add();\
      _sse_store(rs[3])





/* gamma 3 */
#define _sse_42_1_gamma3_minus(sp) \
	  _sse_load((sp)[0]); \
	  _sse_load_up((sp)[2]); \
      _sse_vector_sub()

#define _sse_24_1_gamma3_minus_set() \
	  _sse_load(rs[0]);\
      _sse_vector_add();\
       _sse_store((*rn)[0]);\
      _sse_load(rs[2]);\
      _sse_vector_sub();\
       _sse_store((*rn)[2])

#define _sse_24_1_gamma3_minus_add() \
	  _sse_load(rs[0]);\
      _sse_vector_add();\
       _sse_store(rs[0]);\
      _sse_load(rs[2]);\
      _sse_vector_sub();\
       _sse_store(rs[2])
	  
#define _sse_42_2_gamma3_minus(sp) \
      _sse_load((sp)[1]);\
      _sse_load_up((sp)[3]);\
      _sse_vector_sub()

#define _sse_24_2_gamma3_minus_set() \
	  _sse_load(rs[1]);\
      _sse_vector_add();\
       _sse_store((*rn)[1]);\
      _sse_load(rs[3]);\
      _sse_vector_sub();\
      _sse_store((*rn)[3])

#define _sse_24_2_gamma3_minus_add() \
	  _sse_load(rs[1]);\
      _sse_vector_add();\
       _sse_store(rs[1]);\
      _sse_load(rs[3]);\
      _sse_vector_sub();\
      _sse_store(rs[3])

#define _sse_42_1_gamma3_plus(sm) \
      _sse_load((sm)[0]);\
      _sse_load_up((sm)[2]);\
      _sse_vector_add()

#define _sse_24_1_gamma3_plus_set() \
	  _sse_load(rs[0]);\
      _sse_vector_add();\
       _sse_store((*rn)[0]);\
      _sse_load(rs[2]);\
      _sse_vector_add();\
      _sse_store((*rn)[2])

#define _sse_24_1_gamma3_plus_add() \
	  _sse_load(rs[0]);\
      _sse_vector_add();\
       _sse_store(rs[0]);\
      _sse_load(rs[2]);\
      _sse_vector_add();\
      _sse_store(rs[2])

#define _sse_42_2_gamma3_plus(sm) \
	  _sse_load((sm)[1]);\
      _sse_load_up((sm)[3]);\
      _sse_vector_add()

#define _sse_24_2_gamma3_plus_set() \
       _sse_load(rs[1]);\
      _sse_vector_add();\
       _sse_store((*rn)[1]);\
      _sse_load(rs[3]);\
      _sse_vector_add();\
       _sse_store((*rn)[3])

#define _sse_24_2_gamma3_plus_add() \
       _sse_load(rs[1]);\
       _sse_vector_add();\
       _sse_store(rs[1]);\
       _sse_load(rs[3]);\
       _sse_vector_add();\
       _sse_store(rs[3]);\



void testDecompHvv0Plus::run(void) 
{
  spinor_array spinor ALIGN;
  halfspinor_array hspinor1 ALIGN,  hspinor2 ALIGN;
  u_mat_array matrix ALIGN;

  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	hspinor1[spin2][col][reim] = 0;
	hspinor2[spin2][col][reim] = 0;
      }
    }
  }

  for(int col=0; col < 3; col++) { 
    for(int spin4=0; spin4 < 4; spin4++) {
      for(int reim=0; reim < 2; reim++) { 

	spinor[spin4][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
	
      }
    }
  } // Color

  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) { 
      for(int reim=0; reim < 2; reim++) { 
	matrix[col1][col2][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }
    }
  }

  spinor_array* sm = &spinor;
  u_mat_array*  um = &matrix;
  halfspinor_array *s3 = &hspinor1;

  // Copy from code 
  _sse_42_1_gamma0_plus(*sm);
  _sse_su3_inverse_multiply(*um);
  _sse_store_up((*s3)[0]);

  _sse_42_2_gamma0_plus(*sm);
  _sse_su3_inverse_multiply(*um);
  _sse_store_up((*s3)[1]);


  // My version
  decomp_hvv_gamma0_plus(spinor, matrix, hspinor2);

  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	float diff = hspinor1[spin2][col][reim] - hspinor2[spin2][col][reim];
#if 0
	QDPIO::cout << " col=" << col << " spin= " << spin2 
		    << " rei=" << reim << " diff = " << diff
		    <<std::endl;
#endif
	diff /= (double)(3*2*2);
	assertion( fabs(diff) < 1.0e-17 );
      }
    }
  }
}


void testDecompHvv1Plus::run(void) 
{
  spinor_array spinor;
  halfspinor_array hspinor1, hspinor2;
  u_mat_array matrix;

  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	hspinor1[spin2][col][reim] = 0;
	hspinor2[spin2][col][reim] = 0;
      }
    }
  }

  for(int col=0; col < 3; col++) { 
    for(int spin4=0; spin4 < 4; spin4++) {
      for(int reim=0; reim < 2; reim++) { 

	spinor[spin4][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
	
      }
    }
  } // Color

  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) { 
      for(int reim=0; reim < 2; reim++) { 
	matrix[col1][col2][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }
    }
  }

  spinor_array* sm = &spinor;
  u_mat_array*  um = &matrix;
  halfspinor_array *s3 = &hspinor1;

  // Copy from code 
  _sse_42_1_gamma1_plus(*sm);
  _sse_su3_inverse_multiply(*um);
  _sse_store_up((*s3)[0]);
  
  _sse_42_2_gamma1_plus(*sm);
  _sse_su3_inverse_multiply(*um);
  _sse_store_up((*s3)[1]);

  // My version
  decomp_hvv_gamma1_plus(spinor, matrix, hspinor2);

  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	float diff = hspinor1[spin2][col][reim] - hspinor2[spin2][col][reim];

	diff /= (double)(3*2*2); 
	assertion( fabs(diff) < 1.0e-17 );
      }
    }
  }
}


void testDecompHvv2Plus::run(void) 
{  
  spinor_array spinor;
  halfspinor_array hspinor1, hspinor2;
  u_mat_array matrix;

  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	hspinor1[spin2][col][reim] = 0;
	hspinor2[spin2][col][reim] = 0;
      }
    }
  }

  for(int col=0; col < 3; col++) { 
    for(int spin4=0; spin4 < 4; spin4++) {
      for(int reim=0; reim < 2; reim++) { 

	spinor[spin4][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
	
      }
    }
  } // Color

  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) { 
      for(int reim=0; reim < 2; reim++) { 
	matrix[col1][col2][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }
    }
  }

  spinor_array* sm = &spinor;
  u_mat_array*  um = &matrix;
  halfspinor_array *s3 = &hspinor1;

  // Copy from code 
   _sse_42_1_gamma2_plus(*sm);
    _sse_su3_inverse_multiply(*um);
    _sse_store_up((*s3)[0]);

    _sse_42_2_gamma2_plus(*sm);
    _sse_su3_inverse_multiply(*um);
    _sse_store_up((*s3)[1]);

  // My version
  decomp_hvv_gamma2_plus(spinor, matrix, hspinor2);

  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	float diff = hspinor1[spin2][col][reim] - hspinor2[spin2][col][reim];
	diff /= (double)(3*2*2); 
	assertion( fabs(diff) < 1.0e-17 );
      }
    }
  }

}



void testDecompHvv3Plus::run(void) 
{
  spinor_array spinor;
  halfspinor_array hspinor1, hspinor2;
  u_mat_array matrix;

  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	hspinor1[spin2][col][reim] = 0;
	hspinor2[spin2][col][reim] = 0;
      }
    }
  }

  for(int col=0; col < 3; col++) { 
    for(int spin4=0; spin4 < 4; spin4++) {
      for(int reim=0; reim < 2; reim++) { 

	spinor[spin4][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
	
      }
    }
  } // Color

  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) { 
      for(int reim=0; reim < 2; reim++) { 
	matrix[col1][col2][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }
    }
  }

  spinor_array* sm = &spinor;
  u_mat_array*  um = &matrix;
  halfspinor_array *s3 = &hspinor1;

  // Copy from code 
    _sse_42_1_gamma3_plus(*sm);
    _sse_su3_inverse_multiply(*um);
    _sse_store_up((*s3)[0]);

    _sse_42_2_gamma3_plus(*sm);
    _sse_su3_inverse_multiply(*um);
    _sse_store_up((*s3)[1]);


  // My version
  decomp_hvv_gamma3_plus(spinor, matrix, hspinor2);

  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	float diff = hspinor1[spin2][col][reim] - hspinor2[spin2][col][reim];
	diff /= (double)(3*2*2); 
	assertion( fabs(diff) < 1.0e-17 );
      }
    }
  }

}


void testDecompHvv0Minus::run(void) 
{
  spinor_array spinor;
  halfspinor_array hspinor1, hspinor2;
  u_mat_array matrix;

  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	hspinor1[spin2][col][reim] = 0;
	hspinor2[spin2][col][reim] = 0;
      }
    }
  }

  for(int col=0; col < 3; col++) { 
    for(int spin4=0; spin4 < 4; spin4++) {
      for(int reim=0; reim < 2; reim++) { 

	spinor[spin4][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
	
      }
    }
  } // Color

  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) { 
      for(int reim=0; reim < 2; reim++) { 
	matrix[col1][col2][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }
    }
  }

  spinor_array* sm = &spinor;
  u_mat_array*  um = &matrix;
  halfspinor_array *s3 = &hspinor1;

  // Copy from code 
    _sse_42_1_gamma0_minus(*sm);
    _sse_su3_inverse_multiply(*um);
    _sse_store_up((*s3)[0]);

    _sse_42_2_gamma0_minus(*sm);
    _sse_su3_inverse_multiply(*um);
    _sse_store_up((*s3)[1]);

  // My version
  decomp_hvv_gamma0_minus(spinor, matrix, hspinor2);

  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	float diff = hspinor1[spin2][col][reim] - hspinor2[spin2][col][reim];
	diff /= (double)(3*2*2); 
	assertion( fabs(diff) < 1.0e-17 );
      }
    }
  }

}


void testDecompHvv1Minus::run(void) 
{
  spinor_array spinor;
  halfspinor_array hspinor1, hspinor2;
  u_mat_array matrix;

  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	hspinor1[spin2][col][reim] = 0;
	hspinor2[spin2][col][reim] = 0;
      }
    }
  }

  for(int col=0; col < 3; col++) { 
    for(int spin4=0; spin4 < 4; spin4++) {
      for(int reim=0; reim < 2; reim++) { 

	spinor[spin4][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
	
      }
    }
  } // Color

  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) { 
      for(int reim=0; reim < 2; reim++) { 
	matrix[col1][col2][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }
    }
  }

  spinor_array* sm = &spinor;
  u_mat_array*  um = &matrix;
  halfspinor_array *s3 = &hspinor1;

  // Copy from code 
    _sse_42_1_gamma1_minus(*sm);
    _sse_su3_inverse_multiply(*um);
    _sse_store_up((*s3)[0]);

    _sse_42_2_gamma1_minus(*sm);
    _sse_su3_inverse_multiply(*um);
    _sse_store_up((*s3)[1]);

  // My version
  decomp_hvv_gamma1_minus(spinor, matrix, hspinor2);

  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	float diff = hspinor1[spin2][col][reim] - hspinor2[spin2][col][reim];
	diff /= (double)(3*2*2); 
	assertion( fabs(diff) < 1.0e-17 );
      }
    }
  }
}

void testDecompHvv2Minus::run(void) 
{
  spinor_array spinor;
  halfspinor_array hspinor1, hspinor2;
  u_mat_array matrix;

  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	hspinor1[spin2][col][reim] = 0;
	hspinor2[spin2][col][reim] = 0;
      }
    }
  }

  for(int col=0; col < 3; col++) { 
    for(int spin4=0; spin4 < 4; spin4++) {
      for(int reim=0; reim < 2; reim++) { 

	spinor[spin4][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
	
      }
    }
  } // Color

  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) { 
      for(int reim=0; reim < 2; reim++) { 
	matrix[col1][col2][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }
    }
  }

  spinor_array* sm = &spinor;
  u_mat_array*  um = &matrix;
  halfspinor_array *s3 = &hspinor1;

  // Copy from code 
    _sse_42_1_gamma2_minus(*sm);
    _sse_su3_inverse_multiply(*um);
    _sse_store_up((*s3)[0]);

    _sse_42_2_gamma2_minus(*sm);
    _sse_su3_inverse_multiply(*um);
    _sse_store_up((*s3)[1]);

  // My version
  decomp_hvv_gamma2_minus(spinor, matrix, hspinor2);

  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	float diff = hspinor1[spin2][col][reim] - hspinor2[spin2][col][reim];
	diff /= (double)(3*2*2); 
	assertion( fabs(diff) < 1.0e-17 );
      }
    }
  }
}

void testDecompHvv3Minus::run(void) 
{
  spinor_array spinor;
  halfspinor_array hspinor1, hspinor2;
  u_mat_array matrix;

  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	hspinor1[spin2][col][reim] = 0;
	hspinor2[spin2][col][reim] = 0;
      }
    }
  }

  for(int col=0; col < 3; col++) { 
    for(int spin4=0; spin4 < 4; spin4++) {
      for(int reim=0; reim < 2; reim++) { 

	spinor[spin4][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
	
      }
    }
  } // Color

  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) { 
      for(int reim=0; reim < 2; reim++) { 
	matrix[col1][col2][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }
    }
  }

  spinor_array* sm = &spinor;
  u_mat_array*  um = &matrix;
  halfspinor_array *s3 = &hspinor1;

  // Copy from code 
  _sse_42_1_gamma3_minus(*sm);
  _sse_su3_inverse_multiply(*um);
  _sse_store_up((*s3)[0]);
  
  _sse_42_2_gamma3_minus(*sm);
  _sse_su3_inverse_multiply(*um);
  _sse_store_up((*s3)[1]);
  
  // My version
  decomp_hvv_gamma3_minus(spinor, matrix, hspinor2);

  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	float diff = hspinor1[spin2][col][reim] - hspinor2[spin2][col][reim];
	diff /= (double)(3*2*2); 
	assertion( fabs(diff) < 1.0e-17 );
      }
    }
  }
}
