#include "testSiteDslash.h"

#include "sse_align.h"
#include "types32.h"
#include "sse32.h"

#include <cmath>
#include "site_dslash_32bit_scalar.h"

using namespace QDP;
using namespace Assertions;

  /* gamma 0 */

#define _sse_42_gamma0_minus()   _sse_vector_xch_i_sub()

#define _sse_42_gamma0_plus()     _sse_vector_xch_i_add()

#define _sse_24_gamma0_minus_set()  _sse_vector_xch_i_mul_up()
 
#define _sse_24_gamma0_plus_set()   _sse_vector_xch_i_mul_neg_up()

#define _sse_24_gamma0_minus_add() _sse_vector_xch_i_add()

#define _sse_24_gamma0_plus_add() _sse_vector_xch_i_sub()


/* gamma 1 */

#define _sse_42_gamma1_minus()  \
								_sse_vector_xch();\
								_sse_vector_addsub()
#define _sse_42_gamma1_plus()  \
								_sse_vector_xch();\
								_sse_vector_subadd()
#define _sse_24_gamma1_minus()  \
								_sse_vector_xch();\
								_sse_vector_subadd()
#define _sse_24_gamma1_plus()  \
								_sse_vector_xch();\
								_sse_vector_addsub()



/* gamma 2 */

#define _sse_42_gamma2_minus()   _sse_vector_i_subadd()

#define _sse_42_gamma2_plus()     _sse_vector_i_addsub()

#define _sse_24_gamma2_minus() _sse_vector_i_addsub()

#define _sse_24_gamma2_plus() _sse_vector_i_subadd()

/* gamma 3 */


#define _sse_42_gamma3_minus()   _sse_vector_sub()

#define _sse_42_gamma3_plus()     _sse_vector_add()

#define _sse_24_gamma3_minus() _sse_vector_sub()

#define _sse_24_gamma3_plus() _sse_vector_add()

#define _sse_24_gamma3_minus_rows12() _sse_vector_add()

#define _sse_24_gamma3_plus_rows12() _sse_vector_add()


void 
testSiteDslash0PlusForward::run()
{
  spinor_array input_spinor ALIGN;
  u_mat_array gauge_field ALIGN;
  halfspinor_array upper_sum1 ALIGN, upper_sum2 ALIGN;
  halfspinor_array lower_sum1 ALIGN, lower_sum2 ALIGN;

  // Fill inpur spinor with noise 
  for(int spin4=0; spin4 < 4; spin4++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	input_spinor[spin4][col][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
      }
    }
  }

  for(int spin2=0; spin2 < 2; spin2++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	
	upper_sum1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	lower_sum1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	upper_sum2[col][spin2][reim] = upper_sum1[col][spin2][reim];
	lower_sum2[col][spin2][reim] = lower_sum1[col][spin2][reim];
      }
    }
  }

  // Fill the gauge field with noise 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) {
      for(int reim=0; reim < 2; reim++) { 
	gauge_field[col1][col2][reim] =  (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
      }
    }
  }

  /* Stuff goes here */
  _sse_pair_load(input_spinor[0],input_spinor[1]);
  _sse_pair_load_up(input_spinor[2],input_spinor[3]);
  _sse_42_gamma0_minus();
  _sse_su3_multiply(gauge_field);
  
  _sse_vector_store_up(upper_sum1);
  _sse_24_gamma0_minus_set();
  _sse_vector_store_up(lower_sum1);

  dslash_plus_dir0_forward(input_spinor, gauge_field, upper_sum2, lower_sum2);

  // Check on the results
  for(int spin2=0; spin2 < 2; spin2++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 

	float diff = upper_sum2[col][spin2][reim] - upper_sum1[col][spin2][reim];
	diff /= (float)(2*3*2);

#if 0
	QDPIO::cout << "\t sp=" << spin2 
		    << " co=" << col 
		    << " re=" << reim 
		    << " diff= " << diff << std::endl;
#endif
	assertion( fabs(diff) < 1.0e-9 );

	diff = lower_sum2[col][spin2][reim] - lower_sum1[col][spin2][reim];
	diff /= (float)(2*3*2);
#if 0
	QDPIO::cout << "\t sp=" << spin2 
		    << " co=" << col 
		    << " re=" << reim 
		    << " diff2 = " << diff << std::endl;
#endif
	assertion( fabs(diff) < 1.0e-9 );

      }
    }
  }


}

void 
testSiteDslash0PlusBackwardAdd::run()
{
  spinor_array input_spinor ALIGN;
  u_mat_array gauge_field ALIGN;
  halfspinor_array upper_sum1 ALIGN, upper_sum2 ALIGN;
  halfspinor_array lower_sum1 ALIGN, lower_sum2 ALIGN;

  // Fill inpur spinor with noise 
  for(int spin4=0; spin4 < 4; spin4++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	input_spinor[spin4][col][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
      }
    }
  }

  for(int spin2=0; spin2 < 2; spin2++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	
	upper_sum1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	lower_sum1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	upper_sum2[col][spin2][reim] = upper_sum1[col][spin2][reim];
	lower_sum2[col][spin2][reim] = lower_sum1[col][spin2][reim];
      }
    }
  }

  // Fill the gauge field with noise 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) {
      for(int reim=0; reim < 2; reim++) { 
	gauge_field[col1][col2][reim] =  (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
      }
    }
  }

  /* Stuff goes here */
  _sse_pair_load(input_spinor[0], input_spinor[1]);
  _sse_pair_load_up(input_spinor[2],input_spinor[3]);
  _sse_42_gamma0_plus();
  _sse_su3_inverse_multiply(gauge_field);
  _sse_vector_load(upper_sum1);
  _sse_vector_add();
  _sse_vector_store(upper_sum1);
  _sse_vector_load(lower_sum1);
  _sse_24_gamma0_plus_add();
  _sse_vector_store(lower_sum1);  



  dslash_plus_dir0_backward_add(input_spinor, gauge_field, upper_sum2, lower_sum2);

  // Check on the results
  for(int spin2=0; spin2 < 2; spin2++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 

	float diff = upper_sum2[col][spin2][reim] - upper_sum1[col][spin2][reim];
	diff /= (float)(2*3*2);

#if 0
	QDPIO::cout << "\t sp=" << spin2 
		    << " co=" << col 
		    << " re=" << reim 
		    << " diff= " << diff << std::endl;
#endif
	assertion( fabs(diff) < 1.0e-9 );

	diff = lower_sum2[col][spin2][reim] - lower_sum1[col][spin2][reim];
	diff /= (float)(2*3*2);
#if 0
	QDPIO::cout << "\t sp=" << spin2 
		    << " co=" << col 
		    << " re=" << reim 
		    << " diff2 = " << diff << std::endl;
#endif
	assertion( fabs(diff) < 1.0e-9 );

      }
    }
  }


}

void 
testSiteDslash1PlusForwardAdd::run()
{
  spinor_array input_spinor ALIGN;
  u_mat_array gauge_field ALIGN;
  halfspinor_array upper_sum1 ALIGN, upper_sum2 ALIGN;
  halfspinor_array lower_sum1 ALIGN, lower_sum2 ALIGN;

  // Fill inpur spinor with noise 
  for(int spin4=0; spin4 < 4; spin4++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	input_spinor[spin4][col][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
      }
    }
  }

  for(int spin2=0; spin2 < 2; spin2++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	
	upper_sum1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	lower_sum1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	upper_sum2[col][spin2][reim] = upper_sum1[col][spin2][reim];
	lower_sum2[col][spin2][reim] = lower_sum1[col][spin2][reim];
      }
    }
  }

  // Fill the gauge field with noise 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) {
      for(int reim=0; reim < 2; reim++) { 
	gauge_field[col1][col2][reim] =  (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
      }
    }
  }

  /* Stuff goes here */
  _sse_pair_load(input_spinor[0], input_spinor[1]);
  _sse_pair_load_up(input_spinor[2],input_spinor[3]);
  _sse_42_gamma1_minus();
  _sse_su3_multiply(gauge_field);
  _sse_vector_load(upper_sum1);
  _sse_vector_add();
  _sse_vector_store(upper_sum1);
  _sse_vector_load(lower_sum1);
  _sse_24_gamma1_minus();
  _sse_vector_store(lower_sum1);  



  dslash_plus_dir1_forward_add(input_spinor, gauge_field, upper_sum2, lower_sum2);

  // Check on the results
  for(int spin2=0; spin2 < 2; spin2++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 

	float diff = upper_sum2[col][spin2][reim] - upper_sum1[col][spin2][reim];
	diff /= (float)(2*3*2);

#if 0
	QDPIO::cout << "\t sp=" << spin2 
		    << " co=" << col 
		    << " re=" << reim 
		    << " diff= " << diff << std::endl;
#endif
	assertion( fabs(diff) < 1.0e-9 );

	diff = lower_sum2[col][spin2][reim] - lower_sum1[col][spin2][reim];
	diff /= (float)(2*3*2);
#if 0
	QDPIO::cout << "\t sp=" << spin2 
		    << " co=" << col 
		    << " re=" << reim 
		    << " diff2 = " << diff << std::endl;
#endif
	assertion( fabs(diff) < 1.0e-9 );

      }
    }
  }
 
}

void 
testSiteDslash1PlusBackwardAdd::run()
{
  spinor_array input_spinor ALIGN;
  u_mat_array gauge_field ALIGN;
  halfspinor_array upper_sum1 ALIGN, upper_sum2 ALIGN;
  halfspinor_array lower_sum1 ALIGN, lower_sum2 ALIGN;

  // Fill inpur spinor with noise 
  for(int spin4=0; spin4 < 4; spin4++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	input_spinor[spin4][col][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
      }
    }
  }

  for(int spin2=0; spin2 < 2; spin2++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	
	upper_sum1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	lower_sum1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	upper_sum2[col][spin2][reim] = upper_sum1[col][spin2][reim];
	lower_sum2[col][spin2][reim] = lower_sum1[col][spin2][reim];
      }
    }
  }

  // Fill the gauge field with noise 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) {
      for(int reim=0; reim < 2; reim++) { 
	gauge_field[col1][col2][reim] =  (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
      }
    }
  }

  /* Stuff goes here */
  _sse_pair_load(input_spinor[0], input_spinor[1]);
  _sse_pair_load_up(input_spinor[2],input_spinor[3]);
  _sse_42_gamma1_plus();
  _sse_su3_inverse_multiply(gauge_field);
  _sse_vector_load(upper_sum1);
  _sse_vector_add();
  _sse_vector_store(upper_sum1);
  _sse_vector_load(lower_sum1);
  _sse_24_gamma1_plus();
  _sse_vector_store(lower_sum1);  



  dslash_plus_dir1_backward_add(input_spinor, gauge_field, upper_sum2, lower_sum2);

  // Check on the results
  for(int spin2=0; spin2 < 2; spin2++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 

	float diff = upper_sum2[col][spin2][reim] - upper_sum1[col][spin2][reim];
	diff /= (float)(2*3*2);

#if 0
	QDPIO::cout << "\t sp=" << spin2 
		    << " co=" << col 
		    << " re=" << reim 
		    << " diff= " << diff << std::endl;
#endif
	assertion( fabs(diff) < 1.0e-9 );

	diff = lower_sum2[col][spin2][reim] - lower_sum1[col][spin2][reim];
	diff /= (float)(2*3*2);
#if 0
	QDPIO::cout << "\t sp=" << spin2 
		    << " co=" << col 
		    << " re=" << reim 
		    << " diff2 = " << diff << std::endl;
#endif
	assertion( fabs(diff) < 1.0e-9 );

      }
    }
  }
}


void 
testSiteDslash2PlusForwardAdd::run()
{

  spinor_array input_spinor ALIGN;
  u_mat_array gauge_field ALIGN;
  halfspinor_array upper_sum1 ALIGN, upper_sum2 ALIGN;
  halfspinor_array lower_sum1 ALIGN, lower_sum2 ALIGN;

  // Fill inpur spinor with noise 
  for(int spin4=0; spin4 < 4; spin4++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	input_spinor[spin4][col][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
      }
    }
  }

  for(int spin2=0; spin2 < 2; spin2++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	
	upper_sum1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	lower_sum1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	upper_sum2[col][spin2][reim] = upper_sum1[col][spin2][reim];
	lower_sum2[col][spin2][reim] = lower_sum1[col][spin2][reim];
      }
    }
  }

  // Fill the gauge field with noise 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) {
      for(int reim=0; reim < 2; reim++) { 
	gauge_field[col1][col2][reim] =  (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
      }
    }
  }

  /* Stuff goes here */
  _sse_pair_load(input_spinor[0], input_spinor[1]);
  _sse_pair_load_up(input_spinor[2],input_spinor[3]);
  _sse_42_gamma2_minus();
  _sse_su3_multiply(gauge_field);
  _sse_vector_load(upper_sum1);
  _sse_vector_add();
  _sse_vector_store(upper_sum1);
  _sse_vector_load(lower_sum1);
  _sse_24_gamma2_minus();
  _sse_vector_store(lower_sum1);  



  dslash_plus_dir2_forward_add(input_spinor, gauge_field, upper_sum2, lower_sum2);

  // Check on the results
  for(int spin2=0; spin2 < 2; spin2++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 

	float diff = upper_sum2[col][spin2][reim] - upper_sum1[col][spin2][reim];
	diff /= (float)(2*3*2);

#if 0
	QDPIO::cout << "\t sp=" << spin2 
		    << " co=" << col 
		    << " re=" << reim 
		    << " diff= " << diff << std::endl;
#endif
	assertion( fabs(diff) < 1.0e-9 );

	diff = lower_sum2[col][spin2][reim] - lower_sum1[col][spin2][reim];
	diff /= (float)(2*3*2);
#if 0
	QDPIO::cout << "\t sp=" << spin2 
		    << " co=" << col 
		    << " re=" << reim 
		    << " diff2 = " << diff << std::endl;
#endif
	assertion( fabs(diff) < 1.0e-9 );

      }
    }
  }
 
}

void 
testSiteDslash2PlusBackwardAdd::run()
{
  spinor_array input_spinor ALIGN;
  u_mat_array gauge_field ALIGN;
  halfspinor_array upper_sum1 ALIGN, upper_sum2 ALIGN;
  halfspinor_array lower_sum1 ALIGN, lower_sum2 ALIGN;

  // Fill inpur spinor with noise 
  for(int spin4=0; spin4 < 4; spin4++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	input_spinor[spin4][col][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
      }
    }
  }

  for(int spin2=0; spin2 < 2; spin2++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	
	upper_sum1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	lower_sum1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	upper_sum2[col][spin2][reim] = upper_sum1[col][spin2][reim];
	lower_sum2[col][spin2][reim] = lower_sum1[col][spin2][reim];
      }
    }
  }

  // Fill the gauge field with noise 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) {
      for(int reim=0; reim < 2; reim++) { 
	gauge_field[col1][col2][reim] =  (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
      }
    }
  }

  /* Stuff goes here */
  _sse_pair_load(input_spinor[0], input_spinor[1]);
  _sse_pair_load_up(input_spinor[2],input_spinor[3]);
  _sse_42_gamma2_plus();
  _sse_su3_inverse_multiply(gauge_field);
  _sse_vector_load(upper_sum1);
  _sse_vector_add();
  _sse_vector_store(upper_sum1);
  _sse_vector_load(lower_sum1);
  _sse_24_gamma2_plus();
  _sse_vector_store(lower_sum1);  



  dslash_plus_dir2_backward_add(input_spinor, gauge_field, upper_sum2, lower_sum2);

  // Check on the results
  for(int spin2=0; spin2 < 2; spin2++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 

	float diff = upper_sum2[col][spin2][reim] - upper_sum1[col][spin2][reim];
	diff /= (float)(2*3*2);

#if 0
	QDPIO::cout << "\t sp=" << spin2 
		    << " co=" << col 
		    << " re=" << reim 
		    << " diff= " << diff << std::endl;
#endif
	assertion( fabs(diff) < 1.0e-9 );

	diff = lower_sum2[col][spin2][reim] - lower_sum1[col][spin2][reim];
	diff /= (float)(2*3*2);
#if 0
	QDPIO::cout << "\t sp=" << spin2 
		    << " co=" << col 
		    << " re=" << reim 
		    << " diff2 = " << diff << std::endl;
#endif
	assertion( fabs(diff) < 1.0e-9 );

      }
    }
  }

}


void 
testSiteDslash3PlusForwardAdd::run()
{

  spinor_array input_spinor ALIGN;
  u_mat_array gauge_field ALIGN;
  halfspinor_array upper_sum1 ALIGN, upper_sum2 ALIGN;
  halfspinor_array lower_sum1 ALIGN, lower_sum2 ALIGN;

  // Fill inpur spinor with noise 
  for(int spin4=0; spin4 < 4; spin4++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	input_spinor[spin4][col][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
      }
    }
  }

  for(int spin2=0; spin2 < 2; spin2++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	
	upper_sum1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	lower_sum1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	upper_sum2[col][spin2][reim] = upper_sum1[col][spin2][reim];
	lower_sum2[col][spin2][reim] = lower_sum1[col][spin2][reim];
      }
    }
  }

  // Fill the gauge field with noise 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) {
      for(int reim=0; reim < 2; reim++) { 
	gauge_field[col1][col2][reim] =  (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
      }
    }
  }

  /* Stuff goes here */
  _sse_pair_load(input_spinor[0], input_spinor[1]);
  _sse_pair_load_up(input_spinor[2],input_spinor[3]);
  _sse_42_gamma3_minus();
  _sse_su3_multiply(gauge_field);
  _sse_vector_load(upper_sum1);
  _sse_vector_add();
  _sse_vector_store(upper_sum1);
  _sse_vector_load(lower_sum1);
  _sse_24_gamma3_minus();
  _sse_vector_store(lower_sum1);  



  dslash_plus_dir3_forward_add(input_spinor, gauge_field, upper_sum2, lower_sum2);

  // Check on the results
  for(int spin2=0; spin2 < 2; spin2++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 

	float diff = upper_sum2[col][spin2][reim] - upper_sum1[col][spin2][reim];
	diff /= (float)(2*3*2);

#if 0
	QDPIO::cout << "\t sp=" << spin2 
		    << " co=" << col 
		    << " re=" << reim 
		    << " diff= " << diff << std::endl;
#endif
	assertion( fabs(diff) < 1.0e-9 );

	diff = lower_sum2[col][spin2][reim] - lower_sum1[col][spin2][reim];
	diff /= (float)(2*3*2);
#if 0
	QDPIO::cout << "\t sp=" << spin2 
		    << " co=" << col 
		    << " re=" << reim 
		    << " diff2 = " << diff << std::endl;
#endif
	assertion( fabs(diff) < 1.0e-9 );

      }
    }
  }

}



void 
testSiteDslash3PlusBackwardAddStore::run()
{

  spinor_array input_spinor ALIGN;
  spinor_array output_spinor ALIGN, output_spinor2;
  u_mat_array gauge_field ALIGN;
  halfspinor_array upper_sum1 ALIGN, upper_sum2 ALIGN;
  halfspinor_array lower_sum1 ALIGN, lower_sum2 ALIGN;
  
  // Fill inpur spinor with noise 
  for(int spin4=0; spin4 < 4; spin4++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	input_spinor[spin4][col][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
      }
    }
  }

  for(int spin2=0; spin2 < 2; spin2++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	
	upper_sum1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	lower_sum1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	upper_sum2[col][spin2][reim] = upper_sum1[col][spin2][reim];
	lower_sum2[col][spin2][reim] = lower_sum1[col][spin2][reim];
      }
    }
  }

  // Fill the gauge field with noise 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) {
      for(int reim=0; reim < 2; reim++) { 
	gauge_field[col1][col2][reim] =  (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
      }
    }
  }

  /* Stuff goes here */
  _sse_pair_load(input_spinor[0], input_spinor[1]);
  _sse_pair_load_up(input_spinor[2],input_spinor[3]);
  _sse_42_gamma3_plus();
  _sse_su3_inverse_multiply(gauge_field);
  _sse_vector_load(upper_sum1);
  _sse_24_gamma3_plus_rows12();
  _sse_pair_store(output_spinor[0], output_spinor[1]);
  _sse_vector_load(lower_sum1);
  _sse_24_gamma3_plus();
  _sse_pair_store(output_spinor[2], output_spinor[3]);



  dslash_plus_dir3_backward_add_store(input_spinor, gauge_field, upper_sum2, lower_sum2, output_spinor2);

  // Check on the results
  for(int spin4=0; spin4 < 4; spin4++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 

	float diff = output_spinor[spin4][col][reim] - output_spinor2[spin4][col][reim];
	diff /= (float)(4*3*2);

#if 0
	QDPIO::cout << "\t sp=" << spin4
		    << " co=" << col 
		    << " re=" << reim 
		    << " diff= " << diff << std::endl;
#endif
	assertion( fabs(diff) < 1.0e-9 );


      }
    }
  }
}


void 
testSiteDslash0MinusForward::run()
{
  spinor_array input_spinor ALIGN;
  u_mat_array gauge_field ALIGN;
  halfspinor_array upper_sum1 ALIGN, upper_sum2 ALIGN;
  halfspinor_array lower_sum1 ALIGN, lower_sum2 ALIGN;

  // Fill inpur spinor with noise 
  for(int spin4=0; spin4 < 4; spin4++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	input_spinor[spin4][col][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
      }
    }
  }

  for(int spin2=0; spin2 < 2; spin2++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	
	upper_sum1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	lower_sum1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	upper_sum2[col][spin2][reim] = upper_sum1[col][spin2][reim];
	lower_sum2[col][spin2][reim] = lower_sum1[col][spin2][reim];
      }
    }
  }

  // Fill the gauge field with noise 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) {
      for(int reim=0; reim < 2; reim++) { 
	gauge_field[col1][col2][reim] =  (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
      }
    }
  }

  /* Stuff goes here */
  _sse_pair_load(input_spinor[0],input_spinor[1]);
  _sse_pair_load_up(input_spinor[2],input_spinor[3]);
  _sse_42_gamma0_plus();
  _sse_su3_multiply(gauge_field);
  
  _sse_vector_store_up(upper_sum1);
  _sse_24_gamma0_plus_set();
  _sse_vector_store_up(lower_sum1);

  dslash_minus_dir0_forward(input_spinor, gauge_field, upper_sum2, lower_sum2);

  // Check on the results
  for(int spin2=0; spin2 < 2; spin2++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 

	float diff = upper_sum2[col][spin2][reim] - upper_sum1[col][spin2][reim];
	diff /= (float)(2*3*2);

#if 0
	QDPIO::cout << "\t sp=" << spin2 
		    << " co=" << col 
		    << " re=" << reim 
		    << " diff= " << diff << std::endl;
#endif
	assertion( fabs(diff) < 1.0e-9 );

	diff = lower_sum2[col][spin2][reim] - lower_sum1[col][spin2][reim];
	diff /= (float)(2*3*2);
#if 0
	QDPIO::cout << "\t sp=" << spin2 
		    << " co=" << col 
		    << " re=" << reim 
		    << " diff2 = " << diff << std::endl;
#endif
	assertion( fabs(diff) < 1.0e-9 );

      }
    }
  }


}


void 
testSiteDslash0MinusBackwardAdd::run()
{
  spinor_array input_spinor ALIGN;
  u_mat_array gauge_field ALIGN;
  halfspinor_array upper_sum1 ALIGN, upper_sum2 ALIGN;
  halfspinor_array lower_sum1 ALIGN, lower_sum2 ALIGN;

  // Fill inpur spinor with noise 
  for(int spin4=0; spin4 < 4; spin4++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	input_spinor[spin4][col][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
      }
    }
  }

  for(int spin2=0; spin2 < 2; spin2++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	
	upper_sum1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	lower_sum1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	upper_sum2[col][spin2][reim] = upper_sum1[col][spin2][reim];
	lower_sum2[col][spin2][reim] = lower_sum1[col][spin2][reim];
      }
    }
  }

  // Fill the gauge field with noise 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) {
      for(int reim=0; reim < 2; reim++) { 
	gauge_field[col1][col2][reim] =  (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
      }
    }
  }

  /* Stuff goes here */
  _sse_pair_load(input_spinor[0], input_spinor[1]);
  _sse_pair_load_up(input_spinor[2],input_spinor[3]);
  _sse_42_gamma0_minus();
  _sse_su3_inverse_multiply(gauge_field);
  _sse_vector_load(upper_sum1);
  _sse_vector_add();
  _sse_vector_store(upper_sum1);
  _sse_vector_load(lower_sum1);
  _sse_24_gamma0_minus_add();
  _sse_vector_store(lower_sum1);  



  dslash_minus_dir0_backward_add(input_spinor, gauge_field, upper_sum2, lower_sum2);

  // Check on the results
  for(int spin2=0; spin2 < 2; spin2++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 

	float diff = upper_sum2[col][spin2][reim] - upper_sum1[col][spin2][reim];
	diff /= (float)(2*3*2);

#if 0
	QDPIO::cout << "\t sp=" << spin2 
		    << " co=" << col 
		    << " re=" << reim 
		    << " diff= " << diff << std::endl;
#endif
	assertion( fabs(diff) < 1.0e-9 );

	diff = lower_sum2[col][spin2][reim] - lower_sum1[col][spin2][reim];
	diff /= (float)(2*3*2);
#if 0
	QDPIO::cout << "\t sp=" << spin2 
		    << " co=" << col 
		    << " re=" << reim 
		    << " diff2 = " << diff << std::endl;
#endif
	assertion( fabs(diff) < 1.0e-9 );

      }
    }
  }


}


void 
testSiteDslash1MinusForwardAdd::run()
{
  spinor_array input_spinor ALIGN;
  u_mat_array gauge_field ALIGN;
  halfspinor_array upper_sum1 ALIGN, upper_sum2 ALIGN;
  halfspinor_array lower_sum1 ALIGN, lower_sum2 ALIGN;

  // Fill inpur spinor with noise 
  for(int spin4=0; spin4 < 4; spin4++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	input_spinor[spin4][col][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
      }
    }
  }

  for(int spin2=0; spin2 < 2; spin2++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	
	upper_sum1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	lower_sum1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	upper_sum2[col][spin2][reim] = upper_sum1[col][spin2][reim];
	lower_sum2[col][spin2][reim] = lower_sum1[col][spin2][reim];
      }
    }
  }

  // Fill the gauge field with noise 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) {
      for(int reim=0; reim < 2; reim++) { 
	gauge_field[col1][col2][reim] =  (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
      }
    }
  }

  /* Stuff goes here */
  _sse_pair_load(input_spinor[0], input_spinor[1]);
  _sse_pair_load_up(input_spinor[2],input_spinor[3]);
  _sse_42_gamma1_plus();
  _sse_su3_multiply(gauge_field);
  _sse_vector_load(upper_sum1);
  _sse_vector_add();
  _sse_vector_store(upper_sum1);
  _sse_vector_load(lower_sum1);
  _sse_24_gamma1_plus();
  _sse_vector_store(lower_sum1);  



  dslash_minus_dir1_forward_add(input_spinor, gauge_field, upper_sum2, lower_sum2);

  // Check on the results
  for(int spin2=0; spin2 < 2; spin2++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 

	float diff = upper_sum2[col][spin2][reim] - upper_sum1[col][spin2][reim];
	diff /= (float)(2*3*2);

#if 0
	QDPIO::cout << "\t sp=" << spin2 
		    << " co=" << col 
		    << " re=" << reim 
		    << " diff= " << diff << std::endl;
#endif
	assertion( fabs(diff) < 1.0e-9 );

	diff = lower_sum2[col][spin2][reim] - lower_sum1[col][spin2][reim];
	diff /= (float)(2*3*2);
#if 0
	QDPIO::cout << "\t sp=" << spin2 
		    << " co=" << col 
		    << " re=" << reim 
		    << " diff2 = " << diff << std::endl;
#endif
	assertion( fabs(diff) < 1.0e-9 );

      }
    }
  }
 
}

void 
testSiteDslash1MinusBackwardAdd::run()
{
  spinor_array input_spinor ALIGN;
  u_mat_array gauge_field ALIGN;
  halfspinor_array upper_sum1 ALIGN, upper_sum2 ALIGN;
  halfspinor_array lower_sum1 ALIGN, lower_sum2 ALIGN;

  // Fill inpur spinor with noise 
  for(int spin4=0; spin4 < 4; spin4++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	input_spinor[spin4][col][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
      }
    }
  }

  for(int spin2=0; spin2 < 2; spin2++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	
	upper_sum1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	lower_sum1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	upper_sum2[col][spin2][reim] = upper_sum1[col][spin2][reim];
	lower_sum2[col][spin2][reim] = lower_sum1[col][spin2][reim];
      }
    }
  }

  // Fill the gauge field with noise 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) {
      for(int reim=0; reim < 2; reim++) { 
	gauge_field[col1][col2][reim] =  (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
      }
    }
  }

  /* Stuff goes here */
  _sse_pair_load(input_spinor[0], input_spinor[1]);
  _sse_pair_load_up(input_spinor[2],input_spinor[3]);
  _sse_42_gamma1_minus();
  _sse_su3_inverse_multiply(gauge_field);
  _sse_vector_load(upper_sum1);
  _sse_vector_add();
  _sse_vector_store(upper_sum1);
  _sse_vector_load(lower_sum1);
  _sse_24_gamma1_minus();
  _sse_vector_store(lower_sum1);  



  dslash_minus_dir1_backward_add(input_spinor, gauge_field, upper_sum2, lower_sum2);

  // Check on the results
  for(int spin2=0; spin2 < 2; spin2++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 

	float diff = upper_sum2[col][spin2][reim] - upper_sum1[col][spin2][reim];
	diff /= (float)(2*3*2);

#if 0
	QDPIO::cout << "\t sp=" << spin2 
		    << " co=" << col 
		    << " re=" << reim 
		    << " diff= " << diff << std::endl;
#endif
	assertion( fabs(diff) < 1.0e-9 );

	diff = lower_sum2[col][spin2][reim] - lower_sum1[col][spin2][reim];
	diff /= (float)(2*3*2);
#if 0
	QDPIO::cout << "\t sp=" << spin2 
		    << " co=" << col 
		    << " re=" << reim 
		    << " diff2 = " << diff << std::endl;
#endif
	assertion( fabs(diff) < 1.0e-9 );

      }
    }
  }
}



void 
testSiteDslash2MinusForwardAdd::run()
{
  spinor_array input_spinor ALIGN;
  u_mat_array gauge_field ALIGN;
  halfspinor_array upper_sum1 ALIGN, upper_sum2 ALIGN;
  halfspinor_array lower_sum1 ALIGN, lower_sum2 ALIGN;

  // Fill inpur spinor with noise 
  for(int spin4=0; spin4 < 4; spin4++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	input_spinor[spin4][col][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
      }
    }
  }

  for(int spin2=0; spin2 < 2; spin2++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	
	upper_sum1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	lower_sum1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	upper_sum2[col][spin2][reim] = upper_sum1[col][spin2][reim];
	lower_sum2[col][spin2][reim] = lower_sum1[col][spin2][reim];
      }
    }
  }

  // Fill the gauge field with noise 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) {
      for(int reim=0; reim < 2; reim++) { 
	gauge_field[col1][col2][reim] =  (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
      }
    }
  }

  /* Stuff goes here */
  _sse_pair_load(input_spinor[0], input_spinor[1]);
  _sse_pair_load_up(input_spinor[2],input_spinor[3]);
  _sse_42_gamma2_plus();
  _sse_su3_multiply(gauge_field);
  _sse_vector_load(upper_sum1);
  _sse_vector_add();
  _sse_vector_store(upper_sum1);
  _sse_vector_load(lower_sum1);
  _sse_24_gamma2_plus();
  _sse_vector_store(lower_sum1);  



  dslash_minus_dir2_forward_add(input_spinor, gauge_field, upper_sum2, lower_sum2);

  // Check on the results
  for(int spin2=0; spin2 < 2; spin2++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 

	float diff = upper_sum2[col][spin2][reim] - upper_sum1[col][spin2][reim];
	diff /= (float)(2*3*2);

#if 0
	QDPIO::cout << "\t sp=" << spin2 
		    << " co=" << col 
		    << " re=" << reim 
		    << " diff= " << diff << std::endl;
#endif
	assertion( fabs(diff) < 1.0e-9 );

	diff = lower_sum2[col][spin2][reim] - lower_sum1[col][spin2][reim];
	diff /= (float)(2*3*2);
#if 0
	QDPIO::cout << "\t sp=" << spin2 
		    << " co=" << col 
		    << " re=" << reim 
		    << " diff2 = " << diff << std::endl;
#endif
	assertion( fabs(diff) < 1.0e-9 );

      }
    }
  }
 
}

void 
testSiteDslash2MinusBackwardAdd::run()
{
  spinor_array input_spinor ALIGN;
  u_mat_array gauge_field ALIGN;
  halfspinor_array upper_sum1 ALIGN, upper_sum2 ALIGN;
  halfspinor_array lower_sum1 ALIGN, lower_sum2 ALIGN;

  // Fill inpur spinor with noise 
  for(int spin4=0; spin4 < 4; spin4++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	input_spinor[spin4][col][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
      }
    }
  }

  for(int spin2=0; spin2 < 2; spin2++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	
	upper_sum1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	lower_sum1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	upper_sum2[col][spin2][reim] = upper_sum1[col][spin2][reim];
	lower_sum2[col][spin2][reim] = lower_sum1[col][spin2][reim];
      }
    }
  }

  // Fill the gauge field with noise 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) {
      for(int reim=0; reim < 2; reim++) { 
	gauge_field[col1][col2][reim] =  (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
      }
    }
  }

  /* Stuff goes here */
  _sse_pair_load(input_spinor[0], input_spinor[1]);
  _sse_pair_load_up(input_spinor[2],input_spinor[3]);
  _sse_42_gamma2_minus();
  _sse_su3_inverse_multiply(gauge_field);
  _sse_vector_load(upper_sum1);
  _sse_vector_add();
  _sse_vector_store(upper_sum1);
  _sse_vector_load(lower_sum1);
  _sse_24_gamma2_minus();
  _sse_vector_store(lower_sum1);  



  dslash_minus_dir2_backward_add(input_spinor, gauge_field, upper_sum2, lower_sum2);

  // Check on the results
  for(int spin2=0; spin2 < 2; spin2++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 

	float diff = upper_sum2[col][spin2][reim] - upper_sum1[col][spin2][reim];
	diff /= (float)(2*3*2);

#if 0
	QDPIO::cout << "\t sp=" << spin2 
		    << " co=" << col 
		    << " re=" << reim 
		    << " diff= " << diff << std::endl;
#endif
	assertion( fabs(diff) < 1.0e-9 );

	diff = lower_sum2[col][spin2][reim] - lower_sum1[col][spin2][reim];
	diff /= (float)(2*3*2);
#if 0
	QDPIO::cout << "\t sp=" << spin2 
		    << " co=" << col 
		    << " re=" << reim 
		    << " diff2 = " << diff << std::endl;
#endif
	assertion( fabs(diff) < 1.0e-9 );

      }
    }
  }
}


void 
testSiteDslash3MinusForwardAdd::run()
{

  spinor_array input_spinor ALIGN;
  u_mat_array gauge_field ALIGN;
  halfspinor_array upper_sum1 ALIGN, upper_sum2 ALIGN;
  halfspinor_array lower_sum1 ALIGN, lower_sum2 ALIGN;

  // Fill inpur spinor with noise 
  for(int spin4=0; spin4 < 4; spin4++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	input_spinor[spin4][col][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
      }
    }
  }

  for(int spin2=0; spin2 < 2; spin2++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	
	upper_sum1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	lower_sum1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	upper_sum2[col][spin2][reim] = upper_sum1[col][spin2][reim];
	lower_sum2[col][spin2][reim] = lower_sum1[col][spin2][reim];
      }
    }
  }

  // Fill the gauge field with noise 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) {
      for(int reim=0; reim < 2; reim++) { 
	gauge_field[col1][col2][reim] =  (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
      }
    }
  }

  /* Stuff goes here */
  _sse_pair_load(input_spinor[0], input_spinor[1]);
  _sse_pair_load_up(input_spinor[2],input_spinor[3]);
  _sse_42_gamma3_plus();
  _sse_su3_multiply(gauge_field);
  _sse_vector_load(upper_sum1);
  _sse_vector_add();
  _sse_vector_store(upper_sum1);
  _sse_vector_load(lower_sum1);
  _sse_24_gamma3_plus();
  _sse_vector_store(lower_sum1);  



  dslash_minus_dir3_forward_add(input_spinor, gauge_field, upper_sum2, lower_sum2);

  // Check on the results
  for(int spin2=0; spin2 < 2; spin2++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 

	float diff = upper_sum2[col][spin2][reim] - upper_sum1[col][spin2][reim];
	diff /= (float)(2*3*2);

#if 0
	QDPIO::cout << "\t sp=" << spin2 
		    << " co=" << col 
		    << " re=" << reim 
		    << " diff= " << diff << std::endl;
#endif
	assertion( fabs(diff) < 1.0e-9 );

	diff = lower_sum2[col][spin2][reim] - lower_sum1[col][spin2][reim];
	diff /= (float)(2*3*2);
#if 0
	QDPIO::cout << "\t sp=" << spin2 
		    << " co=" << col 
		    << " re=" << reim 
		    << " diff2 = " << diff << std::endl;
#endif
	assertion( fabs(diff) < 1.0e-9 );

      }
    }
  }

}



void 
testSiteDslash3MinusBackwardAddStore::run()
{

  spinor_array input_spinor ALIGN;
  spinor_array output_spinor ALIGN, output_spinor2;
  u_mat_array gauge_field ALIGN;
  halfspinor_array upper_sum1 ALIGN, upper_sum2 ALIGN;
  halfspinor_array lower_sum1 ALIGN, lower_sum2 ALIGN;
  
  // Fill inpur spinor with noise 
  for(int spin4=0; spin4 < 4; spin4++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	input_spinor[spin4][col][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
      }
    }
  }

  for(int spin2=0; spin2 < 2; spin2++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	
	upper_sum1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	lower_sum1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	upper_sum2[col][spin2][reim] = upper_sum1[col][spin2][reim];
	lower_sum2[col][spin2][reim] = lower_sum1[col][spin2][reim];
      }
    }
  }

  // Fill the gauge field with noise 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) {
      for(int reim=0; reim < 2; reim++) { 
	gauge_field[col1][col2][reim] =  (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
      }
    }
  }

  /* Stuff goes here */
  _sse_pair_load(input_spinor[0], input_spinor[1]);
  _sse_pair_load_up(input_spinor[2],input_spinor[3]);
  _sse_42_gamma3_minus();
  _sse_su3_inverse_multiply(gauge_field);
  _sse_vector_load(upper_sum1);
  _sse_24_gamma3_minus_rows12();
  _sse_pair_store(output_spinor[0], output_spinor[1]);
  _sse_vector_load(lower_sum1);
  _sse_24_gamma3_minus();
  _sse_pair_store(output_spinor[2], output_spinor[3]);



  dslash_minus_dir3_backward_add_store(input_spinor, gauge_field, upper_sum2, lower_sum2, output_spinor2);

  // Check on the results
  for(int spin4=0; spin4 < 4; spin4++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 

	float diff = output_spinor[spin4][col][reim] - output_spinor2[spin4][col][reim];
	diff /= (float)(4*3*2);

#if 0
	QDPIO::cout << "\t sp=" << spin4
		    << " co=" << col 
		    << " re=" << reim 
		    << " diff= " << diff << std::endl;
#endif
	assertion( fabs(diff) < 1.0e-9 );


      }
    }
  }
}




