#include "testSiteDslash.h"

#include "sse_align.h"
#include "types64.h"
#include "sse64.h"

#include <cmath>
#include "site_dslash_64bit_scalar.h"

using namespace QDP;
using namespace Assertions;

  /* gamma 0 */

  /* Project */
#define _sse_42_1_gamma0_minus(sp) \
      _sse_load((sp)[0]); \
      _sse_load_up((sp)[3]);\
      _sse_vector_i_mul();\
      _sse_vector_sub()

  /* Recons with Store */
#define _sse_24_1_gamma0_minus_set() \
      _sse_store_up(rs[0]);\
      _sse_vector_i_mul_up();\
      _sse_store_up(rs[3]) 

  /* Recons with Accumulate */
#define _sse_24_1_gamma0_minus_add() \
      _sse_load(rs[0]);\
      _sse_vector_add();\
      _sse_store(rs[0]);\
      _sse_load(rs[3]);\
      _sse_vector_i_mul();\
      _sse_vector_add();\
      _sse_store(rs[3]) 
	
  /* Project 2 */
#define _sse_42_2_gamma0_minus(sp) \
      _sse_load((sp)[1]);\
      _sse_load_up((sp)[2]);\
      _sse_vector_i_mul();\
      _sse_vector_sub()

  /* Reconst store */
#define _sse_24_2_gamma0_minus_set() \
	  _sse_store_up(rs[1]);\
      _sse_vector_i_mul_up();\
      _sse_store_up(rs[2])

  /* Recons add ... etc */
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


void 
testSiteDslash0PlusForward::run()
{

  spinor_array input_spinor ALIGN;
  u_mat_array gauge_field ALIGN;
  spinor_array rs ALIGN;
  spinor_array output_spinor2 ALIGN;

  // Fill inpur spinor with noise 
  for(int spin4=0; spin4 < 4; spin4++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	input_spinor[spin4][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }
    }
  }

  for(int spin=0; spin < 4; spin++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	
	rs[spin][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
	output_spinor2[spin][col][reim] = rs[spin][col][reim];
      }
    }
  }

  // Fill the gauge field with noise 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) {
      for(int reim=0; reim < 2; reim++) { 
	gauge_field[col1][col2][reim] =  (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }
    }
  }


 _sse_42_1_gamma0_minus(input_spinor);
 _sse_su3_multiply(gauge_field);
 _sse_24_1_gamma0_minus_set();
 _sse_42_2_gamma0_minus(input_spinor);
 _sse_su3_multiply(gauge_field);
 _sse_24_2_gamma0_minus_set();

  dslash_plus_dir0_forward(input_spinor, gauge_field, output_spinor2);


  // Check on the results
  for(int spin=0; spin < 4; spin++) {

    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 

	double diff = rs[spin][col][reim] - output_spinor2[spin][col][reim];
	diff /= (double)(4*3*2);

#if 0
	QDPIO::cout << "\t sp=" << spin
		    << " co=" << col 
		    << " re=" << reim 
		    << " diff= " << diff << std::endl;
#endif
	assertion( fabs(diff) < 1.0e-17 );

      }
    }
  }


}

void 
testSiteDslash0PlusBackwardAdd::run()
{
  spinor_array input_spinor ALIGN;
  u_mat_array gauge_field ALIGN;
  spinor_array rs ALIGN;
  spinor_array output_spinor2 ALIGN;

  // Fill inpur spinor with noise 
  for(int spin4=0; spin4 < 4; spin4++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	input_spinor[spin4][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }
    }
  }

  for(int spin=0; spin < 4; spin++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	
	rs[spin][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
	output_spinor2[spin][col][reim] = rs[spin][col][reim];
      }
    }
  }

  // Fill the gauge field with noise 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) {
      for(int reim=0; reim < 2; reim++) { 
	gauge_field[col1][col2][reim] =  (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }
    }
  }


 _sse_42_1_gamma0_plus(input_spinor);
 _sse_su3_inverse_multiply(gauge_field);
 _sse_24_1_gamma0_plus_add();
 _sse_42_2_gamma0_plus(input_spinor);
 _sse_su3_inverse_multiply(gauge_field);
 _sse_24_2_gamma0_plus_add();

  dslash_plus_dir0_backward_add(input_spinor, gauge_field, output_spinor2);


  // Check on the results
  for(int spin=0; spin < 4; spin++) {

    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 

	double diff = rs[spin][col][reim] - output_spinor2[spin][col][reim];
	diff /= (double)(4*3*2);

#if 0
	QDPIO::cout << "\t sp=" << spin
		    << " co=" << col 
		    << " re=" << reim 
		    << " diff= " << diff << std::endl;
#endif
	assertion( fabs(diff) < 1.0e-17 );

      }
    }
  }


}

void 
testSiteDslash1PlusForwardAdd::run()
{
  
  spinor_array input_spinor ALIGN;
  u_mat_array gauge_field ALIGN;
  spinor_array rs ALIGN;
  spinor_array output_spinor2 ALIGN;

  // Fill inpur spinor with noise 
  for(int spin4=0; spin4 < 4; spin4++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	input_spinor[spin4][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }
    }
  }

  for(int spin=0; spin < 4; spin++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	
	rs[spin][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
	output_spinor2[spin][col][reim] = rs[spin][col][reim];
      }
    }
  }

  // Fill the gauge field with noise 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) {
      for(int reim=0; reim < 2; reim++) { 
	gauge_field[col1][col2][reim] =  (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }
    }
  }


 _sse_42_1_gamma1_minus(input_spinor);
 _sse_su3_multiply(gauge_field);
 _sse_24_1_gamma1_minus();
 _sse_42_2_gamma1_minus(input_spinor);
 _sse_su3_multiply(gauge_field);
 _sse_24_2_gamma1_minus();

  dslash_plus_dir1_forward_add(input_spinor, gauge_field, output_spinor2);


  // Check on the results
  for(int spin=0; spin < 4; spin++) {

    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 

	double diff = rs[spin][col][reim] - output_spinor2[spin][col][reim];
	diff /= (double)(4*3*2);

#if 0
	QDPIO::cout << "\t sp=" << spin
		    << " co=" << col 
		    << " re=" << reim 
		    << " diff= " << diff << std::endl;
#endif
	assertion( fabs(diff) < 1.0e-17 );

      }
    }
  }


}

void 
testSiteDslash1PlusBackwardAdd::run()
{

  spinor_array input_spinor ALIGN;
  u_mat_array gauge_field ALIGN;
  spinor_array rs ALIGN;
  spinor_array output_spinor2 ALIGN;

  // Fill inpur spinor with noise 
  for(int spin4=0; spin4 < 4; spin4++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	input_spinor[spin4][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }
    }
  }

  for(int spin=0; spin < 4; spin++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	
	rs[spin][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
	output_spinor2[spin][col][reim] = rs[spin][col][reim];
      }
    }
  }

  // Fill the gauge field with noise 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) {
      for(int reim=0; reim < 2; reim++) { 
	gauge_field[col1][col2][reim] =  (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }
    }
  }


 _sse_42_1_gamma1_plus(input_spinor);
 _sse_su3_inverse_multiply(gauge_field);
 _sse_24_1_gamma1_plus();
 _sse_42_2_gamma1_plus(input_spinor);
 _sse_su3_inverse_multiply(gauge_field);
 _sse_24_2_gamma1_plus();

  dslash_plus_dir1_backward_add(input_spinor, gauge_field, output_spinor2);


  // Check on the results
  for(int spin=0; spin < 4; spin++) {

    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 

	double diff = rs[spin][col][reim] - output_spinor2[spin][col][reim];
	diff /= (double)(4*3*2);

#if 0
	QDPIO::cout << "\t sp=" << spin
		    << " co=" << col 
		    << " re=" << reim 
		    << " diff= " << diff << std::endl;
#endif
	assertion( fabs(diff) < 1.0e-17 );

      }
    }
  }


}


void 
testSiteDslash2PlusForwardAdd::run()
{

  spinor_array input_spinor ALIGN;
  u_mat_array gauge_field ALIGN;
  spinor_array rs ALIGN;
  spinor_array output_spinor2 ALIGN;

  // Fill inpur spinor with noise 
  for(int spin4=0; spin4 < 4; spin4++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	input_spinor[spin4][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }
    }
  }

  for(int spin=0; spin < 4; spin++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	
	rs[spin][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
	output_spinor2[spin][col][reim] = rs[spin][col][reim];
      }
    }
  }

  // Fill the gauge field with noise 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) {
      for(int reim=0; reim < 2; reim++) { 
	gauge_field[col1][col2][reim] =  (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }
    }
  }


 _sse_42_1_gamma2_minus(input_spinor);
 _sse_su3_multiply(gauge_field);
 _sse_24_1_gamma2_minus();
 _sse_42_2_gamma2_minus(input_spinor);
 _sse_su3_multiply(gauge_field);
 _sse_24_2_gamma2_minus();

  dslash_plus_dir2_forward_add(input_spinor, gauge_field, output_spinor2);


  // Check on the results
  for(int spin=0; spin < 4; spin++) {

    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 

	double diff = rs[spin][col][reim] - output_spinor2[spin][col][reim];
	diff /= (double)(4*3*2);

#if 0
	QDPIO::cout << "\t sp=" << spin
		    << " co=" << col 
		    << " re=" << reim 
		    << " diff= " << diff << std::endl;
#endif
	assertion( fabs(diff) < 1.0e-17 );

      }
    }
  }


}

void 
testSiteDslash2PlusBackwardAdd::run()
{

  spinor_array input_spinor ALIGN;
  u_mat_array gauge_field ALIGN;
  spinor_array rs ALIGN;
  spinor_array output_spinor2 ALIGN;

  // Fill inpur spinor with noise 
  for(int spin4=0; spin4 < 4; spin4++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	input_spinor[spin4][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }
    }
  }

  for(int spin=0; spin < 4; spin++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	
	rs[spin][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
	output_spinor2[spin][col][reim] = rs[spin][col][reim];
      }
    }
  }

  // Fill the gauge field with noise 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) {
      for(int reim=0; reim < 2; reim++) { 
	gauge_field[col1][col2][reim] =  (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }
    }
  }


 _sse_42_1_gamma2_plus(input_spinor);
 _sse_su3_inverse_multiply(gauge_field);
 _sse_24_1_gamma2_plus();
 _sse_42_2_gamma2_plus(input_spinor);
 _sse_su3_inverse_multiply(gauge_field);
 _sse_24_2_gamma2_plus();

  dslash_plus_dir2_backward_add(input_spinor, gauge_field, output_spinor2);


  // Check on the results
  for(int spin=0; spin < 4; spin++) {

    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 

	double diff = rs[spin][col][reim] - output_spinor2[spin][col][reim];
	diff /= (double)(4*3*2);

#if 0
	QDPIO::cout << "\t sp=" << spin
		    << " co=" << col 
		    << " re=" << reim 
		    << " diff= " << diff << std::endl;
#endif
	assertion( fabs(diff) < 1.0e-17 );

      }
    }
  }


}


void 
testSiteDslash3PlusForwardAdd::run()
{

  spinor_array input_spinor ALIGN;
  u_mat_array gauge_field ALIGN;
  spinor_array rs ALIGN;
  spinor_array output_spinor2 ALIGN;

  // Fill inpur spinor with noise 
  for(int spin4=0; spin4 < 4; spin4++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	input_spinor[spin4][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }
    }
  }

  for(int spin=0; spin < 4; spin++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	
	rs[spin][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
	output_spinor2[spin][col][reim] = rs[spin][col][reim];
      }
    }
  }

  // Fill the gauge field with noise 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) {
      for(int reim=0; reim < 2; reim++) { 
	gauge_field[col1][col2][reim] =  (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }
    }
  }


 _sse_42_1_gamma3_minus(input_spinor);
 _sse_su3_multiply(gauge_field);
 _sse_24_1_gamma3_minus_add();
 _sse_42_2_gamma3_minus(input_spinor);
 _sse_su3_multiply(gauge_field);
 _sse_24_2_gamma3_minus_add();

  dslash_plus_dir3_forward_add(input_spinor, gauge_field, output_spinor2);


  // Check on the results
  for(int spin=0; spin < 4; spin++) {

    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 

	double diff = rs[spin][col][reim] - output_spinor2[spin][col][reim];
	diff /= (double)(4*3*2);

#if 0
	QDPIO::cout << "\t sp=" << spin
		    << " co=" << col 
		    << " re=" << reim 
		    << " diff= " << diff << std::endl;
#endif
	assertion( fabs(diff) < 1.0e-17 );

      }
    }
  }


}



void 
testSiteDslash3PlusBackwardAddStore::run()
{
  spinor_array input_spinor ALIGN;
  u_mat_array gauge_field ALIGN;
  spinor_array rs ALIGN;
  spinor_array output_spinor ALIGN;
  spinor_array output_spinor2 ALIGN;
  spinor_array* rn = &output_spinor;

  // Fill inpur spinor with noise 
  for(int spin4=0; spin4 < 4; spin4++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	input_spinor[spin4][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }
    }
  }

  for(int spin=0; spin < 4; spin++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	
	rs[spin][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
	output_spinor2[spin][col][reim] = rs[spin][col][reim];
      }
    }
  }

  // Fill the gauge field with noise 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) {
      for(int reim=0; reim < 2; reim++) { 
	gauge_field[col1][col2][reim] =  (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }
    }
  }


 _sse_42_1_gamma3_plus(input_spinor);
 _sse_su3_inverse_multiply(gauge_field);
 _sse_24_1_gamma3_plus_set();
 _sse_42_2_gamma3_plus(input_spinor);
 _sse_su3_inverse_multiply(gauge_field);
 _sse_24_2_gamma3_plus_set();

  dslash_plus_dir3_backward_add_store(input_spinor, gauge_field, output_spinor2);


  // Check on the results
  for(int spin=0; spin < 4; spin++) {

    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 

	double diff = output_spinor[spin][col][reim] - output_spinor2[spin][col][reim];
	diff /= (double)(4*3*2);

#if 0
	QDPIO::cout << "\t sp=" << spin
		    << " co=" << col 
		    << " re=" << reim 
		    << " diff= " << diff << std::endl;
#endif
	assertion( fabs(diff) < 1.0e-17 );

      }
    }
  }


}


void 
testSiteDslash0MinusForward::run()
{
  spinor_array input_spinor ALIGN;
  u_mat_array gauge_field ALIGN;
  spinor_array rs ALIGN;
  spinor_array output_spinor2 ALIGN;

  // Fill inpur spinor with noise 
  for(int spin4=0; spin4 < 4; spin4++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	input_spinor[spin4][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }
    }
  }

  for(int spin=0; spin < 4; spin++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	
	rs[spin][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
	output_spinor2[spin][col][reim] = rs[spin][col][reim];
      }
    }
  }

  // Fill the gauge field with noise 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) {
      for(int reim=0; reim < 2; reim++) { 
	gauge_field[col1][col2][reim] =  (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }
    }
  }


 _sse_42_1_gamma0_plus(input_spinor);
 _sse_su3_multiply(gauge_field);
 _sse_24_1_gamma0_plus_set();
 _sse_42_2_gamma0_plus(input_spinor);
 _sse_su3_multiply(gauge_field);
 _sse_24_2_gamma0_plus_set();

  dslash_minus_dir0_forward(input_spinor, gauge_field, output_spinor2);


  // Check on the results
  for(int spin=0; spin < 4; spin++) {

    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 

	double diff = rs[spin][col][reim] - output_spinor2[spin][col][reim];
	diff /= (double)(4*3*2);

#if 0
	QDPIO::cout << "\t sp=" << spin
		    << " co=" << col 
		    << " re=" << reim 
		    << " diff= " << diff << std::endl;
#endif
	assertion( fabs(diff) < 1.0e-17 );

      }
    }
  }

}

void 
testSiteDslash0MinusBackwardAdd::run()
{
  spinor_array input_spinor ALIGN;
  u_mat_array gauge_field ALIGN;
  spinor_array rs ALIGN;
  spinor_array output_spinor2 ALIGN;

  // Fill inpur spinor with noise 
  for(int spin4=0; spin4 < 4; spin4++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	input_spinor[spin4][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }
    }
  }

  for(int spin=0; spin < 4; spin++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	
	rs[spin][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
	output_spinor2[spin][col][reim] = rs[spin][col][reim];
      }
    }
  }

  // Fill the gauge field with noise 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) {
      for(int reim=0; reim < 2; reim++) { 
	gauge_field[col1][col2][reim] =  (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }
    }
  }


 _sse_42_1_gamma0_minus(input_spinor);
 _sse_su3_inverse_multiply(gauge_field);
 _sse_24_1_gamma0_minus_add();
 _sse_42_2_gamma0_minus(input_spinor);
 _sse_su3_inverse_multiply(gauge_field);
 _sse_24_2_gamma0_minus_add();

  dslash_minus_dir0_backward_add(input_spinor, gauge_field, output_spinor2);


  // Check on the results
  for(int spin=0; spin < 4; spin++) {

    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 

	double diff = rs[spin][col][reim] - output_spinor2[spin][col][reim];
	diff /= (double)(4*3*2);

#if 0
	QDPIO::cout << "\t sp=" << spin
		    << " co=" << col 
		    << " re=" << reim 
		    << " diff= " << diff << std::endl;
#endif
	assertion( fabs(diff) < 1.0e-17 );

      }
    }
  }


}



void 
testSiteDslash1MinusForwardAdd::run()
{
  
  spinor_array input_spinor ALIGN;
  u_mat_array gauge_field ALIGN;
  spinor_array rs ALIGN;
  spinor_array output_spinor2 ALIGN;

  // Fill inpur spinor with noise 
  for(int spin4=0; spin4 < 4; spin4++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	input_spinor[spin4][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }
    }
  }

  for(int spin=0; spin < 4; spin++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	
	rs[spin][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
	output_spinor2[spin][col][reim] = rs[spin][col][reim];
      }
    }
  }

  // Fill the gauge field with noise 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) {
      for(int reim=0; reim < 2; reim++) { 
	gauge_field[col1][col2][reim] =  (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }
    }
  }


 _sse_42_1_gamma1_plus(input_spinor);
 _sse_su3_multiply(gauge_field);
 _sse_24_1_gamma1_plus();
 _sse_42_2_gamma1_plus(input_spinor);
 _sse_su3_multiply(gauge_field);
 _sse_24_2_gamma1_plus();

  dslash_minus_dir1_forward_add(input_spinor, gauge_field, output_spinor2);


  // Check on the results
  for(int spin=0; spin < 4; spin++) {

    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 

	double diff = rs[spin][col][reim] - output_spinor2[spin][col][reim];
	diff /= (double)(4*3*2);

#if 0
	QDPIO::cout << "\t sp=" << spin
		    << " co=" << col 
		    << " re=" << reim 
		    << " diff= " << diff << std::endl;
#endif
	assertion( fabs(diff) < 1.0e-17 );

      }
    }
  }



}

void 
testSiteDslash1MinusBackwardAdd::run()
{

  spinor_array input_spinor ALIGN;
  u_mat_array gauge_field ALIGN;
  spinor_array rs ALIGN;
  spinor_array output_spinor2 ALIGN;

  // Fill inpur spinor with noise 
  for(int spin4=0; spin4 < 4; spin4++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	input_spinor[spin4][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }
    }
  }

  for(int spin=0; spin < 4; spin++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	
	rs[spin][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
	output_spinor2[spin][col][reim] = rs[spin][col][reim];
      }
    }
  }

  // Fill the gauge field with noise 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) {
      for(int reim=0; reim < 2; reim++) { 
	gauge_field[col1][col2][reim] =  (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }
    }
  }


 _sse_42_1_gamma1_minus(input_spinor);
 _sse_su3_inverse_multiply(gauge_field);
 _sse_24_1_gamma1_minus();
 _sse_42_2_gamma1_minus(input_spinor);
 _sse_su3_inverse_multiply(gauge_field);
 _sse_24_2_gamma1_minus();

  dslash_minus_dir1_backward_add(input_spinor, gauge_field, output_spinor2);


  // Check on the results
  for(int spin=0; spin < 4; spin++) {

    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 

	double diff = rs[spin][col][reim] - output_spinor2[spin][col][reim];
	diff /= (double)(4*3*2);

#if 0
	QDPIO::cout << "\t sp=" << spin
		    << " co=" << col 
		    << " re=" << reim 
		    << " diff= " << diff << std::endl;
#endif
	assertion( fabs(diff) < 1.0e-17 );

      }
    }
  }


}



void 
testSiteDslash2MinusForwardAdd::run()
{

  spinor_array input_spinor ALIGN;
  u_mat_array gauge_field ALIGN;
  spinor_array rs ALIGN;
  spinor_array output_spinor2 ALIGN;

  // Fill inpur spinor with noise 
  for(int spin4=0; spin4 < 4; spin4++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	input_spinor[spin4][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }
    }
  }

  for(int spin=0; spin < 4; spin++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	
	rs[spin][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
	output_spinor2[spin][col][reim] = rs[spin][col][reim];
      }
    }
  }

  // Fill the gauge field with noise 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) {
      for(int reim=0; reim < 2; reim++) { 
	gauge_field[col1][col2][reim] =  (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }
    }
  }


 _sse_42_1_gamma2_plus(input_spinor);
 _sse_su3_multiply(gauge_field);
 _sse_24_1_gamma2_plus();
 _sse_42_2_gamma2_plus(input_spinor);
 _sse_su3_multiply(gauge_field);
 _sse_24_2_gamma2_plus();

  dslash_minus_dir2_forward_add(input_spinor, gauge_field, output_spinor2);


  // Check on the results
  for(int spin=0; spin < 4; spin++) {

    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 

	double diff = rs[spin][col][reim] - output_spinor2[spin][col][reim];
	diff /= (double)(4*3*2);

#if 0
	QDPIO::cout << "\t sp=" << spin
		    << " co=" << col 
		    << " re=" << reim 
		    << " diff= " << diff << std::endl;
#endif
	assertion( fabs(diff) < 1.0e-17 );

      }
    }
  }


}

void 
testSiteDslash2MinusBackwardAdd::run()
{

  spinor_array input_spinor ALIGN;
  u_mat_array gauge_field ALIGN;
  spinor_array rs ALIGN;
  spinor_array output_spinor2 ALIGN;

  // Fill inpur spinor with noise 
  for(int spin4=0; spin4 < 4; spin4++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	input_spinor[spin4][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }
    }
  }

  for(int spin=0; spin < 4; spin++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	
	rs[spin][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
	output_spinor2[spin][col][reim] = rs[spin][col][reim];
      }
    }
  }

  // Fill the gauge field with noise 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) {
      for(int reim=0; reim < 2; reim++) { 
	gauge_field[col1][col2][reim] =  (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }
    }
  }


 _sse_42_1_gamma2_minus(input_spinor);
 _sse_su3_inverse_multiply(gauge_field);
 _sse_24_1_gamma2_minus();
 _sse_42_2_gamma2_minus(input_spinor);
 _sse_su3_inverse_multiply(gauge_field);
 _sse_24_2_gamma2_minus();

  dslash_minus_dir2_backward_add(input_spinor, gauge_field, output_spinor2);


  // Check on the results
  for(int spin=0; spin < 4; spin++) {

    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 

	double diff = rs[spin][col][reim] - output_spinor2[spin][col][reim];
	diff /= (double)(4*3*2);

#if 0
	QDPIO::cout << "\t sp=" << spin
		    << " co=" << col 
		    << " re=" << reim 
		    << " diff= " << diff << std::endl;
#endif
	assertion( fabs(diff) < 1.0e-17 );

      }
    }
  }


}



void 
testSiteDslash3MinusForwardAdd::run()
{

  spinor_array input_spinor ALIGN;
  u_mat_array gauge_field ALIGN;
  spinor_array rs ALIGN;
  spinor_array output_spinor2 ALIGN;

  // Fill inpur spinor with noise 
  for(int spin4=0; spin4 < 4; spin4++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	input_spinor[spin4][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }
    }
  }

  for(int spin=0; spin < 4; spin++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	
	rs[spin][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
	output_spinor2[spin][col][reim] = rs[spin][col][reim];
      }
    }
  }

  // Fill the gauge field with noise 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) {
      for(int reim=0; reim < 2; reim++) { 
	gauge_field[col1][col2][reim] =  (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }
    }
  }


 _sse_42_1_gamma3_plus(input_spinor);
 _sse_su3_multiply(gauge_field);
 _sse_24_1_gamma3_plus_add();
 _sse_42_2_gamma3_plus(input_spinor);
 _sse_su3_multiply(gauge_field);
 _sse_24_2_gamma3_plus_add();

  dslash_minus_dir3_forward_add(input_spinor, gauge_field, output_spinor2);


  // Check on the results
  for(int spin=0; spin < 4; spin++) {

    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 

	double diff = rs[spin][col][reim] - output_spinor2[spin][col][reim];
	diff /= (double)(4*3*2);

#if 0
	QDPIO::cout << "\t sp=" << spin
		    << " co=" << col 
		    << " re=" << reim 
		    << " diff= " << diff << std::endl;
#endif
	assertion( fabs(diff) < 1.0e-17 );

      }
    }
  }


}



void 
testSiteDslash3MinusBackwardAddStore::run()
{
  spinor_array input_spinor ALIGN;
  u_mat_array gauge_field ALIGN;
  spinor_array rs ALIGN;
  spinor_array output_spinor ALIGN;
  spinor_array output_spinor2 ALIGN;
  spinor_array* rn = &output_spinor;

  // Fill inpur spinor with noise 
  for(int spin4=0; spin4 < 4; spin4++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	input_spinor[spin4][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }
    }
  }

  for(int spin=0; spin < 4; spin++) {
    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 
	
	rs[spin][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
	output_spinor2[spin][col][reim] = rs[spin][col][reim];
      }
    }
  }

  // Fill the gauge field with noise 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) {
      for(int reim=0; reim < 2; reim++) { 
	gauge_field[col1][col2][reim] =  (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }
    }
  }


 _sse_42_1_gamma3_minus(input_spinor);
 _sse_su3_inverse_multiply(gauge_field);
 _sse_24_1_gamma3_minus_set();
 _sse_42_2_gamma3_minus(input_spinor);
 _sse_su3_inverse_multiply(gauge_field);
 _sse_24_2_gamma3_minus_set();

  dslash_minus_dir3_backward_add_store(input_spinor, gauge_field, output_spinor2);


  // Check on the results
  for(int spin=0; spin < 4; spin++) {

    for(int col=0; col < 3; col++) { 
      for(int reim=0; reim <2 ; reim++) { 

	double diff = output_spinor[spin][col][reim] - output_spinor2[spin][col][reim];
	diff /= (double)(4*3*2);

#if 0
	QDPIO::cout << "\t sp=" << spin
		    << " co=" << col 
		    << " re=" << reim 
		    << " diff= " << diff << std::endl;
#endif
	assertion( fabs(diff) < 1.0e-17 );

      }
    }
  }


}




