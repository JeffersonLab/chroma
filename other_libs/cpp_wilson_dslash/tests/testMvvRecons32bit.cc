#include "testMvvRecons.h"
#include "sse32.h"
#include "types32.h"
#include "mvv_recons_32bit.h"
#include "unittest.h"

using namespace QDP;
using namespace Assertions;

#define _sse_24_gamma0_minus_set() _sse_vector_xch_i_mul_up()
#define _sse_24_gamma0_plus_set()  _sse_vector_xch_i_mul_neg_up()
#define _sse_24_gamma0_minus_add() _sse_vector_xch_i_add()
#define _sse_24_gamma0_plus_add()  _sse_vector_xch_i_sub()

#define _sse_24_gamma1_minus()     _sse_vector_xch();\
				   _sse_vector_subadd()
#define _sse_24_gamma1_plus()      _sse_vector_xch();\
		 		   _sse_vector_addsub()

#define _sse_24_gamma2_minus()     _sse_vector_i_addsub()
#define _sse_24_gamma2_plus()      _sse_vector_i_subadd()

#define _sse_24_gamma3_minus()     _sse_vector_sub()
#define _sse_24_gamma3_plus()      _sse_vector_add()
#define _sse_24_gamma3_minus_rows12() _sse_vector_add()
#define _sse_24_gamma3_plus_rows12() _sse_vector_add()

void testMvvRecons0Plus::run(void) 
{
  halfspinor_array hspinor;
  u_mat_array matrix;

  halfspinor_array r12_1 ALIGN, r34_1 ALIGN;
  halfspinor_array r12_2 ALIGN,r34_2 ALIGN;

  /* Random numbers in halfspinors */
  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	hspinor[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;

      }
    }
  }
  

  // Random matrix 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) { 
      for(int reim=0; reim < 2; reim++) { 
	matrix[col1][col2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
      }
    }
  }

  _sse_vector_load(hspinor);
  _sse_su3_multiply((matrix));
  _sse_vector_store_up(r12_1);
  _sse_24_gamma0_minus_set(); 
  _sse_vector_store_up(r34_1);

  mvv_recons_gamma0_plus(hspinor, matrix, r12_2, r34_2);
  
  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	float diff = r12_1[col][spin2][reim] - r12_2[col][spin2][reim];
	diff /= (float)(3*2*2);

#if 0
	QDPIO::cout << "  col=" << col 
		    << " sp=" << spin2
		    << " re=" << reim 
		    << " diff upper = " << diff 
		    << std::endl;
#endif

	assertion( fabs(diff) < 1.0e-9 );

	float diff2 = r34_1[col][spin2][reim] - r34_2[col][spin2][reim];
	diff2 /= (float)(3*2*2);
#if 0
	QDPIO::cout << "  col=" << col 
		    << " sp=" << spin2
		    << " re=" << reim 
		    << " diff lower = " << diff2 
		    << std::endl;
#endif 
	assertion( fabs(diff2) < 1.0e-9 );

      }
    }
  }

}

void testMvvRecons1PlusAdd::run(void) 
{

  halfspinor_array hspinor;
  u_mat_array matrix;

  halfspinor_array r12_1 ALIGN, r34_1 ALIGN;
  halfspinor_array r12_2 ALIGN,r34_2 ALIGN;

  /* Random numbers in halfspinors */
  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	hspinor[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	r12_1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	r34_1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;

	/* Set up 'existing partial sum' */
	r12_2[col][spin2][reim] = r12_1[col][spin2][reim];
	r34_2[col][spin2][reim] = r34_1[col][spin2][reim];
	

      }
    }
  }
  

  // Random matrix 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) { 
      for(int reim=0; reim < 2; reim++) { 
	matrix[col1][col2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
      }
    }
  }

  _sse_vector_load(hspinor);
  _sse_su3_multiply((matrix));
  _sse_vector_load(r12_1);
  _sse_vector_add();
  _sse_vector_store(r12_1);
  _sse_vector_load(r34_1);
  _sse_24_gamma1_minus();
  _sse_vector_store(r34_1);




  mvv_recons_gamma1_plus_add(hspinor, matrix, r12_2, r34_2);
  

  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	float diff = r12_1[col][spin2][reim] - r12_2[col][spin2][reim];
	diff /= (float)(3*2*2);

#if 0
	QDPIO::cout << "  col=" << col 
		    << " sp=" << spin2
		    << " re=" << reim 
		    << " diff upper = " << diff 
		    << std::endl;
#endif

	assertion( fabs(diff) < 1.0e-9 );

	float diff2 = r34_1[col][spin2][reim] - r34_2[col][spin2][reim];
	diff2 /= (float)(3*2*2);

#if 0
	QDPIO::cout << "  col=" << col 
		    << " sp=" << spin2
		    << " re=" << reim 
		    << " diff lower = " << diff2 
		    << std::endl;
#endif
	assertion( fabs(diff2) < 1.0e-9 );

      }
    }
  }

}

void testMvvRecons2PlusAdd::run(void) 
{

  halfspinor_array hspinor;
  u_mat_array matrix;

  halfspinor_array r12_1 ALIGN, r34_1 ALIGN;
  halfspinor_array r12_2 ALIGN,r34_2 ALIGN;

  /* Random numbers in halfspinors */
  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	hspinor[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	r12_1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	r34_1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;

	/* Set up 'existing partial sum' */
	r12_2[col][spin2][reim] = r12_1[col][spin2][reim];
	r34_2[col][spin2][reim] = r34_1[col][spin2][reim];
	

      }
    }
  }
  

  // Random matrix 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) { 
      for(int reim=0; reim < 2; reim++) { 
	matrix[col1][col2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
      }
    }
  }

  _sse_vector_load(hspinor);
  _sse_su3_multiply((matrix));
  _sse_vector_load(r12_1);
  _sse_vector_add();
  _sse_vector_store(r12_1);
  _sse_vector_load(r34_1);
  _sse_24_gamma2_minus();
  _sse_vector_store(r34_1);




  mvv_recons_gamma2_plus_add(hspinor, matrix, r12_2, r34_2);
  

  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	float diff = r12_1[col][spin2][reim] - r12_2[col][spin2][reim];
	diff /= (float)(3*2*2);

#if 0
	QDPIO::cout << "  col=" << col 
		    << " sp=" << spin2
		    << " re=" << reim 
		    << " diff upper = " << diff 
		    << std::endl;
#endif

	assertion( fabs(diff) < 1.0e-9 );

	float diff2 = r34_1[col][spin2][reim] - r34_2[col][spin2][reim];
	diff2 /= (float)(3*2*2);

#if 0
	QDPIO::cout << "  col=" << col 
		    << " sp=" << spin2
		    << " re=" << reim 
		    << " diff lower = " << diff2 
		    << std::endl;
#endif
	assertion( fabs(diff2) < 1.0e-9 );

      }
    }
  }

}

void testMvvRecons2PlusAddStore::run(void) 
{

  halfspinor_array hspinor;
  u_mat_array matrix;

  halfspinor_array r12_1 ALIGN, r34_1 ALIGN;
  halfspinor_array r12_2 ALIGN,r34_2 ALIGN;

  halfspinor_array hs1 ALIGN;
  halfspinor_array hs2 ALIGN;

  spinor_array spinor1;
  spinor_array spinor2;


  /* Random numbers in halfspinors */
  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	hspinor[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	r12_1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	r34_1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;

	/* Set up 'existing partial sum' */
	r12_2[col][spin2][reim] = r12_1[col][spin2][reim];
	r34_2[col][spin2][reim] = r34_1[col][spin2][reim];
	

      }
    }
  }
  

  // Random matrix 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) { 
      for(int reim=0; reim < 2; reim++) { 
	matrix[col1][col2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
      }
    }
  }

  _sse_vector_load(hspinor);
  _sse_su3_multiply((matrix));
  _sse_vector_load(r12_1);
  _sse_vector_add();
  _sse_vector_store(hs1);
  _sse_vector_load(r34_1);
  _sse_24_gamma2_minus();
  _sse_vector_store(hs2);



  mvv_recons_gamma2_plus_add_store(hspinor, matrix, r12_2, r34_2, spinor2);
  
  for(int j=0; j < 3*2*2; j++) { 
    float diff = *((float *)hs1+j) - *((float *)spinor2[0]+j);
#if 0
    QDPIO::cout << " diff lower = " << diff << std::endl;
#endif
    assertion( fabs(diff) < 1.0e-9 );

    float diff2 = *((float *)hs2+j) - *((float *)spinor2[2]+j);
#if 0
    QDPIO::cout << " diff upper = " << diff2 << std::endl;
#endif
    assertion( fabs(diff2) < 1.0e-9 );

  }
}

void testMvvRecons3PlusAddStore::run(void) 
{

  halfspinor_array hspinor;
  u_mat_array matrix;

  halfspinor_array r12_1 ALIGN, r34_1 ALIGN;
  halfspinor_array r12_2 ALIGN,r34_2 ALIGN;

  halfspinor_array hs1 ALIGN;
  halfspinor_array hs2 ALIGN;

  spinor_array spinor1;
  spinor_array spinor2;


  /* Random numbers in halfspinors */
  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	hspinor[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	r12_1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	r34_1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;

	/* Set up 'existing partial sum' */
	r12_2[col][spin2][reim] = r12_1[col][spin2][reim];
	r34_2[col][spin2][reim] = r34_1[col][spin2][reim];
	

      }
    }
  }
  

  // Random matrix 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) { 
      for(int reim=0; reim < 2; reim++) { 
	matrix[col1][col2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
      }
    }
  }

  _sse_vector_load(hspinor);
  _sse_su3_multiply((matrix));
  _sse_vector_load(r12_1);
  _sse_vector_add();
  _sse_vector_store(hs1);
  _sse_vector_load(r34_1);
  _sse_24_gamma3_minus();
  _sse_vector_store(hs2);



  mvv_recons_gamma3_plus_add_store(hspinor, matrix, r12_2, r34_2, spinor2);
  
  for(int j=0; j < 3*2*2; j++) { 
    float diff = *((float *)hs1+j) - *((float *)spinor2[0]+j);
#if 0
    QDPIO::cout << " diff lower = " << diff << std::endl;
#endif
    assertion( fabs(diff) < 1.0e-9 );

    float diff2 = *((float *)hs2+j) - *((float *)spinor2[2]+j);
#if 0
    QDPIO::cout << " diff upper = " << diff2 << std::endl;
#endif
    assertion( fabs(diff2) < 1.0e-9 );

  }
}

void testMvvRecons0Minus::run(void) 
{
  halfspinor_array hspinor;
  u_mat_array matrix;

  halfspinor_array r12_1 ALIGN, r34_1 ALIGN;
  halfspinor_array r12_2 ALIGN,r34_2 ALIGN;

  /* Random numbers in halfspinors */
  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	hspinor[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;

      }
    }
  }
  

  // Random matrix 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) { 
      for(int reim=0; reim < 2; reim++) { 
	matrix[col1][col2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
      }
    }
  }

  _sse_vector_load(hspinor);
  _sse_su3_multiply((matrix));
  _sse_vector_store_up(r12_1);
  _sse_24_gamma0_plus_set(); 
  _sse_vector_store_up(r34_1);

  mvv_recons_gamma0_minus(hspinor, matrix, r12_2, r34_2);
  
  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	float diff = r12_1[col][spin2][reim] - r12_2[col][spin2][reim];
	diff /= (float)(3*2*2);

#if 0
	QDPIO::cout << "  col=" << col 
		    << " sp=" << spin2
		    << " re=" << reim 
		    << " diff upper = " << diff 
		    << std::endl;
#endif
	assertion( fabs(diff) < 1.0e-9 );

	float diff2 = r34_1[col][spin2][reim] - r34_2[col][spin2][reim];
	diff2 /= (float)(3*2*2);

#if 0
	QDPIO::cout << "  col=" << col 
		    << " sp=" << spin2
		    << " re=" << reim 
		    << " diff lower = " << diff2 
		    << std::endl;
#endif 
	assertion( fabs(diff2) < 1.0e-9 );

      }
    }
  }

}

void testMvvRecons1MinusAdd::run(void) 
{

  halfspinor_array hspinor;
  u_mat_array matrix;

  halfspinor_array r12_1 ALIGN, r34_1 ALIGN;
  halfspinor_array r12_2 ALIGN,r34_2 ALIGN;

  /* Random numbers in halfspinors */
  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	hspinor[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	r12_1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	r34_1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;

	/* Set up 'existing partial sum' */
	r12_2[col][spin2][reim] = r12_1[col][spin2][reim];
	r34_2[col][spin2][reim] = r34_1[col][spin2][reim];
	

      }
    }
  }
  

  // Random matrix 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) { 
      for(int reim=0; reim < 2; reim++) { 
	matrix[col1][col2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
      }
    }
  }

  _sse_vector_load(hspinor);
  _sse_su3_multiply((matrix));
  _sse_vector_load(r12_1);
  _sse_vector_add();
  _sse_vector_store(r12_1);
  _sse_vector_load(r34_1);
  _sse_24_gamma1_plus();
  _sse_vector_store(r34_1);




  mvv_recons_gamma1_minus_add(hspinor, matrix, r12_2, r34_2);
  

  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	float diff = r12_1[col][spin2][reim] - r12_2[col][spin2][reim];
	diff /= (float)(3*2*2);

#if 0
	QDPIO::cout << "  col=" << col 
		    << " sp=" << spin2
		    << " re=" << reim 
		    << " diff upper = " << diff 
		    << std::endl;
#endif

	assertion( fabs(diff) < 1.0e-9 );

	float diff2 = r34_1[col][spin2][reim] - r34_2[col][spin2][reim];
	diff2 /= (float)(3*2*2);

#if 0
	QDPIO::cout << "  col=" << col 
		    << " sp=" << spin2
		    << " re=" << reim 
		    << " diff lower = " << diff2 
		    << std::endl;
#endif
	assertion( fabs(diff2) < 1.0e-9 );

      }
    }
  }

}

void testMvvRecons2MinusAdd::run(void) 
{

  halfspinor_array hspinor;
  u_mat_array matrix;

  halfspinor_array r12_1 ALIGN, r34_1 ALIGN;
  halfspinor_array r12_2 ALIGN,r34_2 ALIGN;

  /* Random numbers in halfspinors */
  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	hspinor[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	r12_1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	r34_1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;

	/* Set up 'existing partial sum' */
	r12_2[col][spin2][reim] = r12_1[col][spin2][reim];
	r34_2[col][spin2][reim] = r34_1[col][spin2][reim];
	

      }
    }
  }
  

  // Random matrix 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) { 
      for(int reim=0; reim < 2; reim++) { 
	matrix[col1][col2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
      }
    }
  }

  _sse_vector_load(hspinor);
  _sse_su3_multiply((matrix));
  _sse_vector_load(r12_1);
  _sse_vector_add();
  _sse_vector_store(r12_1);
  _sse_vector_load(r34_1);
  _sse_24_gamma2_plus();
  _sse_vector_store(r34_1);




  mvv_recons_gamma2_minus_add(hspinor, matrix, r12_2, r34_2);
  

  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	float diff = r12_1[col][spin2][reim] - r12_2[col][spin2][reim];
	diff /= (float)(3*2*2);

#if 0
	QDPIO::cout << "  col=" << col 
		    << " sp=" << spin2
		    << " re=" << reim 
		    << " diff upper = " << diff 
		    << std::endl;
#endif

	assertion( fabs(diff) < 1.0e-9 );

	float diff2 = r34_1[col][spin2][reim] - r34_2[col][spin2][reim];
	diff2 /= (float)(3*2*2);

#if 0
	QDPIO::cout << "  col=" << col 
		    << " sp=" << spin2
		    << " re=" << reim 
		    << " diff lower = " << diff2 
		    << std::endl;
#endif
	assertion( fabs(diff2) < 1.0e-9 );

      }
    }
  }

}


void testMvvRecons2MinusAddStore::run(void) 
{

  halfspinor_array hspinor;
  u_mat_array matrix;

  halfspinor_array r12_1 ALIGN, r34_1 ALIGN;
  halfspinor_array r12_2 ALIGN,r34_2 ALIGN;

  halfspinor_array hs1 ALIGN;
  halfspinor_array hs2 ALIGN;

  spinor_array spinor2;


  /* Random numbers in halfspinors */
  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	hspinor[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	r12_1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	r34_1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;

	/* Set up 'existing partial sum' */
	r12_2[col][spin2][reim] = r12_1[col][spin2][reim];
	r34_2[col][spin2][reim] = r34_1[col][spin2][reim];
	

      }
    }
  }
  

  // Random matrix 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) { 
      for(int reim=0; reim < 2; reim++) { 
	matrix[col1][col2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
      }
    }
  }

  _sse_vector_load(hspinor);
  _sse_su3_multiply((matrix));
  _sse_vector_load(r12_1);
  _sse_vector_add();
  _sse_vector_store(hs1);  /* Store sum but unswizzled */
  _sse_vector_load(r34_1);
  _sse_24_gamma2_plus();
  _sse_vector_store(hs2);  /* Store sum but unswizzeld */


  mvv_recons_gamma2_minus_add_store(hspinor, matrix, r12_2, r34_2, spinor2);
    
  for(int j=0; j < 3*2*2; j++) { 
    float diff = *((float *)hs1+j) - *((float *)spinor2[0]+j);
#if 0
    QDPIO::cout << " diff lower = " << diff << std::endl;
#endif
    assertion( fabs(diff) < 1.0e-9 );

    float diff2 = *((float *)hs2+j) - *((float *)spinor2[2]+j);
#if 0
    QDPIO::cout << " diff upper = " << diff2 << std::endl;
#endif
    assertion( fabs(diff2) < 1.0e-9 );

  }

}

void testMvvRecons3MinusAddStore::run(void) 
{

  halfspinor_array hspinor;
  u_mat_array matrix;

  halfspinor_array r12_1 ALIGN, r34_1 ALIGN;
  halfspinor_array r12_2 ALIGN,r34_2 ALIGN;

  halfspinor_array hs1 ALIGN;
  halfspinor_array hs2 ALIGN;

  spinor_array spinor2;


  /* Random numbers in halfspinors */
  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	hspinor[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	r12_1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	r34_1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;

	/* Set up 'existing partial sum' */
	r12_2[col][spin2][reim] = r12_1[col][spin2][reim];
	r34_2[col][spin2][reim] = r34_1[col][spin2][reim];
	

      }
    }
  }
  

  // Random matrix 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) { 
      for(int reim=0; reim < 2; reim++) { 
	matrix[col1][col2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
      }
    }
  }

  _sse_vector_load(hspinor);
  _sse_su3_multiply((matrix));
  _sse_vector_load(r12_1);
  _sse_vector_add();
  _sse_vector_store(hs1);  /* Store sum but unswizzled */
  _sse_vector_load(r34_1);
  _sse_24_gamma3_plus();
  _sse_vector_store(hs2);  /* Store sum but unswizzeld */


  mvv_recons_gamma3_minus_add_store(hspinor, matrix, r12_2, r34_2, spinor2);
    
  for(int j=0; j < 3*2*2; j++) { 
    float diff = *((float *)hs1+j) - *((float *)spinor2[0]+j);
#if 0
    QDPIO::cout << " diff lower = " << diff << std::endl;
#endif
    assertion( fabs(diff) < 1.0e-9 );

    float diff2 = *((float *)hs2+j) - *((float *)spinor2[2]+j);
#if 0
    QDPIO::cout << " diff upper = " << diff2 << std::endl;
#endif
    assertion( fabs(diff2) < 1.0e-9 );

  }

}
