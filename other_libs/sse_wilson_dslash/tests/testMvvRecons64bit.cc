#include "testMvvRecons.h"
#include "sse64.h"
#include "types64.h"
#include "mvv_recons_64bit.h"
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

#define _sse_24_1_gamma2_minus_set() \
      _sse_load(rs[0]);\
      _sse_vector_add();\
      _sse_store((*rn)[0]);			\
      _sse_load(rs[2]);\
      _sse_vector_i_mul();   \
      _sse_vector_add();\
      _sse_store((*rn)[2])
	  
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

#define _sse_24_2_gamma2_minus_set() \
	   _sse_load(rs[1]);\
      _sse_vector_add();\
      _sse_store((*rn)[1]);			\
      _sse_load(rs[3]);\
      _sse_vector_i_mul();   \
      _sse_vector_sub();\
      _sse_store((*rn)[3])

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

#define _sse_24_1_gamma2_plus_set() \
      _sse_load(rs[0]);\
      _sse_vector_add();\
      _sse_store((*rn)[0]);			\
      _sse_load(rs[2]);\
      _sse_vector_i_mul();   \
      _sse_vector_sub();\
      _sse_store((*rn)[2]);

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

#define _sse_24_2_gamma2_plus_set() \
      _sse_load(rs[1]);\
      _sse_vector_add();\
      _sse_store((*rn)[1]);			\
      _sse_load(rs[3]);\
      _sse_vector_i_mul();     \
      _sse_vector_add();\
      _sse_store((*rn)[3])





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
      _sse_load(rs[1]);		     \
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



void testMvvRecons0Plus::run(void) 
{
  halfspinor_array hspinor ALIGN;
  u_mat_array matrix;
  spinor_array rs; // Has to be the goddamn result for the
                   // assembler

  spinor_array result2 ALIGN;


  /* Random numbers in halfspinors */
  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	hspinor[spin2][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;

      }
    }
  }

/* Random numbers in halfspinors */
  for(int col=0; col < 3; col++) { 
    for(int spin4=0; spin4 < 4; spin4++) { 
      for(int reim=0; reim < 2; reim++) { 
	rs[spin4][col][reim] =0;
	result2[spin4][col][reim] =0;
      }
    }
  }
  
  
  // Random matrix 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) { 
      for(int reim=0; reim < 2; reim++) { 
	matrix[col1][col2][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }
    }
  }

  halfspinor_array *s3=&hspinor;
  u_mat_array* up=&matrix;

  _sse_load((*s3)[0]);
  _sse_su3_multiply(*up);
  _sse_24_1_gamma0_minus_set();

   _sse_load(hspinor[1]);
   _sse_su3_multiply(*up);
   _sse_24_2_gamma0_minus_set();


  mvv_recons_gamma0_plus(hspinor, matrix, result2);
  
  for(int col=0; col < 3; col++) { 
      for(int reim=0; reim < 2; reim++) { 
	for(int spin4=0; spin4 < 4; spin4++) { 
	  double diff = rs[spin4][col][reim] - result2[spin4][col][reim];
	  diff /= (double)(4*3*2) ;
#if 0
	  QDPIO::cout << "  col=" << col 
		    << " sp=" << spin4
		    << " re=" << reim 
		    << " diff upper = " << diff 
		    << std::endl;
#endif
	  assertion( fabs(diff) < 1.0e-17 );
	}
      }
  }
}

void testMvvRecons1PlusAdd::run(void) 
{
  halfspinor_array hspinor ALIGN;
  u_mat_array matrix ALIGN;
  spinor_array rs ALIGN; // Has to be the goddamn result for the
                   // assembler

  spinor_array result2 ALIGN;


  /* Random numbers in halfspinors */
  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	hspinor[spin2][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;

      }
    }
  }

/* Random numbers in halfspinors */
  for(int col=0; col < 3; col++) { 
    for(int spin4=0; spin4 < 4; spin4++) { 
      for(int reim=0; reim < 2; reim++) { 
	rs[spin4][col][reim] =  (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
	result2[spin4][col][reim] = rs[spin4][col][reim];
      }
    }
  }
  
  
  // Random matrix 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) { 
      for(int reim=0; reim < 2; reim++) { 
	matrix[col1][col2][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }
    }
  }

  halfspinor_array *s3=&hspinor;
  u_mat_array* up=&matrix;

  _sse_load((*s3)[0]);
  _sse_su3_multiply(*up);      
  _sse_24_1_gamma1_minus();

  _sse_load((*s3)[1]);
  _sse_su3_multiply(*up);
  _sse_24_2_gamma1_minus();


  mvv_recons_gamma1_plus_add(hspinor, matrix, result2);
  
  for(int col=0; col < 3; col++) { 
      for(int reim=0; reim < 2; reim++) { 
	for(int spin4=0; spin4 < 4; spin4++) { 
	  double diff = rs[spin4][col][reim] - result2[spin4][col][reim];
	  diff /= (double)(4*3*2) ;
#if 0
	  QDPIO::cout << "  col=" << col 
		    << " sp=" << spin4
		    << " re=" << reim 
		    << " diff upper = " << diff 
		    << std::endl;
#endif
	  assertion( fabs(diff) < 1.0e-17 );
	}
      }
  }
}

void testMvvRecons2PlusAdd::run(void) 
{
  halfspinor_array hspinor ALIGN;
  u_mat_array matrix ALIGN;
  spinor_array rs ALIGN; // Has to be the goddamn result for the
                   // assembler

  spinor_array result2 ALIGN;


  /* Random numbers in halfspinors */
  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	hspinor[spin2][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;

      }
    }
  }

/* Random numbers in halfspinors */
  for(int col=0; col < 3; col++) { 
    for(int spin4=0; spin4 < 4; spin4++) { 
      for(int reim=0; reim < 2; reim++) { 
	rs[spin4][col][reim] =  (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
	result2[spin4][col][reim] = rs[spin4][col][reim];
      }
    }
  }
  
  
  // Random matrix 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) { 
      for(int reim=0; reim < 2; reim++) { 
	matrix[col1][col2][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }
    }
  }

  halfspinor_array *s3=&hspinor;
  u_mat_array* up=&matrix;

  _sse_load((*s3)[0]);
  _sse_su3_multiply(*up);      
  _sse_24_1_gamma2_minus();

  _sse_load((*s3)[1]);
  _sse_su3_multiply(*up);
  _sse_24_2_gamma2_minus();


  mvv_recons_gamma2_plus_add(hspinor, matrix, result2);
  
  for(int col=0; col < 3; col++) { 
      for(int reim=0; reim < 2; reim++) { 
	for(int spin4=0; spin4 < 4; spin4++) { 
	  double diff = rs[spin4][col][reim] - result2[spin4][col][reim];
	  diff /= (double)(4*3*2) ;
#if 0
	  QDPIO::cout << "  col=" << col 
		    << " sp=" << spin4
		    << " re=" << reim 
		    << " diff upper = " << diff 
		    << std::endl;
#endif
	  assertion( fabs(diff) < 1.0e-17 );
	}
      }
  }
}

void testMvvRecons2PlusAddStore::run(void) 
{
  halfspinor_array hspinor ALIGN;
  u_mat_array matrix ALIGN;
  spinor_array rs ALIGN; // Has to be the goddamn result for the
                   // assembler
  spinor_array result1 ALIGN;
  spinor_array *rn = &result1;
  
  spinor_array result2 ALIGN;


  /* Random numbers in halfspinors */
  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	hspinor[spin2][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;

      }
    }
  }

/* Random numbers in halfspinors */
  for(int col=0; col < 3; col++) { 
    for(int spin4=0; spin4 < 4; spin4++) { 
      for(int reim=0; reim < 2; reim++) { 
	rs[spin4][col][reim] =  (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }
    }
  }
  
  
  // Random matrix 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) { 
      for(int reim=0; reim < 2; reim++) { 
	matrix[col1][col2][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }
    }
  }

  halfspinor_array *s3=&hspinor;
  u_mat_array* up=&matrix;

  _sse_load((*s3)[0]);
  _sse_su3_multiply(*up);      
  _sse_24_1_gamma2_minus_set();


  _sse_load((*s3)[1]);
  _sse_su3_multiply(*up);
  _sse_24_2_gamma2_minus_set();


  mvv_recons_gamma2_plus_add_store(hspinor, matrix, rs, result2);
  
  for(int col=0; col < 3; col++) { 
      for(int reim=0; reim < 2; reim++) { 
	for(int spin4=0; spin4 < 4; spin4++) { 
	  double diff = result1[spin4][col][reim] - result2[spin4][col][reim];
	  diff /= (double)(4*3*2) ;
#if 0
	  QDPIO::cout << "  col=" << col 
		    << " sp=" << spin4
		    << " re=" << reim 
		    << " diff upper = " << diff 
		    << std::endl;
#endif
	  assertion( fabs(diff) < 1.0e-17 );
	}
      }
  }
}
void testMvvRecons3PlusAddStore::run(void) 
{
  halfspinor_array hspinor ALIGN;
  u_mat_array matrix ALIGN;
  spinor_array rs ALIGN; // Has to be the goddamn result for the
                   // assembler
  spinor_array result1 ALIGN;
  spinor_array *rn = &result1;
  
  spinor_array result2 ALIGN;


  /* Random numbers in halfspinors */
  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	hspinor[spin2][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;

      }
    }
  }

/* Random numbers in halfspinors */
  for(int col=0; col < 3; col++) { 
    for(int spin4=0; spin4 < 4; spin4++) { 
      for(int reim=0; reim < 2; reim++) { 
	rs[spin4][col][reim] =  (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }
    }
  }
  
  
  // Random matrix 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) { 
      for(int reim=0; reim < 2; reim++) { 
	matrix[col1][col2][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }
    }
  }

  halfspinor_array *s3=&hspinor;
  u_mat_array* up=&matrix;

  _sse_load((*s3)[0]);
  _sse_su3_multiply(*up);      
  _sse_24_1_gamma3_minus_set();


  _sse_load((*s3)[1]);
  _sse_su3_multiply(*up);
  _sse_24_2_gamma3_minus_set();


  mvv_recons_gamma3_plus_add_store(hspinor, matrix, rs, result2);
  
  for(int col=0; col < 3; col++) { 
      for(int reim=0; reim < 2; reim++) { 
	for(int spin4=0; spin4 < 4; spin4++) { 
	  double diff = result1[spin4][col][reim] - result2[spin4][col][reim];
	  diff /= (double)(4*3*2) ;
#if 0
	  QDPIO::cout << "  col=" << col 
		    << " sp=" << spin4
		    << " re=" << reim 
		    << " diff upper = " << diff 
		    << std::endl;
#endif
	  assertion( fabs(diff) < 1.0e-17 );
	}
      }
  }
}



void testMvvRecons0Minus::run(void) 
{
  halfspinor_array hspinor ALIGN;
  u_mat_array matrix;
  spinor_array rs; // Has to be the goddamn result for the
                   // assembler

  spinor_array result2 ALIGN;


  /* Random numbers in halfspinors */
  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	hspinor[spin2][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;

      }
    }
  }

/* Random numbers in halfspinors */
  for(int col=0; col < 3; col++) { 
    for(int spin4=0; spin4 < 4; spin4++) { 
      for(int reim=0; reim < 2; reim++) { 
	rs[spin4][col][reim] =0;
	result2[spin4][col][reim] =0;
      }
    }
  }
  
  
  // Random matrix 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) { 
      for(int reim=0; reim < 2; reim++) { 
	matrix[col1][col2][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }
    }
  }

  halfspinor_array *s3=&hspinor;
  u_mat_array* up=&matrix;

  _sse_load((*s3)[0]);
  _sse_su3_multiply(*up);
  _sse_24_1_gamma0_plus_set();

   _sse_load(hspinor[1]);
   _sse_su3_multiply(*up);
   _sse_24_2_gamma0_plus_set();


  mvv_recons_gamma0_minus(hspinor, matrix, result2);
  
  for(int col=0; col < 3; col++) { 
      for(int reim=0; reim < 2; reim++) { 
	for(int spin4=0; spin4 < 4; spin4++) { 
	  double diff = rs[spin4][col][reim] - result2[spin4][col][reim];
	  diff /= (double)(4*3*2) ;
#if 0
	  QDPIO::cout << "  col=" << col 
		    << " sp=" << spin4
		    << " re=" << reim 
		    << " diff upper = " << diff 
		    << std::endl;
#endif
	  assertion( fabs(diff) < 1.0e-17 );
	}
      }
  }
}

void testMvvRecons1MinusAdd::run(void) 
{
  halfspinor_array hspinor ALIGN;
  u_mat_array matrix ALIGN;
  spinor_array rs ALIGN; // Has to be the goddamn result for the
                   // assembler

  spinor_array result2 ALIGN;


  /* Random numbers in halfspinors */
  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	hspinor[spin2][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;

      }
    }
  }

/* Random numbers in halfspinors */
  for(int col=0; col < 3; col++) { 
    for(int spin4=0; spin4 < 4; spin4++) { 
      for(int reim=0; reim < 2; reim++) { 
	rs[spin4][col][reim] =  (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
	result2[spin4][col][reim] = rs[spin4][col][reim];
      }
    }
  }
  
  
  // Random matrix 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) { 
      for(int reim=0; reim < 2; reim++) { 
	matrix[col1][col2][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }
    }
  }

  halfspinor_array *s3=&hspinor;
  u_mat_array* up=&matrix;

  _sse_load((*s3)[0]);
  _sse_su3_multiply(*up);      
  _sse_24_1_gamma1_plus();

  _sse_load((*s3)[1]);
  _sse_su3_multiply(*up);
  _sse_24_2_gamma1_plus();


  mvv_recons_gamma1_minus_add(hspinor, matrix, result2);
  
  for(int col=0; col < 3; col++) { 
      for(int reim=0; reim < 2; reim++) { 
	for(int spin4=0; spin4 < 4; spin4++) { 
	  double diff = rs[spin4][col][reim] - result2[spin4][col][reim];
	  diff /= (double)(4*3*2) ;
#if 0
	  QDPIO::cout << "  col=" << col 
		    << " sp=" << spin4
		    << " re=" << reim 
		    << " diff upper = " << diff 
		    << std::endl;
#endif
	  assertion( fabs(diff) < 1.0e-17 );
	}
      }
  }
}

void testMvvRecons2MinusAdd::run(void) 
{
  halfspinor_array hspinor ALIGN;
  u_mat_array matrix ALIGN;
  spinor_array rs ALIGN; // Has to be the goddamn result for the
                   // assembler

  spinor_array result2 ALIGN;


  /* Random numbers in halfspinors */
  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	hspinor[spin2][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;

      }
    }
  }

/* Random numbers in halfspinors */
  for(int col=0; col < 3; col++) { 
    for(int spin4=0; spin4 < 4; spin4++) { 
      for(int reim=0; reim < 2; reim++) { 
	rs[spin4][col][reim] =  (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
	result2[spin4][col][reim] = rs[spin4][col][reim];
      }
    }
  }
  
  
  // Random matrix 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) { 
      for(int reim=0; reim < 2; reim++) { 
	matrix[col1][col2][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }
    }
  }

  halfspinor_array *s3=&hspinor;
  u_mat_array* up=&matrix;

  _sse_load((*s3)[0]);
  _sse_su3_multiply(*up);      
  _sse_24_1_gamma2_plus();

  _sse_load((*s3)[1]);
  _sse_su3_multiply(*up);
  _sse_24_2_gamma2_plus();


  mvv_recons_gamma2_minus_add(hspinor, matrix, result2);
  
  for(int col=0; col < 3; col++) { 
      for(int reim=0; reim < 2; reim++) { 
	for(int spin4=0; spin4 < 4; spin4++) { 
	  double diff = rs[spin4][col][reim] - result2[spin4][col][reim];
	  diff /= (double)(4*3*2) ;
#if 0
	  QDPIO::cout << "  col=" << col 
		    << " sp=" << spin4
		    << " re=" << reim 
		    << " diff upper = " << diff 
		    << std::endl;
#endif
	  assertion( fabs(diff) < 1.0e-17 );
	}
      }
  }
}

void testMvvRecons2MinusAddStore::run(void) 
{
  halfspinor_array hspinor ALIGN;
  u_mat_array matrix ALIGN;
  spinor_array rs ALIGN; // Has to be the goddamn result for the
                   // assembler
  spinor_array result1 ALIGN;
  spinor_array *rn = &result1;
  
  spinor_array result2 ALIGN;


  /* Random numbers in halfspinors */
  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	hspinor[spin2][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;

      }
    }
  }

/* Random numbers in halfspinors */
  for(int col=0; col < 3; col++) { 
    for(int spin4=0; spin4 < 4; spin4++) { 
      for(int reim=0; reim < 2; reim++) { 
	rs[spin4][col][reim] =  (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }
    }
  }
  
  
  // Random matrix 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) { 
      for(int reim=0; reim < 2; reim++) { 
	matrix[col1][col2][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }
    }
  }

  halfspinor_array *s3=&hspinor;
  u_mat_array* up=&matrix;

  _sse_load((*s3)[0]);
  _sse_su3_multiply(*up);      
  _sse_24_1_gamma2_plus_set();


  _sse_load((*s3)[1]);
  _sse_su3_multiply(*up);
  _sse_24_2_gamma2_plus_set();


  mvv_recons_gamma2_minus_add_store(hspinor, matrix, rs, result2);
  
  for(int col=0; col < 3; col++) { 
      for(int reim=0; reim < 2; reim++) { 
	for(int spin4=0; spin4 < 4; spin4++) { 
	  double diff = result1[spin4][col][reim] - result2[spin4][col][reim];
	  diff /= (double)(4*3*2) ;
#if 0
	  QDPIO::cout << "  col=" << col 
		    << " sp=" << spin4
		    << " re=" << reim 
		    << " diff upper = " << diff 
		    << std::endl;
#endif
	  assertion( fabs(diff) < 1.0e-17 );
	}
      }
  }
}

void testMvvRecons3MinusAddStore::run(void) 
{
  halfspinor_array hspinor ALIGN;
  u_mat_array matrix ALIGN;
  spinor_array rs ALIGN; // Has to be the goddamn result for the
                   // assembler
  spinor_array result1 ALIGN;
  spinor_array *rn = &result1;
  
  spinor_array result2 ALIGN;


  /* Random numbers in halfspinors */
  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	hspinor[spin2][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;

      }
    }
  }

/* Random numbers in halfspinors */
  for(int col=0; col < 3; col++) { 
    for(int spin4=0; spin4 < 4; spin4++) { 
      for(int reim=0; reim < 2; reim++) { 
	rs[spin4][col][reim] =  (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }
    }
  }
  
  
  // Random matrix 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) { 
      for(int reim=0; reim < 2; reim++) { 
	matrix[col1][col2][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }
    }
  }

  halfspinor_array *s3=&hspinor;
  u_mat_array* up=&matrix;

  _sse_load((*s3)[0]);
  _sse_su3_multiply(*up);      
  _sse_24_1_gamma3_plus_set();


  _sse_load((*s3)[1]);
  _sse_su3_multiply(*up);
  _sse_24_2_gamma3_plus_set();


  mvv_recons_gamma3_minus_add_store(hspinor, matrix, rs, result2);
  
  for(int col=0; col < 3; col++) { 
      for(int reim=0; reim < 2; reim++) { 
	for(int spin4=0; spin4 < 4; spin4++) { 
	  double diff = result1[spin4][col][reim] - result2[spin4][col][reim];
	  diff /= (double)(4*3*2) ;
#if 0
	  QDPIO::cout << "  col=" << col 
		    << " sp=" << spin4
		    << " re=" << reim 
		    << " diff upper = " << diff 
		    << std::endl;
#endif
	  assertion( fabs(diff) < 1.0e-17 );
	}
      }
  }
}
