#include "testRecons.h"
#include "unittest.h"
#include <cmath>
#include "recons.h"
#include "sse32.h"
#include "sse_align.h"

using namespace QDP;
using namespace Assertions;

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

void 
testRecons4DirPlus::run(void)
{
  halfspinor_array hs0, hs1, hs2, hs3;
  spinor_array spinor1, spinor2;

  for(int col=0; col < 3; col++) {
    for(int reim=0; reim < 2; reim++) { 
      for(int spin2=0; spin2 < 2; spin2++) {
	hs0[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	hs1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	hs2[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	hs3[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
      }

      for(int spin4=0; spin4 < 4; spin4++) { 
	spinor1[spin4][col][reim]=0;
	spinor2[spin4][col][reim]=0;
      }
    }
  }

  _sse_pair_load(spinor1[0], spinor1[1]);

  _sse_vector_load_up(hs0);
  _sse_vector_add();

  _sse_vector_load_up(hs1);
  _sse_vector_add();
  _sse_vector_load_up(hs2);
  _sse_vector_add();
  _sse_vector_load_up(hs3);
  _sse_vector_add();

  _sse_pair_store(spinor1[0], spinor1[1]);

  _sse_pair_load(spinor1[2], spinor1[3]);
  _sse_vector_load_up(hs0);
  _sse_24_gamma0_plus_add();
  _sse_vector_load_up(hs1);
  _sse_24_gamma1_plus();
  _sse_vector_load_up(hs2);
  _sse_24_gamma2_plus();
  _sse_vector_load_up(hs3);
  _sse_24_gamma3_plus();
  _sse_pair_store(spinor1[2], spinor1[3]);


   recons_4dir_plus(hs0, hs1, hs2, hs3, spinor2);

   

     for(int spin=0; spin < 4; spin++) { 
       for(int col=0; col < 3; col++) { 
       for(int reim=0; reim < 2; reim++) { 
	 float diff = spinor1[spin][col][reim] - spinor2[spin][col][reim];
	 diff /= (float)(4*3*2);

#if 0
	 QDPIO::cout << "  col=" << col 
		     << " sp=" << spin
		     << " re=" << reim 
		     << " diff upper = " << diff 
		     << std::endl;
#endif
	 assertion( diff < 1.0e-9 );
	 
       }
     }
   }
}

void
testRecons4DirMinus::run(void)
{
  halfspinor_array hs0, hs1, hs2, hs3;
  spinor_array spinor1, spinor2;

  for(int col=0; col < 3; col++) {
    for(int reim=0; reim < 2; reim++) { 
      for(int spin2=0; spin2 < 2; spin2++) {
	hs0[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	hs1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	hs2[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	hs3[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
      }

      for(int spin4=0; spin4 < 4; spin4++) { 
	spinor1[spin4][col][reim]=0;
	spinor2[spin4][col][reim]=0;
      }
    }
  }

  _sse_pair_load(spinor1[0], spinor1[1]);

  _sse_vector_load_up(hs0);
  _sse_vector_add();

  _sse_vector_load_up(hs1);
  _sse_vector_add();
  _sse_vector_load_up(hs2);
  _sse_vector_add();
  _sse_vector_load_up(hs3);
  _sse_vector_add();

  _sse_pair_store(spinor1[0], spinor1[1]);

  _sse_pair_load(spinor1[2], spinor1[3]);
  _sse_vector_load_up(hs0);
  _sse_24_gamma0_minus_add();
  _sse_vector_load_up(hs1);
  _sse_24_gamma1_minus();
  _sse_vector_load_up(hs2);
  _sse_24_gamma2_minus();
  _sse_vector_load_up(hs3);
  _sse_24_gamma3_minus();
  _sse_pair_store(spinor1[2], spinor1[3]);


   recons_4dir_minus(hs0, hs1, hs2, hs3, spinor2);

   
   for(int col=0; col < 3; col++) { 
     for(int spin=0; spin < 4; spin++) { 
       for(int reim=0; reim < 2; reim++) { 
	 float diff = spinor1[spin][col][reim] - spinor2[spin][col][reim];
	 diff /= (float)(4*3*2);
	 
#if 0
	 QDPIO::cout << "  col=" << col 
		     << " sp=" << spin
		     << " re=" << reim 
		     << " diff upper = " << diff 
		     << std::endl;
#endif
	 assertion( diff < 1.0e-9 );
	 
       }
     }
   }
   
}
