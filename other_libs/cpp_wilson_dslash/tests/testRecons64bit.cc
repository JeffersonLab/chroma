#include "testRecons.h"
#include "unittest.h"
#include <cmath>
#include "recons.h"
#include "sse64.h"
#include "types64.h"
#include "sse_align.h"

using namespace QDP;
using namespace Assertions;



void 
testRecons4DirPlus::run(void)
{
  halfspinor_array hs0, hs1, hs2, hs3;
  spinor_array spinor1, spinor2;

  for(int col=0; col < 3; col++) {
    for(int reim=0; reim < 2; reim++) { 
      for(int spin2=0; spin2 < 2; spin2++) {
	hs0[spin2][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
	hs1[spin2][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
	hs2[spin2][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
	hs3[spin2][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }

      for(int spin4=0; spin4 < 4; spin4++) { 
	spinor1[spin4][col][reim]=0;
	spinor2[spin4][col][reim]=0;
      }
    }
  }



    _sse_load_up(hs0[0]);
    _sse_load(spinor1[0]);
    _sse_vector_add();
    _sse_load_up(hs1[0]);
    _sse_vector_add();
    _sse_load_up(hs2[0]);
    _sse_vector_add();
    _sse_load_up(hs3[0]);
    _sse_vector_add();
    _sse_store(spinor1[0]);

    _sse_load_up(hs0[1]);
    _sse_load(spinor1[1]);
    _sse_vector_add();
    _sse_load_up(hs1[1]);
    _sse_vector_add();
    _sse_load_up(hs2[1]);
    _sse_vector_add();
    _sse_load_up(hs3[1]);
    _sse_vector_add();
    _sse_store(spinor1[1]);

    _sse_load_up(hs0[1]);
    _sse_load(spinor1[2]);
    _sse_vector_i_mul_up();
    _sse_vector_sub();
    _sse_load_up(hs1[1]);
    _sse_vector_add();
    _sse_load_up(hs2[0]);
    _sse_vector_i_mul();
    _sse_vector_sub();
    _sse_load_up(hs3[0]);
    _sse_vector_add();
    _sse_store(spinor1[2]);

    _sse_load_up(hs0[0]);
    _sse_load(spinor1[3]);
    _sse_vector_i_mul_up();
    _sse_vector_sub();
    _sse_load_up(hs1[0]);
    _sse_vector_sub();
    _sse_load_up(hs2[1]);
    _sse_vector_i_mul_up();
    _sse_vector_add();
    _sse_load_up(hs3[1]);
    _sse_vector_add();
    _sse_store(spinor1[3]);


   recons_4dir_plus(hs0, hs1, hs2, hs3, spinor2);

   
   for(int col=0; col < 3; col++) { 
     for(int spin=0; spin < 4; spin++) { 
       for(int reim=0; reim < 2; reim++) { 
	 double diff = spinor1[spin][col][reim] - spinor2[spin][col][reim];
	 diff /= (double)(4*3*2);
	 
#if 0
	 QDPIO::cout << "  col=" << col 
		     << " sp=" << spin
		     << " re=" << reim 
		     << " diff upper = " << diff 
		     << std::endl;
#endif
	 assertion( fabs(diff) < 1.0e-18 );
	 
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
	hs0[spin2][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
	hs1[spin2][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
	hs2[spin2][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
	hs3[spin2][col][reim] = (double)(rand() - RAND_MAX/2)/(double)(RAND_MAX/2)*2.0;
      }

      for(int spin4=0; spin4 < 4; spin4++) { 
	spinor1[spin4][col][reim]=0;
	spinor2[spin4][col][reim]=0;
      }
    }
  }


  _sse_load_up(hs0[0]);
  _sse_load(spinor1[0]);
  _sse_vector_add();
  _sse_load_up(hs1[0]);
  _sse_vector_add();
  _sse_load_up(hs2[0]);
  _sse_vector_add();
  _sse_load_up(hs3[0]);
  _sse_vector_add();
  _sse_store(spinor1[0]);

  _sse_load_up(hs0[1]);
  _sse_load(spinor1[1]);
  _sse_vector_add();
  _sse_load_up(hs1[1]);
  _sse_vector_add();
  _sse_load_up(hs2[1]);
  _sse_vector_add();
  _sse_load_up(hs3[1]);
  _sse_vector_add();
  _sse_store(spinor1[1]);

  _sse_load_up(hs0[1]);
  _sse_load(spinor1[2]);
  _sse_vector_i_mul_up();
  _sse_vector_add();
  _sse_load_up(hs1[1]);
  _sse_vector_sub();
  _sse_load_up(hs2[0]);
  _sse_vector_i_mul();
  _sse_vector_add();
  _sse_load_up(hs3[0]);
  _sse_vector_sub();
  _sse_store(spinor1[2]);
	  
  _sse_load_up(hs0[0]);
  _sse_load(spinor1[3]);
  _sse_vector_i_mul_up();
  _sse_vector_add();
  _sse_load_up(hs1[0]);
  _sse_vector_add();
  _sse_load_up(hs2[1]);
  _sse_vector_i_mul_up();
  _sse_vector_sub();
  _sse_load_up(hs3[1]);
  _sse_vector_sub();
  _sse_store(spinor1[3]);




  recons_4dir_minus(hs0, hs1, hs2, hs3, spinor2);

   
   for(int col=0; col < 3; col++) { 
     for(int spin=0; spin < 4; spin++) { 
       for(int reim=0; reim < 2; reim++) { 
	 double diff = spinor1[spin][col][reim] - spinor2[spin][col][reim];
	 diff /= (double)(4*3*2);
#if 0
	 QDPIO::cout << "  col=" << col 
		     << " sp=" << spin
		     << " re=" << reim 
		     << " diff upper = " << diff 
		     << std::endl;
#endif
	 assertion( fabs(diff) < 1.0e-18 );
	 
       }
     }
   }
   
}
