#ifndef RECONS
#define RECONS

#include "sse_config.h"

#if SSE_PRECISION == 32
#include "types32.h"
#else 
#include "types64.h"
#endif

#ifdef __cplusplus
extern "C" { 
#endif

void recons_4dir_plus( halfspinor_array hs0,
		       halfspinor_array hs1,
		       halfspinor_array hs2,
		       halfspinor_array hs3,
		      spinor_array spinor);

void recons_3dir_plus( halfspinor_array hs0,
		       halfspinor_array hs1,
		       halfspinor_array hs2,
		       spinor_array spinor);

void recons_4dir_minus( halfspinor_array hs0,
		        halfspinor_array hs1,
		        halfspinor_array hs2,
		        halfspinor_array hs3,
		       spinor_array spinor);

void recons_3dir_minus( halfspinor_array hs0,
		        halfspinor_array hs1,
		        halfspinor_array hs2,
			spinor_array spinor);



#ifdef __cplusplus
};
#endif

#endif
