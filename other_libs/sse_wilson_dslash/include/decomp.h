#ifndef DECOMP_H
#define DECOMP_H

#include "sse_config.h"

#if SSE_PRECISION == 32
#include "types32.h"
#else 
#include "types64.h"
#endif

#ifdef __cplusplus
extern "C" { 
#endif

void decomp_gamma0_minus(  spinor_array src, halfspinor_array dst);
void decomp_gamma1_minus(  spinor_array src, halfspinor_array dst);
void decomp_gamma2_minus(  spinor_array src, halfspinor_array dst);  
void decomp_gamma3_minus(  spinor_array src, halfspinor_array dst);

void decomp_gamma0_plus(  spinor_array src, halfspinor_array dst);
void decomp_gamma1_plus(  spinor_array src, halfspinor_array dst);
void decomp_gamma2_plus(  spinor_array src, halfspinor_array dst);  
void decomp_gamma3_plus(  spinor_array src, halfspinor_array dst);


#ifdef __cplusplus
};
#endif

#endif
