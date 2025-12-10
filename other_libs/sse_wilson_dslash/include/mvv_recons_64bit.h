#ifndef RECONS_MVV_64BIT
#define RECONS_MVV_64BIT

#include <sse_config.h>
#include "types64.h"

#ifdef __cplusplus
extern "C" { 
#endif

void mvv_recons_gamma0_plus( halfspinor_array src, 
			     u_mat_array u,
			     spinor_array dst);

void mvv_recons_gamma1_plus_add( halfspinor_array src, 
				 u_mat_array u,
				spinor_array dst);

void mvv_recons_gamma2_plus_add( halfspinor_array src, 
				 u_mat_array u,
				spinor_array dst);

void mvv_recons_gamma2_plus_add_store( halfspinor_array src, 
				       u_mat_array u,
				       spinor_array sum,
				       spinor_array dst);

void mvv_recons_gamma3_plus_add_store( halfspinor_array src, 
			     u_mat_array u,
			     spinor_array sum,
			    spinor_array dst);


void mvv_recons_gamma0_minus( halfspinor_array src, 
			      u_mat_array u,
			      spinor_array dst);

void mvv_recons_gamma1_minus_add( halfspinor_array src, 
				 u_mat_array u,
				spinor_array dst);

void mvv_recons_gamma2_minus_add( halfspinor_array src, 
				 u_mat_array u,
				spinor_array dst);

void mvv_recons_gamma2_minus_add_store( halfspinor_array src, 
					u_mat_array u,
					spinor_array sum,
					spinor_array dst);

void mvv_recons_gamma3_minus_add_store( halfspinor_array src, 
			     u_mat_array u,
			     spinor_array sum,
			    spinor_array dst);

#ifdef __cplusplus
};
#endif

#endif
