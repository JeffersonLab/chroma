#ifndef RECONS_MVV_32BIT
#define RECONS_MVV_32BIT

#include "sse_config.h"
#include "types32.h"

#ifdef __cplusplus
extern "C" { 
#endif

  void mvv_recons_gamma0_plus(  halfspinor_array src, 
			        u_mat_array u,
			      halfspinor_array upper_sum,
			      halfspinor_array lower_sum);

  void mvv_recons_gamma1_plus_add(  halfspinor_array src, 
				    u_mat_array u,
				  halfspinor_array upper_sum,
				  halfspinor_array lower_sum);

  void mvv_recons_gamma2_plus_add(  halfspinor_array src, 
				    u_mat_array u,
				  halfspinor_array upper_sum, 
				  halfspinor_array lower_sum);

  void mvv_recons_gamma2_plus_add_store(  halfspinor_array src, 
					  u_mat_array u,
					  halfspinor_array upper_sum, 
					  halfspinor_array lower_sum,
					spinor_array dst);

  void mvv_recons_gamma3_plus_add_store(  halfspinor_array src, 
					  u_mat_array u,
					  halfspinor_array upper_sum,   halfspinor_array lower_sum,
					spinor_array dst);

  void mvv_recons_gamma0_minus(  halfspinor_array src, 
			         u_mat_array u,
			       halfspinor_array upper_sum,
			       halfspinor_array lower_sum);
  
  void mvv_recons_gamma1_minus_add(  halfspinor_array src, 
				     u_mat_array u,
				   halfspinor_array upper_sum,
				   halfspinor_array lower_sum);
  
  void mvv_recons_gamma2_minus_add(  halfspinor_array src, 
				     u_mat_array u,
				   halfspinor_array upper_sum, halfspinor_array lower_sum);
  
  void mvv_recons_gamma2_minus_add_store(  halfspinor_array src, 
					   u_mat_array u,
					   halfspinor_array upper_sum, 
					   halfspinor_array lower_sum,
					 spinor_array dst);

  void mvv_recons_gamma3_minus_add_store(  halfspinor_array src, 
					   u_mat_array u,
					   halfspinor_array upper_sum, 
					   halfspinor_array lower_sum,
					 spinor_array dst);


#ifdef __cplusplus
};
#endif

#endif
