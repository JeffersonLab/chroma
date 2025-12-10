#ifndef DISPATCH_PARSCALAR_H
#define DISPATCH_PARSCALAR_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h> /* For size_t */
#include <sse_config.h> /* For info about the precision */

#if SSE_PRECISION == 32
#include <types32.h>    /* Sadly we need the overlay types here */
#else
#include <types64.h>    /* Sadly we need the overlay types here */
#endif

  typedef struct {
    spinor_array* spinor;             /*!< Spinor either read or write */
    halfspinor_array *half_spinor;         /*!< Half Spinor - either read or write */
    u_mat_array       (*u)[4];      /*!< Gauge field - suidably packed */
    int cb;                         /*!< Checkerboard (source) */
  } ThreadWorkerArgs;

  
  /*! The thread dispatcher. In a 64 bit ABI it should be cheap to call
   * as a function call.
   \param func  The function to be dispatched to the thread of type conforming
                to qmt_userfunc
   \param the_spinor  The 4 spinor array that is either the source or target
                      of the operation
   \param the_halfspinor The 2-spinor array that is either the target or
                         source of the operation

   \param u              A pointer to the suitably packed gauge field
   \param n_sites        The total number of sites that we want to divvy up
                         between the threads. (In this dslash, 
			 typically this is the number
			 of sites per checkerboard.)
  */

  void dispatch_to_threads(void (*func) (size_t, size_t, int, const void*),
			   spinor_array* the_spinor,
			   halfspinor_array* the_halfspinor, 
			   my_mat_array u,
			   int cb,
			   int n_sites);
  
  
#ifdef __cplusplus
};
#endif

#endif
