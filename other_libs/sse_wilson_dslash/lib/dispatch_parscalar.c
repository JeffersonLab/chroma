#include "dispatch_parscalar.h"

#ifdef __cplusplus
extern "C" { 
#endif


#ifdef DSLASH_USE_QMT_THREADS
  /* Threaded version of the dispatch. We call the qmt_call routine
     with our func, n_sites, and argument */
#include <qmt.h>

void dispatch_to_threads(void (*func)(size_t, size_t,int, const void *),
			 spinor_array* the_spinor,
			 halfspinor_array* the_halfspinor, 
			 my_mat_array u,
			 int cb,
			 int n_sites)
{
  ThreadWorkerArgs a;

  a.spinor = the_spinor;
  a.half_spinor = the_halfspinor;
  a.u = u;
  a.cb = cb; 
  qmt_call((qmt_userfunc_t)func, n_sites, &a);
}

#else

#ifdef DSLASH_USE_OMP_THREADS
#warning Using OpenMP dispatch
#include <omp.h>

  /* OpenMP dispatch */
void dispatch_to_threads(void (*func)(size_t, size_t, int, const void*),
			 spinor_array* the_spinor,
			 halfspinor_array* the_halfspinor, 
			 my_mat_array u,
			 int cb,
			 int n_sites)
{
  ThreadWorkerArgs a;
  
  int threads_num;
  int myId;
  int low;
  int high;

  a.spinor = the_spinor;
  a.half_spinor = the_halfspinor;
  a.u = u;
  a.cb = cb; 
  //(*func)(0, n_sites, 0, &a);

  #pragma omp parallel shared(func, n_sites, a) \
      private(threads_num, myId, low, high) default(none)
    {

      threads_num = omp_get_num_threads();
      myId = omp_get_thread_num();
      low = n_sites * myId / threads_num;
      high = n_sites * (myId+1) / threads_num;
      (*func)(low, high, myId, &a);
    }
  
}

#else

  /* Unthreaded dispatch. We call the function directly. The 'low i
ndex'
     is the first site, and the 'thread' should do all the sites */
void dispatch_to_threads(void (*func)(size_t, size_t, int, const void*),
			 spinor_array* the_spinor,
			 halfspinor_array* the_halfspinor, 
			 my_mat_array u,
			 int cb,
			 int n_sites)
{
  ThreadWorkerArgs a;

  a.spinor = the_spinor;
  a.half_spinor = the_halfspinor;
  a.u = u;
  a.cb = cb; 
  (*func)(0, n_sites, 0, &a);
}

#endif
#endif





#ifdef __cplusplus
};
#endif
