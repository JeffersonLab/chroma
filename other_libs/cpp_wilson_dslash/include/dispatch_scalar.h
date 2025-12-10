#ifndef DISPATCH_SCALAR_H
#define DISPATCH_SCALAR_H

#include <cstdlib>         /* for size_t */

#include <shift_table_scalar.h>
namespace CPlusPlusWilsonDslash {


  /* This is absolutely horrible -- There must be a better way */

  /* Thread worker argument structure */
  struct ThreadWorkerArgs {
    void *res;           /*!< Spinor either read */
    void *psi;           /*!< Spinor to  write */
    void *u;             /*!< Gauge field - suitably packed */
    void *s;             /* Shift table */
    int cb;              /*!< Checkerboard (source) */

  };

  /* Functions: Thread Dispatch */
  void dispatchToThreads(void (*func)(size_t, size_t, int, const void *),
			 void* source,
			 void* result, 
			 void* u,
			 void* s,
			 int cb,
			 int n_sites);

}; // namespace

namespace CPlusPlusClover { 

  struct CloverThreadWorkerArgs { 
    void *res;
    void *psi;
    void *u;
    void *invclov_ee;
    void *clov_oo;
    void *t_spinor;
    void *s;
  };
  void dispatchToThreads(void (*func)(size_t, size_t, int, const void *),
			 void* source,
			 void* result, 
			 void* u,
			 void* invclov_ee,
			 void* clov_oo,
			 void* t_spinor,
			 void* s,
			 int n_sites);

};

#endif
