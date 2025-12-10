#ifndef DISPATCH_PARSCALAR_H
#define DISPATCH_PARSCALAR_H

#include <cstdlib>         /* for size_t */

namespace CPlusPlusWilsonDslash {
  /* Thread worker argument structure */
  struct ThreadWorkerArgs {
    void *spinor;           /*!< Spinor either read */
    void *half_spinor;           /*!< Spinor to  write */
    void *u;        /*!< Gauge field - suitably packed */
    void *s;
    int cb;            /*!< Checkerboard (source) */
  };

  /* Functions: Thread Dispatch */
  void dispatchToThreads(void (*func)(size_t, size_t, int, const void *),
			 void* the_spinor,
			 void* the_halfspinor, 
			 void* u,
			 void* s,
			 int cb,
			 int n_sites);

}; // namespace

namespace CPlusPlusClover {
  /* Thread worker argument structure */
  struct CloverThreadWorkerArgs {
    void *spinor;           /*!< Spinor either read */
    void *half_spinor;      /*!< Spinor to  write */
    void *u;                /*!< Gauge field - suitably packed */
    void *clov;       /*!< Inverse Clover */
    void *spinor2;
    void *s;
    int half;
  };


  /* Functions: Thread Dispatch */
  void dispatchToThreads(void (*func)(size_t, size_t, int, const void *),
			 void* the_spinor,
			 void* the_halfspinor, 
			 void* u,
			 void* clov,
			 void* spinor2,
			 void* s,
			 int half,
			 int n_sites);
    
}; // namespace

#endif
