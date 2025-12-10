#include <dispatch_scalar.h>
#include <shift_table_scalar.h>

namespace CPlusPlusWilsonDslash {

void dispatchToThreads(void (*func)(size_t, size_t, int, const void *),
		       void* source,
		       void* result, 
		       void *u,
		       void *s,
		       int cb,
		       int n_sites)
{
  struct ThreadWorkerArgs a;
  a.psi = source;
  a.res =result;
  a.u = u;
  a.cb = cb;
  a.s = s;
  /* Call dispatch function, with lo=0, hi=n_sites */
  (*func)(0, n_sites, 0, &a);
}

  
}; // End Namespace


namespace CPlusPlusClover { 

  void dispatchToThreads(void (*func)(size_t, size_t, int, const void *),
			 void* source,
			 void* result, 
			 void* u,
			 void* invclov_ee,
			 void* clov_oo,
			 void* t_spinor,
			 void* s,
			 int n_sites) 
  {
    struct CloverThreadWorkerArgs a;
    a.psi = source;
    a.res =result;
    a.u = u;
    a.invclov_ee = invclov_ee;
    a.clov_oo = clov_oo;
    a.t_spinor = t_spinor;
    a.s = s;
    /* Call dispatch function, with lo=0, hi=n_sites */
    (*func)(0, n_sites, 0, &a);
    
  }

};
