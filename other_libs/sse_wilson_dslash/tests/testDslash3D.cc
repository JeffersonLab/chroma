#include "unittest.h"
#include "testDslash3D.h"

#include "qdp.h"
using namespace QDP;

#ifndef DSLASH_3D_W_H
#include "dslash_3d_w.h"
#endif

#ifndef REUNIT_H
#include "reunit.h"
#endif

#include "sse_dslash_3d.h"
#include "sse_dslash_qdp_packer_3d.h"

using namespace Assertions;

#ifdef DSLASH_USE_OMP_THREADS
#include <omp.h>
#endif


void
testDslash3D::run(void) 
{

  // If we have openmp then do this
#ifdef DSLASH_USE_OMP_THREADS
  int threads_num;
  int myId;

#pragma omp parallel private(threads_num, myId) default(none)
  {
    threads_num = omp_get_num_threads();
    myId = omp_get_thread_num();
    if ( myId == 0 ) { 
      printf("\nRunning with %d OpenMP threads\n", threads_num);
    }
  }
#endif

  LatticeFermion chi, chi2, psi;

  // What we consider to be small enough...
  Double small;
  if ( sizeof(SSEREAL) == 4 ) { 
    small = Double(1.0e-7);
  }
  else {
    // Adjust this...
    // Adjust this...
#ifdef SSEDSLASH_SLOPPY
    small = Double(1.0e-7);
#else
    small = Double(1.0e-16);
#endif
  }

  // Make a random gauge field 
  multi1d<LatticeColorMatrix> u(4);

  for(int mu=0; mu < 4; mu++) { 
    gaussian(u[mu]);
    reunit(u[mu]);
  }

  // Make a random source
  gaussian(psi);

  
  // Initialize the wilson dslash
  init_sse_su3dslash_3d(Layout::lattSize().slice(),
			Layout::QDPXX_getSiteCoords,
			Layout::QDPXX_getLinearSiteIndex,
			Layout::QDPXX_nodeNumber);

  /// Pack the gauge fields
  multi1d<SSEDslash3D::PrimitiveSU3Matrix> packed_gauge;
  packed_gauge.resize( 4 * Layout::sitesOnNode() );
  SSEDslash3D::qdp_pack_gauge_3d(u, packed_gauge);
 
  QDPIO::cout << std::endl;

  // Go through the test cases -- apply SSE dslash versus, QDP Dslash 
  for(int isign=1; isign >= -1; isign -=2) {
    for(int cb=0; cb < 2; cb++) { 
      int source_cb = 1 - cb;
      int target_cb = cb;
      chi = zero;
      chi2 = zero;

      // Apply SSE Dslash
      sse_su3dslash_wilson_3d((SSEREAL *)&(packed_gauge[0]),
			      (SSEREAL *)&(psi.elem(0).elem(0).elem(0).real()),
			      (SSEREAL *)&(chi.elem(0).elem(0).elem(0).real()),
			      isign, source_cb);
      
      // Apply QDP Dslash
      dslash_3d(chi2,u,psi, isign, target_cb);
      
      // Check the difference per number in chi std::vector
      LatticeFermion diff= zero;
      diff = chi2 -chi;

      Double diff_norm = sqrt( norm2( diff ) )
	/ ( Real(4*3*2*Layout::vol()) / Real(2));
	
      QDPIO::cout << "\t cb = " << source_cb << "  isign = " << isign << "  diff_norm = " << diff_norm << std::endl;      
      // Assert things are OK...
      assertion( toBool( diff_norm < small ) );

    }
  }

  // Finalize the Dslash
  free_sse_su3dslash_3d();

}
