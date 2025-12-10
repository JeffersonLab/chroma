#include "unittest.h"
#include "timeDslash3D.h"

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

/* Cray PAT Performance Analysis tool */
#undef PAT  
#ifdef PAT
#include <pat_api.h>
#endif

void
timeDslash3D::run(void) 
{

  // If we have openmp then do this
#ifdef DSLASH_USE_OMP_THREADS
  int threads_num;
  int myId;
  if ( Layout::nodeNumber() == 0 ) { 
#pragma omp parallel private(threads_num, myId) default(none)
    {
      threads_num = omp_get_num_threads();
      myId = omp_get_thread_num();
      if ( myId == 0 ) { 
	printf("\nRunning with %d OpenMP threads\n", threads_num);
      }
    }
  }
#endif


  LatticeFermion chi, chi2, psi;

  // What we consider to be small enough...
  Double small;
  if ( sizeof(SSEREAL) == 4 ) { 
    small = Double(5.0e-9);
  }
  else {
    // Adjust this...
    small = Double(1.0e-17);
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
  StopWatch swatch;
  double time=0;
  int iters=131000;
  double n_secs=40;

#if 0
  {
    iters=1;
    QDPIO::cout << std::endl << "\t Calibrating for " << n_secs << " seconds " << std::endl;
    do {
      swatch.reset();
      swatch.start();
      for(int i=0; i < iters; i++) { 
	sse_su3dslash_wilson_3d((SSEREAL *)&(packed_gauge[0]),
				(SSEREAL *)&(psi.elem(0).elem(0).elem(0).real()),
				(SSEREAL *)&(chi.elem(0).elem(0).elem(0).real()),
				1, 0);
      }
      swatch.stop();
      time=swatch.getTimeInSeconds();
      
      // Average time over nodes
      QDPInternal::globalSum(time);
      time /= (double)Layout::numNodes();
      
      if (time < n_secs) {
	iters *=2;
	QDPIO::cout << "." << std::flush;
      }
    }
    while ( time < (double)n_secs );
    
    QDPIO::cout << std::endl;
    QDPIO::cout << "\t Timing with " << iters << " counts" << std::endl;
   }
#endif 
   {
    swatch.reset();
    double dummy=5;
    QDPInternal::globalSum(dummy);

#ifdef PAT
    int ierr;
    ierr=PAT_region_begin(21, "DslashLoop3d");
#endif
    swatch.start();

    for(int i=0; i < iters; ++i) {
      sse_su3dslash_wilson_3d((SSEREAL *)&(packed_gauge[0]),
			      (SSEREAL *)&(psi.elem(0).elem(0).elem(0).real()),
			      (SSEREAL *)&(chi.elem(0).elem(0).elem(0).real()),
			      1, 0);
      
    }
    swatch.stop();
#ifdef PAT
  ierr=PAT_region_end(21);
#endif
    time=swatch.getTimeInSeconds();
    
    // Average time over nodes
    QDPInternal::globalSum(time);
    time /= (double)Layout::numNodes();
    
    QDPIO::cout << "\t " << iters << " iterations in " << time << " seconds " << std::endl;
    QDPIO::cout << "\t " << 1.0e6*time/(double)iters << " u sec/iteration" << std::endl;    
    // Full 4D dslash is 1390 Mflops. 3D Dslash is 3/4*1390~1042.5 ? */
    double Mflops = (1043.0f*(double)(iters)*(double)(Layout::vol()/2))/1.0e6;
    double perf = Mflops/time;
    QDPIO::cout << "\t Performance is: " << perf << " Mflops in Total" << std::endl;
    QDPIO::cout << "\t Performance is: " << perf / (double)Layout::numNodes() << " per MPI Process" << std::endl;
    
    // Finalize the Dslash
    free_sse_su3dslash_3d();
  }

}
