#include "unittest.h"
#include "timeDecomp.h"

#include "qdp.h"
using namespace QDP;

#ifndef DSLASH_M_W_H
#include "dslashm_w.h"
#endif

#ifndef REUNIT_H
#include "reunit.h"
#endif

#include "sse_dslash.h"
#include "sse_dslash_qdp_packer.h"


/* Use these for testing */
#include <sse_config.h>
#include <sse_align.h>
#include <types32.h>
#include <dispatch_parscalar.h>

extern "C" {
extern  halfspinor_array* chi1;
extern  halfspinor_array* chi2;
extern  void decomp_plus(size_t lo,size_t hi, int id, const void *ptr);
};
/* This is to call the dispatcher */

using namespace Assertions;

void
timeDecomp::run(void) 
{
  /* The fermion to decompose */
  LatticeFermion psi;
  int subgrid_vol_cb=Layout::sitesOnNode()/2;
  int cb=0;
  // Make a random gauge field 
  // Strictly this won't be used...
  multi1d<LatticeColorMatrix> u(4);

  for(int mu=0; mu < 4; mu++) { 
    gaussian(u[mu]);
    reunit(u[mu]);
  }

  // Make a random source
  gaussian(psi);

  
  // Initialize the wilson dslash
  init_sse_su3dslash(Layout::lattSize().slice(),
		     Layout::QDPXX_getSiteCoords,
		     Layout::QDPXX_getLinearSiteIndex,
		     Layout::QDPXX_nodeNumber);

  /// Pack the gauge fields
  multi1d<SSEDslash::PrimitiveSU3Matrix> packed_gauge;
  packed_gauge.resize( 4 * Layout::sitesOnNode() );
  SSEDslash::qdp_pack_gauge(u, packed_gauge);
 
  QDPIO::cout << std::endl;

  StopWatch swatch;
  double n_secs = 10;
  int iters=1;
  double time=0;
  QDPIO::cout << std::endl << "\t Calibrating for " << n_secs << " seconds " << std::endl;
  do {
    swatch.reset();
    swatch.start();
    for(int i=0; i < iters; i++) { 

      dispatch_to_threads(decomp_plus,
			  (float(*)[4][3][2])&(psi.elem(0).elem(0).elem(0).real()),
			  chi1,
			  (float (*)[4][3][3][2])&(packed_gauge[0]),
			  cb,
			  subgrid_vol_cb);



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

  swatch.reset();
  swatch.start();
  
  for(int i=0; i < iters; ++i) {
    dispatch_to_threads(decomp_plus,
			(float(*)[4][3][2])&(psi.elem(0).elem(0).elem(0).real()),
			chi1,
			(float (*)[4][3][3][2])&(packed_gauge[0]),
			cb,
			subgrid_vol_cb);

  }
  swatch.stop();
  time=swatch.getTimeInSeconds();

  // Average time over nodes
  QDPInternal::globalSum(time);
  time /= (double)Layout::numNodes();

  QDPIO::cout << "\t " << iters << " iterations in " << time << " seconds " << std::endl;
  QDPIO::cout << "\t " << 1.0e6*time/(double)iters << " u sec/iteration" << std::endl;    

  // Finalize the Dslash
  free_sse_su3dslash();

}
