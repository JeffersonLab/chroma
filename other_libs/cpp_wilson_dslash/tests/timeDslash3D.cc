#include "unittest.h"
#include "timeDslash3D.h"

#include "qdp.h"
using namespace QDP;

#undef PAT
#ifdef PAT
#include <pat_api.h>
#endif

#ifndef DSLASH_3D_W_H
#include "dslash_3d_w.h"
#endif

#ifndef REUNIT_H
#include "reunit.h"
#endif

#include "cpp_dslash.h"
#include "cpp_dslash_qdp_packer.h"

using namespace Assertions;
using namespace CPlusPlusWilsonDslash;

void
timeDslash3D::run(void) 
{

  // If we have openmp then do this
  LatticeFermionF chi, psi;
  LatticeFermionD chid, psid;

  // Make a random gauge field 
  multi1d<LatticeColorMatrixF> u(4);
  multi1d<LatticeColorMatrixD> ud(4);

  for(int mu=0; mu < 4; mu++) { 
    gaussian(u[mu]);
    reunit(u[mu]);
    ud[mu]=u[mu];
  }

  // Make a random source
  gaussian(psi);
  gaussian(psid);
  
  // Initialize the wilson dslash
  Dslash3D<float> D32(Layout::lattSize().slice(),
		    Layout::QDPXX_getSiteCoords,
		    Layout::QDPXX_getLinearSiteIndex,
		    Layout::QDPXX_nodeNumber);
  // Initialize the wilson dslash
  Dslash3D<double> D64(Layout::lattSize().slice(),
		    Layout::QDPXX_getSiteCoords,
		    Layout::QDPXX_getLinearSiteIndex,
		    Layout::QDPXX_nodeNumber);
  
  /// Pack the gauge fields
  multi1d<PrimitiveSU3MatrixF> packed_gauge;
  packed_gauge.resize( 4 * Layout::sitesOnNode() );
  qdp_pack_gauge_3d(u, packed_gauge);

  multi1d<PrimitiveSU3MatrixD> packed_gauged;
  packed_gauged.resize( 4 * Layout::sitesOnNode() );
  qdp_pack_gauge_3d(ud, packed_gauged);
 

  QDPIO::cout << std::endl;

  StopWatch swatch;
  double time=0;
  double n_secs = 25;
  int iters=65000;
  QDPIO::cout << std::endl;
  QDPIO::cout << "\t Timing with " << iters << " counts" << std::endl;


  swatch.reset();
  swatch.start();

#ifdef PAT
  int ierr;
  ierr=PAT_region_begin(19, "DslashLoop");
#endif
  for(int i=0; i < iters; ++i) {
    D32( (float *)&(chi.elem(0).elem(0).elem(0).real()),
	 (float *)&(psi.elem(0).elem(0).elem(0).real()),
	 (float *)&(packed_gauge[0]),
	 1, 0);

  }
#ifdef PAT
  ierr=PAT_region_end(19);
#endif

  swatch.stop();
  time=swatch.getTimeInSeconds();

  // Average time over nodes
  QDPInternal::globalSum(time);
  time /= (double)Layout::numNodes();

  QDPIO::cout << "\t " << iters << " iterations in " << time << " seconds " << std::endl;
  QDPIO::cout << "\t " << 1.0e6*time/(double)iters << " u sec/iteration" << std::endl;    
  double Mflops = 984.0f*(double)(iters)*(double)(Layout::vol()/2)/1.0e6;
  double perf = Mflops/time;
  QDPIO::cout << "\t Performance is: " << perf << " Mflops (sp) in Total" << std::endl;
  QDPIO::cout << "\t Performance is: " << perf / (double)Layout::numNodes() << " per MPI Process" << std::endl;
  QDPIO::cout << std::endl;


  QDPIO::cout << "\t Timing with " << iters << " counts" << std::endl;

  swatch.reset();
  swatch.start();
  
  for(int i=0; i < iters; ++i) {
    D64(  (double *)&(chid.elem(0).elem(0).elem(0).real()),
	  (double *)&(psid.elem(0).elem(0).elem(0).real()),
	  (double *)&(packed_gauged[0]),
	   -1, 0);

  }
  swatch.stop();
  time=swatch.getTimeInSeconds();

  // Average time over nodes
  QDPInternal::globalSum(time);
  time /= (double)Layout::numNodes();

  QDPIO::cout << "\t " << iters << " iterations in " << time << " seconds " << std::endl;
  QDPIO::cout << "\t " << 1.0e6*time/(double)iters << " u sec/iteration" << std::endl;    
  Mflops = 984.0f*(double)(iters)*(double)(Layout::vol()/2)/1.0e6;
  perf = Mflops/time;
  QDPIO::cout << "\t Performance is: " << perf << " Mflops (dp) in Total" << std::endl;
  QDPIO::cout << "\t Performance is: " << perf / (double)Layout::numNodes() << " per MPI Process" << std::endl;

  // Finalize the Dslash
}
