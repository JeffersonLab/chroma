#include "unittest.h"
#include "testDslash3D.h"

#include "qdp.h"
using namespace QDP;
using namespace std;

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
testDslash3D::run(void) 
{

  LatticeFermionF3 chi, chi2, psi;
  LatticeFermionD3 chid, chi2d, psid;

  // What we consider to be small enough...
  Double small32,small64;
  small32 = Double(1.0e-7);
  small64 = Double(1.0e-16);


  // Make a random gauge field 
  multi1d<LatticeColorMatrixF> u(4);
  multi1d<LatticeColorMatrixD> ud(4);

  for(int mu=0; mu < 4; mu++) { 
       gaussian(u[mu]);
       reunit(u[mu]);
       ud[mu] = u[mu];
       
  }

  // Make a random source
  gaussian(psi);
  gaussian(psid);

  Dslash3D<float> D32(Layout::lattSize().slice(),
		    Layout::QDPXX_getSiteCoords,
		    Layout::QDPXX_getLinearSiteIndex,
		    Layout::QDPXX_nodeNumber);
  

  //  QDPIO::cout << "Need to allocate the packed gauge: "<< 4*Layout::sitesOnNode()*sizeof(PrimitiveSU3MatrixF) << std::endl << std::flush;
  
  multi1d<PrimitiveSU3MatrixF> packed_gauge;
  packed_gauge.resize(4*Layout::sitesOnNode());
  qdp_pack_gauge_3d(u, packed_gauge);

  std::cout << std::endl;

  // Go through the test cases -- apply SSE dslash versus, QDP Dslash 
  for(int isign=1; isign >= -1; isign -=2) {
    for(int cb=0; cb < 2; cb++) { 
      int source_cb = 1 - cb;
      int target_cb = cb;
      chi = zero;
      chi2 = zero;

      // Apply SSE Dslash
      D32((float *)&(chi.elem(0).elem(0).elem(0).real()),	  
	  (float *)&(psi.elem(0).elem(0).elem(0).real()),
	  (float *)&(packed_gauge[0]),
	  isign, 
	  source_cb);
      
      // Apply QDP Dslash
      dslash_3d(chi2,u,psi, isign, target_cb);
      
      // Check the difference per number in chi std::vector
      LatticeFermionF diff = chi2 -chi;

      Double diff_norm = sqrt( norm2( diff ) ) 
	/ ( Real(4*3*2*Layout::vol()) / Real(2));
	
      QDPIO::cout << "\t cb = " << source_cb << "  isign = " << isign << "  diff_norm = " << diff_norm << std::endl;      
      // Assert things are OK...
      assertion( toBool( diff_norm < small32 ) );

    }
  }


  Dslash3D<double> D64(Layout::lattSize().slice(),
			 Layout::QDPXX_getSiteCoords,
			 Layout::QDPXX_getLinearSiteIndex,
			 Layout::QDPXX_nodeNumber);


   /// Pack the gauge fields
   multi1d<PrimitiveSU3MatrixD> packed_gauged __attribute__((aligned(16)));
   packed_gauged.resize( 4 * Layout::sitesOnNode() );
   qdp_pack_gauge_3d(ud, packed_gauged);
   psi = psid;  // Downcasts 

   QDPIO::cout << std::endl;

   // Go through the test cases -- apply SSE dslash versus, QDP Dslash 
   for(int isign=1; isign >= -1; isign -=2) {
     for(int cb=0; cb < 2; cb++) { 
      int source_cb = 1 - cb;
      int target_cb = cb;
      chid = zero;
      chi2= zero;

      // Apply SSE Dslash
      D64((double *)&(chid.elem(0).elem(0).elem(0).real()),	  
	  (double *)&(psid.elem(0).elem(0).elem(0).real()),
	  (double *)&(packed_gauged[0]),
	  isign, 
	  source_cb);
      
      // Apply QDP Dslash
      chi = chid;
      //      dslash_3d(chi2,u,psi, isign, target_cb);
      D32((float *)&(chi2.elem(0).elem(0).elem(0).real()),	  
	  (float *)&(psi.elem(0).elem(0).elem(0).real()),
	  (float *)&(packed_gauge[0]),
	  isign, 
	  source_cb);

      
      // Check the difference per number in chi std::vector
      LatticeFermionF diff = chi2 -chi;

      Double diff_norm = sqrt( norm2( diff ) ) 
	/ ( Real(4*3*2*Layout::vol()) / Real(2));
	
      QDPIO::cout << "\t cb = " << source_cb << "  isign = " << isign << "  diff_norm = " << diff_norm << std::endl;      
      // Assert things are OK...
      assertion( toBool( diff_norm < small32 ) );

    }
  }

}
