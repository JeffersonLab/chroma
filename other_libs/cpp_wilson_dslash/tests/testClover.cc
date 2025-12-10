#include "unittest.h"
#include "testClover.h"

#include "qdp.h"
#include <iostream>

using namespace QDP;

#ifndef DSLASH_M_W_H
#include "dslashm_w.h"
#endif

#ifndef REUNIT_H
#include "reunit.h"
#endif

#include "cpp_dslash.h"
#include "cpp_dslash_qdp_packer.h"
#include "cpp_clover.h"
#include "cpp_clover_site_apply_32bit.h"

#include "cache.h"
#include <cmath>

using namespace Assertions;
using namespace CPlusPlusWilsonDslash;
using namespace CPlusPlusClover;
using namespace CPlusPlusClover::CPlusPlusClover32Bit;
using namespace Dslash32BitTypes;

void
testClover::run(void) 
{

  LatticeFermionF3 chi, Dpsi, ADpsi;
  LatticeFermionF3 chi2,psi;
  LatticeFermionD3 psid,chid ;

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
  // Downcast
  psid=psi;

  // Reference
  Dslash<float> D32(Layout::lattSize().slice(),
		    Layout::QDPXX_getSiteCoords,
		    Layout::QDPXX_getLinearSiteIndex,
		    Layout::QDPXX_nodeNumber);

  Dslash<double> D64(Layout::lattSize().slice(),
		    Layout::QDPXX_getSiteCoords,
		    Layout::QDPXX_getLinearSiteIndex,
		    Layout::QDPXX_nodeNumber);

  CloverSchur4D<float> Klov32(Layout::lattSize().slice(),
			   Layout::QDPXX_getSiteCoords,
			   Layout::QDPXX_getLinearSiteIndex,
			   Layout::QDPXX_nodeNumber);


  CloverSchur4D<double> Klov64(Layout::lattSize().slice(),
			   Layout::QDPXX_getSiteCoords,
			   Layout::QDPXX_getLinearSiteIndex,
			   Layout::QDPXX_nodeNumber);

  
  multi1d<PrimitiveSU3MatrixF> packed_gauge __attribute__((aligned(16)));
  packed_gauge.resize(4*Layout::sitesOnNode());
  qdp_pack_gauge(u, packed_gauge);

   /// Pack the gauge fields
  multi1d<PrimitiveSU3MatrixD> packed_gauged __attribute__((aligned(16)));
  packed_gauged.resize( 4 * Layout::sitesOnNode() );
  qdp_pack_gauge(ud, packed_gauged);


  multi1d<Clover32BitTypes::CloverTerm> clov __attribute__((aligned(16)));
  clov.resize(Layout::sitesOnNode());
  multi1d<Clover32BitTypes::CloverTerm> invclov __attribute__((aligned(16)));
  invclov.resize(Layout::sitesOnNode());


  Clover64BitTypes::CloverTerm *xclovd, *clovd, *xinvclovd, *invclovd;
  
  xclovd = (Clover64BitTypes::CloverTerm*)malloc(Layout::sitesOnNode()*sizeof(Clover64BitTypes::CloverTerm)+Cache::CacheLineSize);
  
  unsigned long pad=0;

  if( (unsigned long)xclovd % Cache::CacheLineSize != 0 ) { 
    pad = Cache::CacheLineSize - (unsigned long) xclovd % Cache::CacheLineSize;
  }
  clovd =(Clover64BitTypes::CloverTerm*)( (unsigned char*)xclovd + pad );
  
  xinvclovd = (Clover64BitTypes::CloverTerm*)malloc(Layout::sitesOnNode()*sizeof(Clover64BitTypes::CloverTerm)+Cache::CacheLineSize);

  pad = 0;
  if( (unsigned long)xinvclovd % Cache::CacheLineSize != 0 ) { 
    pad = Cache::CacheLineSize - (unsigned long) xinvclovd % Cache::CacheLineSize;
  }
  invclovd =(Clover64BitTypes::CloverTerm*)( (unsigned char*)xinvclovd + pad );
  
  // Randomize clover term....
  for(int site=0; site < Layout::sitesOnNode(); site++) { 
    for(int j=0; j < 2; j++) { 
      for(int d=0; d < 6; d++) { 
	(clov[site])[j].diag[d] = drand48();
	(invclov[site])[j].diag[d] = drand48();
	(clovd[site])[j].diag[d] = (double)(clov[site])[j].diag[d];
	(invclovd[site])[j].diag[d]=(double)(invclov[site])[j].diag[d];
      }
      for(int od=0; od < 15; od++) { 
	(clov[site])[j].off_diag[od][0] = drand48();
	(clov[site])[j].off_diag[od][1] = drand48();
	(invclov[site])[j].off_diag[od][0] = drand48();
	(invclov[site])[j].off_diag[od][1] = drand48();
	
	(clovd[site])[j].off_diag[od][0]= (double)(clov[site])[j].off_diag[od][0]; 

	(clovd[site])[j].off_diag[od][1] = (double)(clov[site])[j].off_diag[od][1]; 
	(invclovd[site])[j].off_diag[od][0] = (double)(invclov[site])[j].off_diag[od][0];
	(invclovd[site])[j].off_diag[od][1] = (double)(invclov[site])[j].off_diag[od][1];
	
      }
    }
  }

  for(int isign=-1; isign < 2; isign+=2) {

    // Apply clover op to psi into chi
    Klov32((float *)&(chi.elem(0).elem(0).elem(0).real()),	  
	   (float *)&(psi.elem(0).elem(0).elem(0).real()),
	   (float *)&(packed_gauge[0]),
	   (float *)&(clov[0]),
	   (float *)&(invclov[0]),
	   isign); 

    
    
    // Now try 
    // Apply SSE Dslash to do it spearately
    
    // D psi
    D32((float *)&(Dpsi.elem(0).elem(0).elem(0).real()),	  
	(float *)&(psi.elem(0).elem(0).elem(0).real()),
	(float *)&(packed_gauge[0]),
	isign, // Isign = 1 (PLUS)  
	1); // source CB = 1 (odd)
    




    const int *rb0tab = rb[0].siteTable().slice();
    for(int j=0; j < rb[0].siteTable().size(); j++) { 
      
      int site = rb0tab[j];
      FourSpinor* res = (FourSpinor *)&(ADpsi.elem(site));
      FourSpinor* src = (FourSpinor *)&(Dpsi.elem(site));
      
      cloverSiteApply(*res,
		      invclov[site], 
		      *src
		      );
    }



    D32((float *)&(Dpsi.elem(0).elem(0).elem(0).real()),	  
	(float *)&(ADpsi.elem(0).elem(0).elem(0).real()),
	(float *)&(packed_gauge[0]),
	isign, 
	0); // source CB = 0 (even)

    
    
    // Apply chi2 = A_oo psi
    const int *rb1tab = rb[1].siteTable().slice();
    for(int j=0; j  < rb[1].siteTable().size(); j++) { 
      int site = rb1tab[j];
      
      cloverSiteApply(*(Dslash32BitTypes::FourSpinor*)&(chi2.elem(site).elem(0).elem(0).real()),
		      clov[site], 
		      *(Dslash32BitTypes::FourSpinor*)&(psi.elem(site).elem(0).elem(0).real()));
    }
    chi2[rb[1]] -= Dpsi;
    LatticeFermionF3 diff_float = chi2 - chi;
    QDPIO::cout << std::endl;
    QDPIO::cout << "\t isign = " << isign << "\t norm2(chi2,rb[1])=" << norm2(chi2, rb[1]) << "\t || diff || = "<< sqrt(norm2(diff_float, rb[1])) / ( Real(4*3*2*Layout::vol()) / Real(2))  << std::endl;
  }
  

  XMLFileWriter xml_out("fdiff");

  //for(int isign=-1; isign < 2; isign+=2) {
  for(int isign=+1; isign > -2; isign-=2) {
    chi=zero;
    chid=zero;

    // Apply clover op to psi into chi
    Klov32((float *)&(chi.elem(0).elem(0).elem(0).real()),	  
	   (float *)&(psi.elem(0).elem(0).elem(0).real()),
	   (float *)&(packed_gauge[0]),
	   (float *)&(clov[0]),
	   (float *)&(invclov[0]),
	   isign); 


    // Apply clover op to psi into chi
    Klov64((double *)&(chid.elem(0).elem(0).elem(0).real()),	  
	   (double *)&(psid.elem(0).elem(0).elem(0).real()),
	   (double *)&(packed_gauged[0]),
	   (double *)&(clovd[0]),
	   (double *)&(invclovd[0]),
	   isign); 
    
    // Downcast
    chi2 = chid;


    LatticeFermionF3 diff_float = chi2 - chi;

    QDPIO::cout << std::endl;
    QDPIO::cout << "\t isign = " << isign << "\t || diff || = "<< sqrt(norm2(diff_float, rb[1])) / ( Real(4*3*2*Layout::vol()) / Real(2))  << std::endl;

  }
  free(xclovd);
  free(xinvclovd);


}


