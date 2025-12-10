#ifndef CPP_DSLASH_QDP_PACKER_H
#define CPP_DSLASH_QDP_PACKER_H


#include "qdp.h"
using namespace QDP;

namespace CPlusPlusWilsonDslash { 

  typedef PColorMatrix<RComplex<REAL32>, 3> PrimitiveSU3MatrixF;
  typedef PColorMatrix<RComplex<REAL64>, 3> PrimitiveSU3MatrixD;
  typedef PSpinVector< PColorVector< RComplex<REAL32>, 3>, 4> PrimitiveSpinorF;
  typedef PSpinVector< PColorVector< RComplex<REAL64>, 3>, 4> PrimitiveSpinorD;
  

  void qdp_pack_gauge(const multi1d<LatticeColorMatrixF>&_u, PrimitiveSU3MatrixF* u_tmp);
  
  void qdp_pack_gauge_3d(const multi1d<LatticeColorMatrixF>&_u, multi1d<PrimitiveSU3MatrixF>& u_tmp);

  void qdp_pack_gauge(const multi1d<LatticeColorMatrixD>&_u, PrimitiveSU3MatrixD* u_tmp);

  void qdp_pack_gauge_3d(const multi1d<LatticeColorMatrixD>&_u, multi1d<PrimitiveSU3MatrixD>& u_tmp);

};

#endif
