#include "state.h"

namespace Chroma {

  template<>
  void recursePQState(const multi1d<LatticeColorMatrixF>& p, 
		      const multi1d<LatticeColorMatrixF>& q,
		      int n,
		      multi1d<LatticeColorMatrixF>& pq)
  {
    if( p.size()!=q.size() ) { 
      QDPIO::cout << "p and q have incompatible sizes" << endl;
      QDP_abort(1);
    }
    
    LatticeColorMatrixF tmp;            
    pq.resize(p.size());
    
    for(int mu=0; mu < Nd; mu++) { 
      tmp = p[mu];
      
      for(int i=0; i < n-1; i++){ 
	pq[mu] = tmp * p[mu];
	tmp = pq[mu];
      }
      
      pq[mu] = tmp*q[mu];
    }
  }

  template<>
  void recursePQState(const multi1d<LatticeColorMatrixD>& p, 
		      const multi1d<LatticeColorMatrixD>& q,
		      int n,
		      multi1d<LatticeColorMatrixD>& pq)
  {
    if( p.size()!=q.size() ) { 
      QDPIO::cout << "p and q have incompatible sizes" << endl;
      QDP_abort(1);
    }
    
    LatticeColorMatrixD tmp;            
    pq.resize(p.size());
    
    for(int mu=0; mu < Nd; mu++) { 
      tmp = p[mu];
      
      for(int i=0; i < n-1; i++){ 
	pq[mu] = tmp * p[mu];
	tmp = pq[mu];
      }
      
      pq[mu] = tmp*q[mu];
    }
  }

}
