// -*- C++ -*-
// $Id: staggered_operators_s.cc,v 1.4 2006-11-19 05:01:34 kostas Exp $
/*! \file
 *  \brief Staggered  operators
 *
 */

#include "chromabase.h"
#include "util/ferm/staggered_operators_s.h"

namespace Chroma 
{

   
  class AntiSymmetricTensor4D {
  private:
    class Datum{
    public:
      multi1d<int> d;
      Real sign ;  
      void init(int i,int j, int k, int l, Real s){
	d.resize(4) ; sign=s ;
	d[0]=i;d[1]=j;d[2]=k;d[3]=l ;
      }
      Datum(){d.resize(4);sign=1.0;}
      Datum(int i,int j, int k, int l, Real s){
	init(i,j,k,l,s);
      }
      ~Datum(){}
    } ;
    
    static Datum eps[24] ;
    
    void InitializeTable(){
      eps[0].init(0,1,2,3,+1.0);
      eps[1].init(0,3,1,2,+1.0);
      eps[2].init(0,2,3,1,+1.0);
      eps[3].init(0,3,2,1,-1.0);
      eps[4].init(0,1,3,2,-1.0);
      eps[5].init(0,2,1,3,-1.0);
      
      eps[6].init(1,0,2,3,-1.0);
      eps[7].init(1,3,0,2,-1.0);
      eps[8].init(1,2,3,0,-1.0);
      eps[9].init(1,3,2,0,+1.0);
      eps[10].init(1,0,3,2,+1.0);
      eps[11].init(1,2,0,3,+1.0);
      
      eps[12].init(2,1,0,3,-1.0);
      eps[13].init(2,3,1,0,-1.0);
      eps[14].init(2,0,3,1,-1.0);
      eps[15].init(2,3,0,1,+1.0);
      eps[16].init(2,1,3,0,+1.0);
      eps[17].init(2,0,1,3,+1.0);
      
      eps[18].init(3,1,2,0,-1.0);
      eps[19].init(3,0,1,2,-1.0);
      eps[20].init(3,2,0,1,-1.0);
      eps[21].init(3,0,2,1,+1.0);
      eps[22].init(3,1,0,2,+1.0);
      eps[23].init(3,2,1,0,+1.0);
    }
    
  public:
    AntiSymmetricTensor4D(){
      InitializeTable();
    }
    
    ~AntiSymmetricTensor4D(){} 
    
    Real sign(int k){return eps[k].sign ;} 
    multi1d<int> indx(int k){return eps[k].d;}
    
  } ;
  


  static AntiSymmetricTensor4D eps ;

  
  typedef LatticeStaggeredPropagator  T ;
  typedef multi1d<LatticeColorMatrix> G ;
  
  void StaggeredZeta(LatticeStaggeredPropagator& dest,int mu){
    LatticeReal sign = 1.0 ;
    
    for(int c(mu+1);c<Nd;c++){
      where(QDP::Layout::latticeCoordinate(c) & 1 == 1, -1.0*sign, sign);
    }
    dest *= sign ;
  }
  
  void StaggeredEta(LatticeStaggeredPropagator& dest,int mu){
    LatticeReal sign = 1.0 ;
    
    for(int c(0);c<mu;c++){
      where(QDP::Layout::latticeCoordinate(c) & 1 == 1, -1.0*sign, sign);
    }

    dest *= sign ;
  }

  void SymShift(T& dest,const T& src,const G& u,const int mu){
    dest = u[mu]*shift(src,FORWARD,mu) + shift(adj(u[mu])*src,BACKWARD,mu) ;
  }
  

  void EtaShift(T& dest, const T& src, const  G& u, const multi1d<int>& mu){
    T tmp(src);
    for(int c(0) ; c<mu.size();c++){
      SymShift(dest, tmp, u, mu[c]);
      StaggeredEta(dest,mu[c]);
      tmp = dest ;
    }
  }

  void EtaShift(T& dest, const T& src, const  G& u, const int mu){
    SymShift(dest, src, u, mu);
    StaggeredEta(dest,mu);
  }

  void ZetaShift(T& dest, const T& src, const  G& u, const multi1d<int>& mu){
    T tmp(src);
    for(int c(0) ; c<mu.size();c++){
      SymShift(dest, tmp, u, mu[c]);
      StaggeredZeta(dest,mu[c]);
      tmp = dest ;
    }
  }

  void ZetaShift(T& dest, const T& src, const  G& u, const int mu){
    SymShift(dest, src, u, mu);
    StaggeredZeta(dest,mu);
  }

  void SpinScalar(T& dest, const T& src, const G& u ){
    dest = src ;
  }
  
  void SpinVector(T& dest, const T& src, const G& u,const int mu){
    EtaShift(dest,src,u,mu) ;
  }

  void SpinTensor(T& dest, const T& src, const G& u,const int mu,const int nu){
    T tmp ;

    multi1d<int> d(2) ;
    d[0]=mu;d[1]=nu;
    EtaShift(dest,src,u,d);
    d[0]=nu;d[1]=mu;
    EtaShift(tmp,src,u,d);
    
    dest -= tmp ;
    dest *= 0.5 ;
  }
  void SpinAxialVector (T& dest, const T& src, const G& u,const int mu){
    T tmp ;
    dest = zero ;
    for(int p(0);p<24;p++)
      if(eps.indx(p)[0]==mu){
	EtaShift(tmp,src,u,eps.indx(p));
	dest += eps.sign(p)/6.0*tmp ;
      } 
  }

  void SpinPseudoScalar(T& dest, const T& src, const G& u){
    T tmp ;
    dest = zero ;
    for(int p(0);p<24;p++){
      EtaShift(tmp,src,u,eps.indx(p));
      dest += eps.sign(p)/24.0*tmp ;
    } 
  }

  void FlavorScalar(T& dest,const T& src,const G& u ){
    dest = src ;
  }
  
  void FlavorVector(T& dest,const T& src,const G& u,const int mu){
    ZetaShift(dest,src,u,mu) ;
  }
  
  void FlavorTensor(T& dest,const T& src,const G& u,const int mu,const int nu){
    T tmp ;

    multi1d<int> d(2) ;
    d[0]=mu;d[1]=nu;
    ZetaShift(dest,src,u,d);
    d[0]=nu;d[1]=mu;
    ZetaShift(tmp,src,u,d);
    
    dest -= tmp ;
    dest *= 0.5 ;
  }
  
  void FlavorAxialVector (T& dest,const T& src,const G& u,const int mu){
    T tmp ;
    dest = zero ;
    for(int p(0);p<24;p++)
      if(eps.indx(p)[0]==mu){
	ZetaShift(tmp,src,u,eps.indx(p));
	dest += eps.sign(p)/6.0*tmp ;
      }
  }

  void FlavorPseudoScalar(T& dest,const T& src,const G& u){
    T tmp ;
    dest = zero ;
    for(int p(0);p<24;p++){
      ZetaShift(tmp,src,u,eps.indx(p));
      dest += eps.sign(p)/24.0*tmp ;
    } 
  }

}  // end namespace Chroma

