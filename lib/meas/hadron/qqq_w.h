// -*- C++ -*-
// $Id: qqq_w.h,v 3.1 2008-09-09 20:30:42 kostas Exp $
/*! \file
 *  \brief constructs 3 quark propagators contracted at the sink
 */

#include "chromabase.h"
#include "util/ft/sftmom.h"

#ifndef __qqq_h__
#define __qqq_h__

//! Baryon-Baryon 4-pt function building block
/* This routine is specific to Wilson fermions!
 *
 */

namespace Chroma {

  class QuarkIndex{
  public:
    int s;
    int c;
    QuarkIndex():s(0),c(0){}
    QuarkIndex(int s_,int c_):s(s_),c(c_){}

    void Zero(){s=0;c=0;}

    //this code may be slow and may need replacement
    //with simple loops... For the momment it reduces
    //typing so I keep it!
    QuarkIndex& operator++(){
      if(c<Nc-1)
	c++;
      else{
	c=0;
	s++;
      }
      return *this ;
    }

    bool NotEnd(){
      return(s<Ns);
    }
    //End of useless code!

  } ;

  class ThreeQuarks{
  private:
    multi1d <DComplex> data ;
    int size ;

  public:
    
    ThreeQuarks(){
      size  = (Nc*Ns)*(Nc*Ns)*(Nc*Ns)*Ns ;
      data.resize(size) ;
    };
 
    ~ThreeQuarks(){}
    
    // index speed in memory: c1 s1 c2 s2 c3 s3 s4
    DComplex& operator()(int s1, int c1,
			 int s2, int c2,
			 int s3, int c3,
			 int s4 ){
      return data[c1+Nc*(s1+Ns*(c2+Nc*(s2+Ns*(c3+Nc*(s3+Ns*s4)))))];
    }
    
    DComplex& operator()(const QuarkIndex& q1,
			 const QuarkIndex& q2,
			 const QuarkIndex& q3,
			 int s4 ){
      return data[q1.c+Nc*(q1.s+Ns*(q2.c+Nc*(q2.s+Ns*(q3.c+Nc*(q3.s+Ns*s4)))))];
    }
	      
    multi1d<DComplex>& handle(){
      return data;
    }

    DComplex& operator[](int i){
      return data[i];
    }

    int Size() const {return size;}
    
};

  void compute_qqq(multi2d<ThreeQuarks>& qqq, 
		   const LatticePropagator& q1,
		   const LatticePropagator& q2,
		   const LatticePropagator& q3,
		   const SftMom& phases,
		   int t0, int bc_spec
		   ) ;

  void compute_qqq(multi2d<ThreeQuarks>& qqq, const int k,
		   const LatticePropagator& q1,
		   const LatticePropagator& q2,
		   const LatticePropagator& q3,
		   const SftMom& phases,
		   int t0, int bc_spec
		   ) ;
    
  void write_qqq(QDPFileWriter& to,
		 multi2d<ThreeQuarks>& qqq, 
		 const SftMom& phases,
		 string type,
		 string sink) ;

};
#endif
