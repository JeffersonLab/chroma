// -*- C++ -*-
// $Id: containers.h,v 1.7 2007-10-09 18:07:14 edwards Exp $

#ifndef _INV_CONTAINERS__H
#define _INV_CONTAINERS__H

#include "chromabase.h"

namespace Chroma
{
  namespace LinAlg
  {
    //! Hold vectors
    /*! \ingroup invert */
    template<class T> class Vectors
    {
    public:
      multi1d<T> vec ;
      int N ; // number of active vectors size of vec must be larger or equal to N

      Vectors():N(0){} 
      Vectors(const multi1d<T>& v):vec(v),N(v.size()){} 
      Vectors(int size){resize(size); } 
    
      ~Vectors(){}
    
      void AddVector(const T& v,const Subset& s){
	if(N<vec.size()){
	  vec[N][s] = v ;
	  N++ ;
	}
      }

      void NormalizeAndAddVector(const T& v,const Double& inorm, 
				 const Subset& s){
	if(N<vec.size()){// inorm is the inverse of the norm
	  vec[N][s] = v  ;
	  vec[N][s] *= inorm  ;
	  N++ ;
	}
      }
      void AddOrReplaceVector(const T& v,const Subset& s){
	if(N<vec.size()){
	  vec[N][s] = v ;
	  N++ ;
	}
	else{// replace the last vector
	  vec[N-1] = v ;
	}
      }

      // This will only add as many vectors as they fit
      void AddVectors(multi1d<T>& v,const Subset& s){
	for(int i(0);i<v.size();i++)
	  AddVector(v[i],s) ;
      }

    
      void resize(int n) {N=0;vec.resize(n);} 
      int size() const { return vec.size();}
      int Nvecs() const { return N; } 
      T& operator[](int i){ return vec[i];}
    };


    //! Holds eigenvalues and eigenvectors
    /*! \ingroup invert */
    template<class T> class RitzPairs
    {
    public:
      Vectors<Double> eval;
      Vectors<T>      evec;
      int Neig;

      RitzPairs() {init(0);}
      RitzPairs(int N) {init(N);}

      void init(int N) {
	eval.resize(N);
	evec.resize(N);
	Neig = 0;
      }

      void AddVector(const Double& e, const T& v,const Subset& s){
	eval.AddVector(e,s);
	evec.AddVector(v,s);
	if (eval.N != evec.N)
	{
	  QDPIO::cerr << __func__ << ": length of value and vector arrays are not the same" << endl;
	  QDP_abort(1);
	}
	Neig = evec.N;
      }

      // This will only add as many vectors as they fit
      void AddVectors(const multi1d<T>& e,const multi1d<T>& v,const Subset& s){
	for(int i(0);i<e.size();i++)
	  AddVector(e[i],s) ;
	for(int i(0);i<v.size();i++)
	  AddVector(v[i],s) ;
	if (eval.N != evec.N)
	{
	  QDPIO::cerr << __func__ << ": length of value and vector arrays are not the same" << endl;
	  QDP_abort(1);
	}
      }
    };

    //! This is a square matrix
    /*! \ingroup invert */
    template<class T> class Matrix
    {
    public:
      multi2d<T> mat ;
      int N ;  // active size 

      Matrix():N(0){} 
      Matrix(const multi2d<T>& v):mat(v),N(v.size1())
	{
	  if(v.size1() != v.size2())
	    QDPIO::cerr<<"WARNING!!! Matrix should be square!! CHECK YOUR CODE!\n";
	}

      Matrix(int size) {resize(size);}
    
      ~Matrix(){}
    
      void resize(int size) {N=0; mat.resize(size,size);}
      int size() const { return mat.size1();}
      int ld() const { return mat.size1();}
      int Nmat() const { return N; } 
      T& operator()(int i,int j){ return mat(i,j);}
    };
  
  }

}

#endif 
