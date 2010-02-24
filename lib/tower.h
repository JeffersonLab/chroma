#ifndef __TOWER_H__
#define __TOWER_H__

#include "chromabase.h"
using namespace QDP;

namespace Chroma { 

template<typename T>
class Tower<T> { 
public: 

  //! Create empty
  Tower(int n) { 
    tower_data.resize(n);
  }
  
  // Destroy
  ~Tower() {}

  void resize(int n) {
    tower_data.resize(n);
  }

  int size(int n) const { 
    return tower_data.size();
  }

  //! get
  const T& operator[](int i) const { 
    return tower_data[i];
  }
  
  //! add item
  T& operator[](int i) {
    return tower_data[i];
  }

  //! Copy
  Tower(const Tower<T>& a) { 
    // Avoid self assignment
    if( this == &a ) return;
    
    this->resize(a.size());
    for(int i=0; i < a.size(); i++) {
      (*this)[i] = a[i];
    }
  }
      
  //! assign 
  Tower<T>& operator=(const Tower<T>&a) 
  {
    // Avoid self assignment
    if( this == &a ) {
      return *this;
    }

    this->resize(a.size());
    for(int i=0; i < a.size(); i++) {
      (*this)[i] = a[i];
    }

    return *this;
  }

  // Add: Assumes types T and T1 have += operation
  Tower<T>& operator+=(const Tower<T>& a)
  {
    for(int i=0; i < a.size(); i++) {
      (*this)[i] += a[i];
    }
    return *this;
  }

  Tower<T>& operator-=(const Tower<T>& a)
  {
    for(int i=0; i < a.size(); i++) {
      (*this)[i] -= a[i];
    }
    return *this;
  }
     
  //! Scalar multiply?
  template<typename T1>
  Tower<T>& operator*=(const OScalar<T1>& alpha)
  {
    for(int i=0; i < a.size(); i++) {
      (*this)[i] *= alpha;
    }
    return *this;
  }
 
  Tower<T>& operator*=(const Tower<T>& b)
  {

    for(int level=a.size()-1; level >=0; --level) { 
      T tmp = (*this)[0]*b[level];
      for(int i=1; i <= level; i++) { 
	T tmp2 = (*this)[i]*b[level-i];
        tmp += Choose(level,i)*tmp2;
      }
      (*this)[level] = tmp;
    }
    return (*thia);
    
  }

  private:
  //! Evaluate (m,n) binomial coefficient recursively
  unsigned int Choose(unsigned int m, unsigned int n) 
  {
    if( m==0 ) return 1;
    if( n==m ) return 1;
    return Choose(n-1,m)+Choose(n-1,m-1);
  }

  multi1d<T> tower_data;
};
 



  Tower<T> operator*(const Tower<T>& a, const Tower<T>& b) const
  {
    Tower<T> ret_val(a);
    ret_val *= b;
    return ret_val;
  }
   

  Tower<T> operator+(const Tower<T>& a, const Tower<T>& b) const
  {
    Tower<T> ret_val(a);
    ret_val += b;
    return ret_val;
  }

  Tower<T> operator-(const Tower<T>& a, const Tower<T>& b) const
  {
    Tower<T> ret_val(a);
    ret_val += b;
    return ret_val;

  }

  template<typename T1>
  Tower<T> operator*(const OScalar<T1>& a, const Tower<T>& b) const
  {
    Tower<T> ret_val(b);
    ret_val *= a;
    return ret_val;
  }

  Tower<T> adj(const Tower<T>& a) 
  {
    Tower<T> ret_val(a.size());
    for(int i=0; i < a.size(); i++) { 
      ret_val[i] = adj(a[i]);
    }
    return ret_val;
  }

  Tower<T> shift(const Tower<T>& a, int mu, int dir) 
  {
    Tower<T> ret_val(a.size());
    for(int i=0; i < a.size(); i++) { 
      ret_val[i] = shift(a[i], mu, dir);
    }
    return ret_val;
  }

  
}

#endif
