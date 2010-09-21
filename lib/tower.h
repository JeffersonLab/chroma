// -*- C++ -*-

#ifndef __TOWER_H__
#define __TOWER_H__

#include "chromabase.h"
using namespace QDP;

namespace Chroma { 

  // Forward Declarationxsxs
  template<typename T>
  class SubTower;

template<typename T>
class Tower { 
public: 

  //! Base
  Tower() {}

  //! Create empty
  Tower(int n) { 
    tower_data.resize(n);
  }
  
  // Destroy
  ~Tower() {}

  void resize(int n) {
    tower_data.resize(n);
  }

  int size() const { 
    return tower_data.size();
  }

  int getHeight() const { 
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

  SubTower<T> operator[](const Subset& s) { 
    SubTower<T> ret_val(*this, s);
    return ret_val;
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


  //! Unary zero:
  Tower(int size_, const QDP::Zero& z)
  {
    tower_data.resize(size_);
    for(int i=0; i < size_; i++) { 
      tower_data[i] = zero;
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


  //! assign zero
  Tower<T>& operator=(const QDP::Zero& z)
  {
    for(int i=0; i < size(); i++) { 
      tower_data[i] = zero;
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
    for(int i=0; i < this->size(); i++) {
      (*this)[i] *= alpha;
    }
    return *this;
  }
 
  Tower<T>& operator*=(const Tower<T>& b)
  {

    for(int level=b.size()-1; level >=0; --level) { 
      T tmp = (*this)[0]*b[level];
      for(int i=1; i <= level; i++) { 
	T tmp2 = (*this)[i]*b[level-i];
        tmp += Choose(level,i)*tmp2;
      }
      (*this)[level] = tmp;
    }
    return (*this);
    
  }

  //! Evaluate (m,n) binomial coefficient recursively
  int Choose(int m, int n) const
  {
    if ( m==0 && n==0 ) return 1;  // Pinnacle

    if ( m==0 ) return 0; 

    return Choose(m-1,n-1)+Choose(m-1,n);
  }

private:


  multi1d<T> tower_data;
};
 


  template<typename T>
  Tower<T> operator*(const Tower<T>& a, const Tower<T>& b) 
  {
    Tower<T> ret_val(a);
    ret_val *= b;
    return ret_val;
  }
   
  template<typename T>
  Tower<T> operator+(const Tower<T>& a, const Tower<T>& b)
  {
    Tower<T> ret_val(a);
    ret_val += b;
    return ret_val;
  }

  template<typename T>
  Tower<T> operator-(const Tower<T>& a, const Tower<T>& b)
  {
    Tower<T> ret_val(a);
    ret_val -= b;
    return ret_val;

  }

  template<typename T, typename T1>
  Tower<T> operator*(const OScalar<T1>& a, const Tower<T>& b)
  {
    Tower<T> ret_val(b);
    ret_val *= a;
    return ret_val;
  }

  template<typename T>
  Tower<T> adj(const Tower<T>& a) 
  {
    Tower<T> ret_val(a.size());
    for(int i=0; i < a.size(); i++) { 
      ret_val[i] = adj(a[i]);
    }
    return ret_val;
  }

  template<typename T>
  Tower<T> shiftTower(const Tower<T>& a, int mu, int dir) 
  {
    Tower<T> ret_val(a.size());
    for(int i=0; i < a.size(); i++) { 
      ret_val[i] = QDP::shift(a[i], mu, dir);
    }
    return ret_val;
  }

  template<typename T>
  void random(Tower<T>& t) 
  {
    for(int i=0; i < t.size(); i++) {
      random(t[i]);
    }
  }

  template<typename T>
  void gaussian(Tower<T>& t) 
  {
    for(int i=0; i < t.size(); i++) {
      gaussian(t[i]);
    }
  }

  template<typename T>
  void print(Tower<T>& t) { 
    for(int i=0; i < t.size(); i++) {
      cout << "Level["<<i<<" ]"<< t[i] << endl;
    }
  }

  template<typename T, typename U>
  void spinTraceOuterProduct(Tower<U>& result, 
			     const Tower<T>& X, 
			     const Tower<T>& Y,
			     const Subset& s)
  {
    // Zero result
    result = zero;
    
    int N=X.size();
    
    // Level 0 -- Always
    result[0][s] = traceSpin(outerProduct(X[0],Y[0]));
    for(int i=1; i < N; i++) {
      result[i][s] = traceSpin(outerProduct(X[i], Y[0]));
      for(int j=1; j <= i; j++) { 
	U tmp;
	tmp[s] = traceSpin(outerProduct(X[i-j],Y[j]));
	result[i][s] += Real(Y.Choose(i,j))*tmp;
      }
    }
  }


  // Sub Tower -- towers with subsets
  template<typename T>
  class SubTower { 
  public: 

    // SubTower just holds a reference to a Tower and a subset
    SubTower(Tower<T>& t_, const Subset& s_) : t(t_), s(s_) {} 

    // Destroy
    ~SubTower() {}

    //! Copy
    SubTower(const SubTower<T>& a) : t(a.t), s(a.s) {     }

    //! assign 
    SubTower<T>& operator=(const Tower<T>&a) 
    {
      
      t.resize(a.size());
      for(int i=0; i < a.size(); i++) {
	// Assign Under Subset
	(t[i])[s] = a[i];
      }
      
      return *this;
    }
    
   //! assign 
    SubTower<T>& operator=(const SubTower<T>&a) 
    {
      
      t.resize(a.t.size());
      for(int i=0; i < a.t.size(); i++) {
	// Assign Under Subset
	(t[i])[s] = a.t[i];
      }
      
      return *this;
    }
    
#if 0
    // This is dangerous... we don't yet know the height of the tower necessarily.
    //! assign zero
    SubTower<T>& operator=(const QDP::Zero& z)
    {
      for(int i=0; i < t.size(); i++) { 
	(t[i])[s] = zero;
      }
      return *this;
    }
#endif
    // Add: Assumes types T and T1 have += operation
    SubTower<T>& operator+=(const Tower<T>& a)
    {
      for(int i=0; i < a.size(); i++) {
	(t[i])[s] += a[i];
      }
      return *this;
    }

    SubTower<T>& operator-=(const Tower<T>& a)
    {
      for(int i=0; i < a.size(); i++) {
	(t[i])[s] -= a[i];
      }
      return *this;
    }
    
    //! Scalar multiply?
    template<typename T1>
    SubTower<T>& operator*=(const OScalar<T1>& alpha)
    {
      for(int i=0; i < this->size(); i++) {
	(t[i])[s] *= alpha;
      }
      return *this;
    }
 
    SubTower<T>& operator*=(const Tower<T>& b)
    {
      
      for(int level=b.size()-1; level >=0; --level) { 
	T tmp;
	tmp[s] = t[0]*b[level];
	for(int i=1; i <= level; i++) { 
	  T tmp2;
	  tmp2[s] =t[i]*b[level-i];
	  tmp[s] += t.Choose(level,i)*tmp2;
	}
	(t[level])[s] = tmp;
      }
      return (*this);
      
    }

  private:
    Tower<T>& t;
    const Subset& s;
  };
 
}


#endif
