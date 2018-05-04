// -*- C++ -*-

#ifndef _OCTAVE_DEBUG_H
#define _OCTAVE_DEBUG_H

#include <sstream>
#include "chromabase.h"



// NEEDS A LOT OF CLEAN UP
namespace Octave 
{
  int count ;
  std::string tag(const std::string& prefix){
    std::stringstream t ;
    t<<prefix<<count ;
    count++;
    return t.str() ;
  }

  void PrintClear( const std::string& fname){
    count=0;
    ofstream foo(fname.c_str(),ios::trunc);
    foo<<"# BEGIN FILE:"<<fname<<"\n" ;
    foo.close();
  }


  void PrintOut(const multi2d<DComplex>& H, const int& N, 
                      const std::string& tag,
                      const std::string& fname)
  {
    ofstream foo(fname.c_str(),ios::app);
    foo<<tag<<" =["<<std::endl; 
    for(int i(0);i<N;i++){
      for(int j(0);j<N;j++)
        foo<<real(H(i,j))<<"+i*"<<imag(H(i,j))<<" ";
      foo<<std::endl ;
    }
    foo<<"];"<<std::endl;
    foo.close();
  }


  void PrintOut(const multi1d<DComplex>& H,
                      const int& N,
                      const std::string& tag, 
                      const std::string& fname)
  {
    ofstream foo(fname.c_str(),ios::app);
    foo<<tag<<" =["<<std::endl; 
    for(int i(0);i<N;i++){
      foo<<real(H[i])<<"+i*"<<imag(H[i])<<" ";
      foo<<std::endl ;
    }
    foo<<"];"<<std::endl;
    foo.close();
  }



} ; // End Namespace Octave

#endif
