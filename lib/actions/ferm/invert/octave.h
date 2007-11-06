// -*- C++ -*-
// $Id: octave.h,v 1.1 2007-11-06 22:04:15 kostas Exp $

#ifndef _OCTAVE_DEBUG_H
#define _OCTAVE_DEBUG_H

#include <sstream>
#include "chromabase.h"



// NEEDS A LOT OF CLEAN UP
namespace Octave 
{
  int count ;
  string tag(const string& prefix){
    stringstream t ;
    t<<prefix<<count ;
    count++;
    return t.str() ;
  }

  void PrintClear( const string& fname){
    count=0;
    ofstream foo(fname.c_str(),ios::trunc);
    foo<<"# BEGIN FILE:"<<fname<<"\n" ;
    foo.close();
  }


  void PrintOut(const multi2d<DComplex>& H, const int& N, 
                      const string& tag,
                      const string& fname)
  {
    ofstream foo(fname.c_str(),ios::app);
    foo<<tag<<" =["<<endl; 
    for(int i(0);i<N;i++){
      for(int j(0);j<N;j++)
        foo<<real(H(i,j))<<"+i*"<<imag(H(i,j))<<" ";
      foo<<endl ;
    }
    foo<<"];"<<endl;
    foo.close();
  }


  void PrintOut(const multi1d<DComplex>& H,
                      const int& N,
                      const string& tag, 
                      const string& fname)
  {
    ofstream foo(fname.c_str(),ios::app);
    foo<<tag<<" =["<<endl; 
    for(int i(0);i<N;i++){
      foo<<real(H[i])<<"+i*"<<imag(H[i])<<" ";
      foo<<endl ;
    }
    foo<<"];"<<endl;
    foo.close();
  }



} ; // End Namespace Octave

#endif
