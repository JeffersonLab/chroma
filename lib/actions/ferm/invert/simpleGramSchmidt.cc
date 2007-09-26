// -*- C++ -*-
// $Id: simpleGramSchmidt.cc,v 1.1 2007-09-26 02:46:00 kostas Exp $

#include "simpleGramSchmidt.h"

namespace Chroma 
{
  typedef LatticeFermion T ; // save sometyping 

  void SimpleGramSchmidt(multi1d<T>& vec, 
			 const int f,
			 const int t,
			 const OrderedSubset& sub){

    for(int i(0);i<f;i++){// normalize the first vectors...
      vec[i][sub] /= Real(sqrt(norm2(vec[i],sub))) ;
    }
    if(!(t<=vec.size())){
      QDPIO::cerr<<"SimpleGramSchmidt:: f="<<f<<" t="<<t<<" vec.size()="<<vec.size()<<endl;
      QDPIO::cerr<<"SimpleGramSchmidt:: Out of bound!\n";
      exit(1);
    }
    for(int i(f);i<t;i++){ // now orthonormalize ther rest
      for(int k(0);k<i;k++){
	DComplex dcc = innerProduct(vec[k], vec[i], sub);
	Complex cc = dcc ;
	//cout<<"GramS: "<<cc<<" "<<dcc<<endl;
	vec[i][sub] = vec[i]  - cc*vec[k] ;
      }
      Double in = 1.0/sqrt(norm2(vec[i],sub)) ;
      vec[i][sub] *= Real(in); 
    }
  }

} // namespace Chroma




