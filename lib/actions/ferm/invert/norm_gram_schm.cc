// -*- C++ -*-
// $Id: norm_gram_schm.cc,v 1.2 2007-10-24 02:38:56 edwards Exp $
/*! \file
 *  \brief Gram-Schmidt with normalization
 */

#include "actions/ferm/invert/norm_gram_schm.h"

namespace Chroma 
{

  //! Gram-Schmidt with normalization
  /*! \ingroup invert */
  template<typename T>
  void normGramSchmidt_T(multi1d<T>& vec, 
			 int f,
			 int t,
			 const Subset& sub)
  {
    for(int i(0);i<f;i++){// normalize the first vectors...
      vec[i][sub] /= Real(sqrt(norm2(vec[i],sub))) ;
    }
    if(!(t<=vec.size())){
      QDPIO::cerr<<"SimpleGramSchmidt:: f="<<f<<" t="<<t<<" vec.size()="<<vec.size()<<endl;
      QDPIO::cerr<<"SimpleGramSchmidt:: Out of bound!\n";
      exit(1);
    }
    for(int i(f);i<t;i++)
    { // now orthonormalize ther rest
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

  //
  // Wrappers
  //
  void normGramSchmidt(multi1d<LatticeFermionF>& vec, 
		       int f,
		       int t,
		       const Subset& sub)
  {
    normGramSchmidt_T<LatticeFermionF>(vec, f, t, sub);
  }

  void normGramSchmidt(multi1d<LatticeFermionD>& vec, 
		       int f,
		       int t,
		       const Subset& sub)
  {
    normGramSchmidt_T<LatticeFermionD>(vec, f, t, sub);
  }

  void normGramSchmidt(multi1d<LatticeStaggeredFermionF>& vec, 
		       int f,
		       int t,
		       const Subset& sub)
  {
    normGramSchmidt_T<LatticeStaggeredFermionF>(vec, f, t, sub);
  }

  void normGramSchmidt(multi1d<LatticeStaggeredFermionD>& vec, 
		       int f,
		       int t,
		       const Subset& sub)
  {
    normGramSchmidt_T<LatticeStaggeredFermionD>(vec, f, t, sub);
  }
  
} // End Namespace Chroma

