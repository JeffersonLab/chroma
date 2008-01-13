// -*- C++ -*-
// $Id: norm_gram_schm.cc,v 1.4 2008-01-13 21:08:13 edwards Exp $
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
    // this need not to be done... we assume 0..f vectors are normal already
    for(int i(0);i<f;i++){// normalize the first vectors...
      vec[i][sub] /= Real(sqrt(norm2(vec[i],sub)));
    }
    if(!(t<=vec.size())){
      QDPIO::cerr << __func__ << ": f="<<f<<" t="<<t<<" vec.size()="<<vec.size()<<endl;
      QDPIO::cerr << __func__ << ": Out of bound!\n";
      exit(1);
    }
    for(int i(f);i<t;i++)
    { // now orthonormalize ther rest
      for(int k(0);k<i;k++){
	DComplex dcc = innerProduct(vec[k], vec[i], sub);
	Complex cc = dcc;
	//QDPIO::cout<<"GramS: "<<i<<" "<<k<<" "<<cc<<" "<<dcc<<endl;
	vec[i][sub] = vec[i]  - cc*vec[k];
      }
      Double in = 1.0/sqrt(norm2(vec[i],sub));
      vec[i][sub] *= Real(in); 
    }
  }


  //! Gram-Schmidt with normalization
  /*! \ingroup invert */
  template<typename T>
  void normGramSchmidtArray_T(multi2d<T>& vec, 
			      int f,
			      int t,
			      const Subset& sub)
  {
    // this need not to be done... we assume 0..f vectors are normal already
    for(int i(0);i<f;i++){// normalize the first vectors...
      Real in = Real(1) / Real(sqrt(norm2(vec[i],sub)));
      for(int s=0; s < vec.size1(); ++s)
	vec[i][s][sub] *= in;
    }
    if(!(t<=vec.size2())){
      QDPIO::cerr << __func__ << ": f="<<f<<" t="<<t<<" vec.size2()="<<vec.size2()<<endl;
      QDPIO::cerr << __func__ << ": Out of bound!\n";
      exit(1);
    }
    for(int i(f);i<t;i++)
    { // now orthonormalize ther rest
      for(int k(0);k<i;k++){
	DComplex dcc = innerProduct(vec[k], vec[i], sub);
	Complex cc = dcc;
	//QDPIO::cout<<"GramS: "<<i<<" "<<k<<" "<<cc<<" "<<dcc<<endl;
	for(int s=0; s < vec.size1(); ++s)
	  vec[i][s][sub] = vec[i][s]  - cc*vec[k][s];
      }
      Real in = Real(1.0) / sqrt(norm2(vec[i],sub));
      for(int s=0; s < vec.size1(); ++s)
	vec[i][s][sub] *= Real(in); 
    }
  }

  //
  // Wrappers
  //
  // 4D versions
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
  

  // 5D versions
  void normGramSchmidt(multi2d<LatticeFermionF>& vec, 
		       int f,
		       int t,
		       const Subset& sub)
  {
    normGramSchmidtArray_T(vec, f, t, sub);
  }

  void normGramSchmidt(multi2d<LatticeFermionD>& vec, 
		       int f,
		       int t,
		       const Subset& sub)
  {
    normGramSchmidtArray_T(vec, f, t, sub);
  }

} // End Namespace Chroma

