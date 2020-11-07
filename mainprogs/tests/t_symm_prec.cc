#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>

#include <cstdio>

#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

#include "chroma.h"
#include "handle.h"
#include "eoprec_linop.h"
#include "seoprec_logdet_linop.h"
#include "eoprec_wilstype_fermact_w.h"
#include "seoprec_wilstype_fermact_w.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "gtest/gtest.h"
#include "chroma_gtest_env.h"

#include "util/gauge/weak_field.h"


using namespace Chroma;



class TestEnvironment : public ::testing::Environment {
public: 

	using T = LatticeFermion;
	using Q = multi1d<LatticeColorMatrix>;
	using P = multi1d<LatticeColorMatrix>;

	using S_asymm_T = EvenOddPrecWilsonTypeFermAct<T,P,Q>;
	using S_symm_T =  SymEvenOddPrecWilsonTypeFermAct<T,P,Q>;
	using LinOpSymm_T = SymEvenOddPrecLinearOperator<T,P,Q>;
	using LinOpAsymm_T = EvenOddPrecLinearOperator<T,P,Q>;

  TestEnvironment()
  {
    const int nrow_in[4] = {4,4,4,8};
    multi1d<int> nrow(4);
    nrow = nrow_in;
    Layout::setLattSize(nrow);
    Layout::create();
  }
  


  ~TestEnvironment() {
  }

  Q u;
  Handle<S_symm_T> S_symm;
  Handle<S_asymm_T> S_asymm;
  Handle<FermState<T,P,Q> > state;

};


int main(int argc, char *argv[]) 
{
  ::testing::InitGoogleTest(&argc, argv);
  ::testing::Environment* const chroma_env = ::testing::AddGlobalTestEnvironment(new ChromaEnvironment(&argc,&argv));
  ::testing::Environment* const test_env = ::testing::AddGlobalTestEnvironment(new TestEnvironment());
  return RUN_ALL_TESTS();
}







