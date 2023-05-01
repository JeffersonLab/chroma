#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>

#include <cstdio>

#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

#include "chroma.h"
#include "gtest/gtest.h"
#include "chroma_gtest_env.h"

#include "util/gauge/weak_field.h"


using namespace Chroma;



class TestEnvironment : public ::testing::Environment {
public: 

	using T = LatticeFermion;
	using Q = multi1d<LatticeColorMatrix>;
	using P = multi1d<LatticeColorMatrix>;

    TestEnvironment()
    {
        const int nrow_in[4] = {2,2,4,4};
        multi1d<int> nrow(4);
        nrow = nrow_in;
        Layout::setLattSize(nrow);
        Layout::create();
    }
  


    ~TestEnvironment() {}

    Q u;
    Handle<FermState<T,P,Q> > state;

};


int main(int argc, char *argv[]) 
{
  ::testing::InitGoogleTest(&argc, argv);
  ::testing::Environment* const chroma_env = ::testing::AddGlobalTestEnvironment(new ChromaEnvironment(&argc,&argv));
  ::testing::Environment* const test_env = ::testing::AddGlobalTestEnvironment(new TestEnvironment());
  return RUN_ALL_TESTS();
}







