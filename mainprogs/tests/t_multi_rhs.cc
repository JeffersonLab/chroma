#include "chroma.h"
#include "gtest/gtest.h"
#include "chroma_gtest_env.h"

using namespace Chroma;

class TestEnvironment : public ::testing::Environment {
public: 

  
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


};


int main(int argc, char *argv[]) 
{
  ::testing::InitGoogleTest(&argc, argv);
  ::testing::Environment* const chroma_env = ::testing::AddGlobalTestEnvironment(new ChromaEnvironment(&argc,&argv));
  ::testing::Environment* const test_env = ::testing::AddGlobalTestEnvironment(new TestEnvironment());
  return RUN_ALL_TESTS();
}







