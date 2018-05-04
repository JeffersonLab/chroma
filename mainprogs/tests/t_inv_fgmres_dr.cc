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

#include "actions/ferm/invert/syssolver_fgmres_dr_params.h"

using namespace Chroma;



std::string fermact_xml =   
  "<?xml version='1.0'?>                              \
   <Param>					      \
   <FermionAction>                                    \
     <FermAct>CLOVER</FermAct>                        \
     <Mass>0.1</Mass>				      \
     <clovCoeff>1</clovCoeff>			      \
     <AnisoParam>				      \
       <anisoP>false</anisoP>			      \
       <t_dir>3</t_dir>				      \
       <xi_0>1</xi_0>				      \
       <nu>1</nu>				      \
     </AnisoParam>				      \
    <FermState>					      \
       <Name>STOUT_FERM_STATE</Name>		      \
       <rho>0.14</rho>				      \
       <n_smear>2</n_smear>			      \
       <orthog_dir>3</orthog_dir>		      \
       <FermionBC>				      \
         <FermBC>SIMPLE_FERMBC</FermBC>		      \
         <boundary>1 1 1 -1</boundary>		      \
       </FermionBC>				      \
     </FermState>				      \
   </FermionAction>				      \
  </Param>";
  	 

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







