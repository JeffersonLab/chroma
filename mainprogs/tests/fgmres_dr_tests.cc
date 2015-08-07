#include "gtest/gtest.h"
#include "chromabase.h"
#include "actions/ferm/invert/syssolver_fgmres_dr_params.h"

using namespace std;
using namespace Chroma;

TEST(FGMRESDRParams, canConstructDefault)
{
  SysSolverFGMRESDRParams p;
  ASSERT_EQ(p.NKrylov, 0);
  ASSERT_EQ(p.NDefl, 0);
  ASSERT_EQ(p.MaxIter, 0);
  ASSERT_EQ(p.PrecondParams.path, "/PrecondParams");
}

std::string xml_for_param = 
  "<?xml version='1.0'?> \
   <Params> \
     <InvertParam>				      \
     <invType>FGMRESDR_INVERTER</invType>	      \
     <RsdTarget>1.0e-7</RsdTarget>		      \
     <NKrylov>30</NKrylov>			      \
     <NDefl>6</NDefl>				      \
     <MaxIter>130</MaxIter>			      \
     <PrecondParams>				      \
       <invType>MR_INVERTER</invType>		      \
       <RsdMR>0.1</RsdMR>			      \
       <MaxMR>10</MaxMR>			      \
     </PrecondParams>				      \
   </InvertParam>				      \
  </Params>";

TEST(FGMRESDRParams, canReadXML)
{
  std::istringstream input(xml_for_param);
  XMLReader xml_in(input);
  SysSolverFGMRESDRParams p( xml_in, "/Params/InvertParam" );
  ASSERT_EQ(p.NKrylov, 30);
  ASSERT_EQ(p.NDefl, 6);
  ASSERT_EQ(p.MaxIter, 130);
  ASSERT_EQ(p.PrecondParams.path, "/PrecondParams");
  ASSERT_EQ(p.PrecondParams.id, "MR_INVERTER");
}
