#include "gtest/gtest.h"
#include "chromabase.h"
#include "mgparams_tests.h"
#include "actions/ferm/invert/quda_solvers/quda_eig_params.h"
#include "actions/ferm/invert/quda_solvers/syssolver_quda_multigrid_clover_params.h"
//using namespace std;
using namespace Chroma;

class MGParamTests : public ::testing::Test {
public:
  MGParamTests()
  {
  }

  void TearDown()
  {
   
  }

  // Virtual destructor
  virtual
  ~MGParamTests() {}

};

TEST_F(MGParamTests, BasicEigRead)
{
    QudaEigParam p = newQudaEigParam();
    std::istringstream input_xml(MGParamTestsEnv::eig_xml);
    std::string outbuf;
    XMLReader xml_in(input_xml);
    read(xml_in, "/Param/EigParams", p);
    XMLBufferWriter xml_out(outbuf);
    write(xml_out, "EigParams", p);
    QDPIO::cout << xml_out.str() << "\n";
}

TEST_F(MGParamTests, ReadMG)
{ 
  std::istringstream stream_in1(MGParamTestsEnv::inv_param_quda_multigrid_asymm_xml);
  std::istringstream stream_in2(MGParamTestsEnv::inv_param_quda_multigrid_asymm_eig_xml);
  XMLReader without_ev(stream_in1);
  XMLReader with_ev(stream_in2);
  SysSolverQUDAMULTIGRIDCloverParams p1(without_ev,"/Param/InvertParam");
  SysSolverQUDAMULTIGRIDCloverParams p2(with_ev, "/Param/InvertParam");

  QDPIO::cout << "Found delfation params in p1: " 
    << (p1.MULTIGRIDParams->got_mg_eig_params ? " yes " : " no ") << "\n";

  QDPIO::cout << "Found delfation params in p2: " 
    << (p2.MULTIGRIDParams->got_mg_eig_params ? " yes " : " no ") << "\n";

  XMLBufferWriter xml_out;
  write(xml_out, "DeflationEigParam", p2.MULTIGRIDParams->mg_eig_params);
  QDPIO::cout << "Deflation Params in p2 are:\n";
  QDPIO::cout << xml_out.str() << "\n";

}