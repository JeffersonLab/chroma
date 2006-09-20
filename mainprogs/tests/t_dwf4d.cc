// $Id: t_dwf4d.cc,v 3.2 2006-09-20 20:28:06 edwards Exp $

#include <iostream>
#include <cstdio>

#include "chroma.h"

#include "qdp_util.h"

using namespace Chroma;


//! To insure linking of code, place the registered code flags here
/*! This is the bit of code that dictates what fermacts are in use */
bool linkage_hack()
{
  bool foo = true;

  // All actions
  foo &= WilsonTypeFermActsEnv::registerAll();

  return foo;
}


/*
 * Input 
 */
// Parameters which must be determined from the XML input
// and written to the XML output
struct Param_t
{
  multi1d<int> nrow;		// Lattice dimension
  GroupXML_t   invParam;   // Inverter parameters
};

//! Mega-structure of all input
struct Test_input_t
{
  Param_t          param;
  Cfg_t            cfg;
  string           action_5d;
  string           action_4d;
};

//! Parameters for running code
void read(XMLReader& xml, const string& path, Param_t& param)
{
  XMLReader paramtop(xml, path);
  read(paramtop, "nrow", param.nrow);
  param.invParam = readXMLGroup(paramtop, "/InvertParam", "invType");
}


// Reader for input parameters
void read(XMLReader& xml, const string& path, Test_input_t& input)
{
  XMLReader inputtop(xml, path);

  // Read all the input groups
  try
  {
    // Read program parameters
    read(inputtop, "Param", input.param);

    // Read in the gauge configuration info
    read(inputtop, "Cfg", input.cfg);

    // Read 4D action stuff
    {
      XMLReader xml_action(inputtop, "Action4D");
      std::ostringstream os;
      xml_action.print(os);
      input.action_4d = os.str();
    }

    // Read 5D action stuff
    {
      XMLReader xml_action(inputtop, "Action5D");
      std::ostringstream os;
      xml_action.print(os);
      input.action_5d = os.str();
    }
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error reading test data: " << e << endl;
    throw;
  }
}




int main(int argc, char **argv)
{
  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);

  QDPIO::cout << "linkage=" << linkage_hack() << endl;

  // Input parameter structure
  Test_input_t  input;

  // Instantiate xml reader for DATA
  XMLReader xml_in(Chroma::getXMLInputFileName());

  // Read data
  read(xml_in, "/t_dwf4d", input);

  // Specify lattice size, shape, etc.
  Layout::setLattSize(input.param.nrow);
  Layout::create();

  QDPIO::cout << "t_dwf4d" << endl;

  // Read in the configuration along with relevant information.
  multi1d<LatticeColorMatrix> u(Nd);
  XMLReader gauge_file_xml, gauge_xml;

  // Start up the gauge field
  gaugeStartup(gauge_file_xml, gauge_xml, u, input.cfg);

  // Instantiate XML writer for XMLDAT
//  XMLFileWriter xml_out("t_dwf4d.xml");
  XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();
  push(xml_out, "t_dwf4d");

  proginfo(xml_out);    // Print out basic program info

  // Write out the input
  write(xml_out, "Input", xml_in);

  // Write out the config header
  write(xml_out, "Config_info", gauge_xml);

  // Check if the gauge field configuration is unitarized
  unitarityCheck(u);

  // Calculate some gauge invariant observables just for info.
  MesPlq(xml_out, "Observables", u);
  xml_out.flush();

  //
  // Initialize fermion action
  //
  std::istringstream  xml_s_5d(input.action_5d);
  std::istringstream  xml_s_4d(input.action_4d);
  XMLReader  fermacttop_5d(xml_s_5d);
  XMLReader  fermacttop_4d(xml_s_4d);
  const string fermact_path_5d = "/Action5D/FermionAction";
  const string fermact_path_4d = "/Action4D/FermionAction";
  string fermact_5d;
  string fermact_4d;

  try
  {
    read(fermacttop_5d, fermact_path_5d + "/FermAct", fermact_5d);
    read(fermacttop_4d, fermact_path_4d + "/FermAct", fermact_4d);
  }
  catch (const std::string& e) 
  {
    QDPIO::cerr << "Error reading fermact: " << e << endl;
    throw;
  }
  catch (const char* e) 
  {
    QDPIO::cerr << "Error reading fermact: " << e << endl;
    throw;
  }

  QDPIO::cout << "FermAct5D = " << fermact_5d << endl;
  QDPIO::cout << "FermAct4D = " << fermact_4d << endl;


  // Deal with auxiliary (and polymorphic) state information
  // eigenvectors, eigenvalues etc. The XML for this should be
  // stored as a string called "stateInfo" in the param struct.

  try {

    // Make a reader for the stateInfo
    const string state_info_path_5d = "/Action5D/StateInfo";
    const string state_info_path_4d = "/Action4D/StateInfo";
    XMLReader state_info_xml_5d(fermacttop_5d,state_info_path_5d);
    XMLReader state_info_xml_4d(fermacttop_4d,state_info_path_4d);

    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    // DWF-like 5D Wilson-Type stuff 
    bool success = false; 

    QDPIO::cerr << "create dwf = " << fermact_5d << endl;

    Handle< WilsonTypeFermAct5D<T,P,Q> >
      S_f_5d(TheWilsonTypeFermAct5DFactory::Instance().createObject(fermact_5d,
								    fermacttop_5d,
								    fermact_path_5d));
  
    Handle< FermState<T,P,Q> > state_5d(S_f_5d->createState(u,
							    state_info_xml_5d,
							    state_info_path_5d));
  
    // Overlap-like stuff
    QDPIO::cerr << "create overlap" << endl;
    Handle< WilsonTypeFermAct<T,P,Q> >
      S_f_4d(TheWilsonTypeFermActFactory::Instance().createObject(fermact_4d,
								  fermacttop_4d,
								  fermact_path_4d));
  
  
    Handle< FermState<T,P,Q> > state_4d(S_f_4d->createState(u,
							    state_info_xml_4d,
							    state_info_path_4d));
  
  
    //-------------------------------------------------------------------------------
    Handle< LinearOperator<T> > A5(S_f_5d->linOp4D(state_5d,S_f_5d->getQuarkMass(),input.param.invParam));
    Handle< LinearOperator<T> > A4(S_f_4d->linOp(state_4d));
  
    LatticeFermion psi, chi;
  
    random(psi);
    random(chi);
    
    LatticeFermion tmp1,tmp2;
    QDPIO::cout << "A5 plus" << endl;
    (*A5)(tmp1, psi, PLUS);
    DComplex nn5_plus  = innerProduct(chi, tmp1);
    QDPIO::cout << "A5 minus" << endl;
    (*A5)(tmp2, chi, MINUS);
    DComplex nn5_minus = innerProduct(tmp2, psi);
  
    LatticeFermion tmp3,tmp4;
    QDPIO::cout << "A4 plus" << endl;
    (*A4)(tmp3, psi, PLUS);
    DComplex nn4_plus  = innerProduct(chi, tmp3);
    QDPIO::cout << "A4 minus" << endl;
    (*A4)(tmp4, chi, MINUS);
    DComplex nn4_minus = innerProduct(tmp4, psi);
  
    push(xml_out,"Innerprods");
    write(xml_out, "nn5_plus", nn5_plus);
    write(xml_out, "nn4_plus", nn4_plus);
    write(xml_out, "norm_diff_plus", Real(norm2(tmp1-tmp3)));
    write(xml_out, "nn5_minus", nn5_minus);
    write(xml_out, "nn4_minus", nn4_minus);
    write(xml_out, "norm_diff_minus", Real(norm2(tmp2-tmp4)));
    pop(xml_out);
  
  }
  catch (const std::string& e) 
  {
    QDPIO::cerr << "Error in t_dwf4d: " << e << endl;
    throw;
  }
  catch (const char* e) 
  {
    QDPIO::cerr << "Error in t_dwf4d: " << e << endl;
    throw;
  }
  
  pop(xml_out);

  // Time to bolt
  Chroma::finalize();

  exit(0);
}
