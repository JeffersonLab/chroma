// $Id: t_precact.cc,v 1.15 2005-02-13 18:15:19 edwards Exp $

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
  foo &= WilsonTypeFermActsEnv::registered;
  return foo;
}

//! Check Qprop
void check_qprop(XMLWriter& xml_out, const string& prefix, 
		 const SystemSolver< LatticeFermion >& PP,
		 const SystemSolver< LatticeFermion >& UP)
{
  LatticeFermion source;
  LatticeFermion sol_unprec;
  LatticeFermion sol_prec;

  // Make some kind of source
  gaussian(source);
  Double s_norm = sqrt(norm2(source));
  Real factor = Real(1)/Real(s_norm);
  source *= factor;

  sol_prec=zero;
  sol_unprec=zero;

  // Do solution one;
  int n_count1 = PP(sol_prec, source);
  int n_count2 = UP(sol_unprec, source);

  LatticeFermion diff = sol_prec - sol_unprec;
  Double norm_diff = sqrt(norm2(diff));

  QDPIO::cout << "Prop Diff: " << norm_diff << endl;
  push(xml_out, prefix);
  write(xml_out, "n_count1", n_count1);
  write(xml_out, "n_count2", n_count2);
  write(xml_out, "prop_diff", norm_diff);
  pop(xml_out);
}


//! Check linops
void check_linops(XMLWriter& xml_out, const string& prefix,
  const EvenOddPrecLinearOperator< multi1d<LatticeFermion>, multi1d<LatticeColorMatrix> >& AP,
  const UnprecLinearOperator< multi1d<LatticeFermion>, multi1d<LatticeColorMatrix> >& AU)
{
  QDPIO::cout << "Check linops" << endl;

  const int N5 = AP.size();
  multi1d<LatticeFermion>  psi(N5), chi(N5);
  multi1d<LatticeFermion>  tmpx(N5), tmpy(N5);
  
  for(int m=0; m < N5; ++m)
  {
    gaussian(psi[m]);
    gaussian(chi[m]);
  }

  multi1d<LatticeFermion>  tmp1(N5), tmp2(N5);
  QDPIO::cout << "AP plus" << endl;
  AP.evenEvenLinOp(tmpx, psi, PLUS);
  AP.evenEvenInvLinOp(tmpy, tmpx, PLUS);
  Double trivial_diff_plus = zero;
  for(int m=0; m < N5; ++m)
    trivial_diff_plus += norm2(tmpy[m]-psi[m],even);

  AP.unprecLinOp(tmp1, psi, PLUS);
  DComplex nnP_plus = zero;
  for(int m=0; m < N5; ++m)
    nnP_plus += innerProduct(chi[m], tmp1[m]);

  QDPIO::cout << "AP minus" << endl;
  AP.evenEvenLinOp(tmpx, psi, MINUS);
  AP.evenEvenInvLinOp(tmpy, tmpx, MINUS);
  Double trivial_diff_minus = zero;
  for(int m=0; m < N5; ++m)
    trivial_diff_minus += norm2(tmpy[m]-psi[m],even);

  AP.unprecLinOp(tmp2, chi, MINUS);
  DComplex nnP_minus = zero;
  for(int m=0; m < N5; ++m)
    nnP_minus += innerProduct(tmp2[m], psi[m]);
  
  multi1d<LatticeFermion> tmp3(N5), tmp4(N5);
  QDPIO::cout << "AU plus" << endl;
  AU(tmp3, psi, PLUS);
  DComplex nnU_plus = zero;
  for(int m=0; m < N5; ++m)
    nnU_plus += innerProduct(chi[m], tmp3[m]);

  QDPIO::cout << "AU minus" << endl;
  AU(tmp4, chi, MINUS);
  DComplex nnU_minus = zero;
  for(int m=0; m < N5; ++m)
    nnU_minus += innerProduct(tmp4[m], psi[m]);
  
  push(xml_out,prefix+"LinOpInnerprods");
  write(xml_out, "trivial_diff_plus", trivial_diff_plus);
  write(xml_out, "nnP_plus", nnP_plus);
  write(xml_out, "nnU_plus", nnU_plus);
  Double norm_diff_plus = zero;
  for(int m=0; m < N5; ++m)
    norm_diff_plus += norm2(tmp1[m]-tmp3[m]);
  write(xml_out, "norm_diff_plus", norm_diff_plus);

  write(xml_out, "trivial_diff_minus", trivial_diff_minus);
  write(xml_out, "nnP_minus", nnP_minus);
  write(xml_out, "nnU_minus", nnU_minus);
  Double norm_diff_minus = zero;
  for(int m=0; m < N5; ++m)
    norm_diff_minus += norm2(tmp2[m]-tmp4[m]);
  write(xml_out, "norm_diff_minus", norm_diff_minus);
  pop(xml_out);
}


//! Apply the operator onto a source vector
/*! User should make sure deriv routines do a resize  */
multi1d<LatticeColorMatrix>
deriv(const EvenOddPrecLinearOperator< multi1d<LatticeFermion>, multi1d<LatticeColorMatrix> >& AP,
      const multi1d<LatticeFermion>& chi, const multi1d<LatticeFermion>& psi, 
      enum PlusMinus isign)
{
  // Need deriv of  (A_oo - D_oe*Ainv_ee*D_eo*psi_e)
  enum PlusMinus msign = (isign == PLUS) ? MINUS : PLUS;

  //
  // Make sure the deriv routines do a resize !!!
  //
  multi1d<LatticeColorMatrix> ds_u, ds_1;
  multi1d<LatticeFermion>  tmp1, tmp2, tmp3;

  //
  // NOTE: even with even-odd decomposition, the ds_u will still have contributions
  // on all cb. So, no adding of ds_1 onto ds_u under a subset
  //
  //  ds_u  =  chi^dag * A'_oo * psi
  AP.derivOddOddLinOp(ds_u, chi, psi, isign);

  //  ds_u  -=  chi^dag * D'_oe * Ainv_ee * D_eo * psi_o
  AP.evenOddLinOp(tmp1, psi, isign);
  AP.evenEvenInvLinOp(tmp2, tmp1, isign);
  AP.derivOddEvenLinOp(ds_1, chi, tmp2, isign);
  ds_u -= ds_1;

  //  ds_u  +=  chi^dag * D_oe * Ainv_ee * A'_ee * Ainv_ee * D_eo * psi_o
  AP.evenOddLinOp(tmp1, psi, isign);
  AP.evenEvenInvLinOp(tmp2, tmp1, isign);
  AP.evenOddLinOp(tmp1, chi, msign);
  AP.evenEvenInvLinOp(tmp3, tmp1, msign);
  AP.derivEvenEvenLinOp(ds_1, tmp3, tmp2, isign);
  ds_u += ds_1;

  //  ds_u  -=  chi^dag * D_oe * Ainv_ee * D'_eo * psi_o
  AP.evenOddLinOp(tmp1, chi, msign);
  AP.evenEvenInvLinOp(tmp3, tmp1, msign);
  AP.derivEvenOddLinOp(ds_1, tmp3, psi, isign);
  ds_u -= ds_1;

  return ds_u;
}



//! Check linops
void check_derivs(XMLWriter& xml_out, const string& prefix,
  const EvenOddPrecLinearOperator< multi1d<LatticeFermion>, multi1d<LatticeColorMatrix> >& AP,
  const UnprecLinearOperator< multi1d<LatticeFermion>, multi1d<LatticeColorMatrix> >& AU)
{
  QDPIO::cout << "Check derivs" << endl;

  const int N5 = AP.size();
  multi1d<LatticeFermion>  psi(N5), chi(N5);
  
  for(int m=0; m < N5; ++m)
  {
    gaussian(psi[m]);
    gaussian(chi[m]);
  }

  multi1d<LatticeColorMatrix>  ds_1, ds_2;
  QDPIO::cout << "AP plus" << endl;
  AP.derivUnprecLinOp(ds_1, chi, psi, PLUS);
  Double nnP_plus = norm2(ds_1);

  QDPIO::cout << "AP minus" << endl;
  AP.derivUnprecLinOp(ds_2, chi, psi, MINUS);
  Double nnP_minus = norm2(ds_2);
  Double norm_diff_prec = zero;
  for(int m=0; m < Nd; ++m)
    norm_diff_prec += norm2(ds_1[m]-ds_2[m]);

  multi1d<LatticeColorMatrix>  ds_tmp;
  AP.deriv(ds_tmp, chi, psi, PLUS);
  ds_tmp -= deriv(AP, chi, psi, PLUS);
  Double norm_diff_check_prec_plus = norm2(ds_tmp);

  AP.deriv(ds_tmp, chi, psi, MINUS);
  ds_tmp -= deriv(AP, chi, psi, MINUS);
  Double norm_diff_check_prec_minus = norm2(ds_tmp);

  multi1d<LatticeColorMatrix>  ds_3, ds_4;
  QDPIO::cout << "AU plus" << endl;
  AU.deriv(ds_3, chi, psi, PLUS);
  Double nnU_plus = norm2(ds_3);

  QDPIO::cout << "AP minus" << endl;
  AU.deriv(ds_4, chi, psi, MINUS);
  Double nnU_minus = norm2(ds_4);
  Double norm_diff_unprec = zero;
  for(int m=0; m < Nd; ++m)
    norm_diff_unprec += norm2(ds_3[m]-ds_4[m]);

  push(xml_out,prefix+"DerivInnerprods");
  write(xml_out, "norm_diff_check_prec_plus", norm_diff_check_prec_plus);
  write(xml_out, "norm_diff_check_prec_minus", norm_diff_check_prec_minus);

  write(xml_out, "nnP_plus", nnP_plus);
  write(xml_out, "nnU_plus", nnU_plus);
  Double norm_diff_plus = zero;
  for(int m=0; m < Nd; ++m)
    norm_diff_plus += norm2(ds_1[m]-ds_3[m]);
  write(xml_out, "norm_diff_plus", norm_diff_plus);
  write(xml_out, "norm_diff_prec", norm_diff_plus);

  write(xml_out, "nnP_minus", nnP_minus);
  write(xml_out, "nnU_minus", nnU_minus);
  Double norm_diff_minus = zero;
  for(int m=0; m < Nd; ++m)
    norm_diff_minus += norm2(ds_2[m]-ds_4[m]);
  write(xml_out, "norm_diff_minus", norm_diff_minus);
  write(xml_out, "norm_diff_unprec", norm_diff_plus);
  pop(xml_out);
}


/*
 * Input 
 */
// Parameters which must be determined from the XML input
// and written to the XML output
struct Param_t
{
  multi1d<int> nrow;		// Lattice dimension
  InvertParam_t   invParam;   // Inverter parameters
};

//! Mega-structure of all input
struct Test_input_t
{
  Param_t          param;
  Cfg_t            cfg;
  string           action_eo;
  string           action_un;
};

//! Parameters for running code
void read(XMLReader& xml, const string& path, Param_t& param)
{
  XMLReader paramtop(xml, path);
  read(paramtop, "nrow", param.nrow);
  read(paramtop, "InvertParam", param.invParam);
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

    // Read unprec action stuff
    {
      XMLReader xml_action(inputtop, "UnprecAction");
      std::ostringstream os;
      xml_action.print(os);
      input.action_un = os.str();
    }

    // Read prec action stuff
    {
      XMLReader xml_action(inputtop, "PrecAction");
      std::ostringstream os;
      xml_action.print(os);
      input.action_eo = os.str();
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
  QDP_initialize(&argc, &argv);

  QDPIO::cout << "linkage=" << linkage_hack() << endl;

  // Input parameter structure
  Test_input_t  input;

  // Instantiate xml reader for DATA
  XMLReader xml_in("DATA");

  // Read data
  read(xml_in, "/t_precact", input);

  // Specify lattice size, shape, etc.
  Layout::setLattSize(input.param.nrow);
  Layout::create();

  QDPIO::cout << "t_precact" << endl;

  // Read in the configuration along with relevant information.
  multi1d<LatticeColorMatrix> u(Nd);
  XMLReader gauge_file_xml, gauge_xml;

  // Start up the gauge field
  gaugeStartup(gauge_file_xml, gauge_xml, u, input.cfg);

  // Instantiate XML writer for XMLDAT
  XMLFileWriter xml_out("t_precact.xml");
  push(xml_out, "t_precact");

  proginfo(xml_out);    // Print out basic program info

  // Write out the input
  write(xml_out, "Input", xml_in);

  // Write out the config header
  write(xml_out, "Config_info", gauge_xml);

  // Check if the gauge field configuration is unitarized
  unitarityCheck(u);

  // Calculate some gauge invariant observables just for info.
  Double w_plaq, s_plaq, t_plaq, link;
  MesPlq(u, w_plaq, s_plaq, t_plaq, link);

  push(xml_out, "Observables");
  write(xml_out, "w_plaq", w_plaq);
  write(xml_out, "s_plaq", s_plaq);
  write(xml_out, "t_plaq", t_plaq);
  write(xml_out, "link", link);
  pop(xml_out);


  //
  // Initialize fermion action
  //
  std::istringstream  xml_s_eo(input.action_eo);
  std::istringstream  xml_s_un(input.action_un);
  XMLReader  fermacttop_eo(xml_s_eo);
  XMLReader  fermacttop_un(xml_s_un);
  const string fermact_path_eo = "/PrecAction/FermionAction";
  const string fermact_path_un = "/UnprecAction/FermionAction";
  string fermact_eo;
  string fermact_un;

  try
  {
    read(fermacttop_eo, fermact_path_eo + "/FermAct", fermact_eo);
    read(fermacttop_un, fermact_path_un + "/FermAct", fermact_un);
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

  QDPIO::cout << "PrecFermAct = " << fermact_eo << endl;
  QDPIO::cout << "UnprecFermAct = " << fermact_un << endl;


  // Deal with auxiliary (and polymorphic) state information
  // eigenvectors, eigenvalues etc. The XML for this should be
  // stored as a string called "stateInfo" in the param struct.

  try {

    // Make a reader for the stateInfo
    const string state_info_path_eo = "/PrecAction/StateInfo";
    const string state_info_path_un = "/UnprecAction/StateInfo";
    XMLReader state_info_xml_eo(fermacttop_eo,state_info_path_eo);
    XMLReader state_info_xml_un(fermacttop_un,state_info_path_un);


    // Preconditioned operator
    bool success = false; 

    QDPIO::cerr << "create prec = " << fermact_eo << endl;
    typedef EvenOddPrecWilsonTypeFermAct5D< LatticeFermion, multi1d<LatticeColorMatrix> >  EO5D;

    Handle<EO5D>
      S_f_eo(dynamic_cast<EO5D*>(TheWilsonTypeFermAct5DFactory::Instance().createObject(fermact_eo,
											fermacttop_eo,
											fermact_path_eo)));
  
    Handle<const ConnectState> state_eo(S_f_eo->createState(u,
							    state_info_xml_eo,
							    state_info_path_eo));
  
    // Unpreconditioned operator
    QDPIO::cerr << "create unprec = " << fermact_un << endl;
    typedef UnprecWilsonTypeFermAct5D< LatticeFermion, multi1d<LatticeColorMatrix> >  U5D;
    Handle<U5D>
      S_f_un(dynamic_cast<U5D*>(TheWilsonTypeFermAct5DFactory::Instance().createObject(fermact_un,
										       fermacttop_un,
										       fermact_path_un)));
  
    Handle<const ConnectState> state_un(S_f_un->createState(u,
							    state_info_xml_un,
							    state_info_path_un));
  

    // Sanity check
    if (S_f_eo->size() != S_f_un->size())
    {
      QDPIO::cout << "Prec and unprec fermacts do not have same size" << endl;
      QDP_abort(1);
    }

    //-------------------------------------------------------------------------------
    {
      Handle<const EvenOddPrecLinearOperator< multi1d<LatticeFermion>, multi1d<LatticeColorMatrix> > > AP(S_f_eo->linOp(state_eo));
      Handle<const UnprecLinearOperator< multi1d<LatticeFermion>, multi1d<LatticeColorMatrix> > > AU(S_f_un->linOp(state_un));
      
      QDPIO::cout << "Check bulk linops" << endl;
      check_linops(xml_out, "Bulk", *AP, *AU);
      xml_out.flush();

      Handle< const SystemSolver<LatticeFermion> > PP = S_f_eo->qprop(state_eo,
								      input.param.invParam);
      Handle< const SystemSolver<LatticeFermion> > UP = S_f_un->qprop(state_un,
								      input.param.invParam);

      QDPIO::cout << "Check qprop" << endl;
      check_qprop(xml_out, "Qprop", *PP, *UP);
      xml_out.flush();

      QDPIO::cout << "Check bulk derivatives" << endl;
      check_derivs(xml_out, "Bulk", *AP, *AU);
      xml_out.flush();
    }


    //-------------------------------------------------------------------------------
    {
      Handle<const EvenOddPrecLinearOperator< multi1d<LatticeFermion>, multi1d<LatticeColorMatrix> > > AP(S_f_eo->linOpPV(state_eo));
      Handle<const UnprecLinearOperator< multi1d<LatticeFermion>, multi1d<LatticeColorMatrix> > > AU(S_f_un->linOpPV(state_un));
      
      QDPIO::cout << "Check PV linops" << endl;
      check_linops(xml_out, "PV", *AP, *AU);
      xml_out.flush();

      QDPIO::cout << "Check PV derivatives" << endl;
      check_derivs(xml_out, "PV", *AP, *AU);
      xml_out.flush();
    }

  }
  catch (const std::string& e) 
  {
    QDPIO::cerr << "Error in t_precact: " << e << endl;
    throw;
  }
  catch (const char* e) 
  {
    QDPIO::cerr << "Error in t_precact: " << e << endl;
    throw;
  }
  
  pop(xml_out);

  // Time to bolt
  QDP_finalize();

  exit(0);
}
