// $Id: t_temp_prec.cc,v 3.9 2009/11/14 20:01:46 eneil Exp $
/*! \file
 *  \brief Test 4d fermion actions
 */

#include <iostream>
#include <cstdio>

#include "chroma.h"


using namespace Chroma;

#include "actions/ferm/linop/eo3dprec_s_cprec_t_wilson_linop_w.h"
#include "actions/ferm/linop/eo3dprec_s_cprec_t_clover_linop_w.h"
//! To insure linking of code, place the registered code flags here
/*! This is the bit of code that dictates what fermacts are in use */
bool linkage_hack()
{
  bool foo = true;
  // All actions
  foo &= WilsonTypeFermActsEnv::registerAll();
  return foo;
}



int main(int argc, char **argv)
{
  // This test has trouble with Nc=2, so make sure we're using 3 colors
#if QDP_NC == 3
  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);

  QDPIO::cout << "linkage=" << linkage_hack() << endl;


  multi1d<int> nrow(Nd);
  for(int i=0; i < Nd-1; i++) { 
    nrow[i] = 4;
  }
  nrow[Nd-1] = 16;

  // Specify lattice size, shape, etc.
  Layout::setLattSize(nrow);
  Layout::create();

  QDPIO::cout << "t_temp_prec" << endl;
  
  // Start up a weak field
  struct Cfg_t config = { CFG_TYPE_WEAK_FIELD, "dummy" };
  multi1d<LatticeColorMatrix> u(Nd);
  XMLReader gauge_file_xml, gauge_xml;
  gaugeStartup(gauge_file_xml, gauge_xml, u, config);
  // Check if the gauge field configuration is unitarized
  unitarityCheck(u);
  
  // Instantiate XML writer for XMLDAT
  XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();
  push(xml_out, "t_temp_prec");
  proginfo(xml_out);    // Print out basic program info


  // Calculate some gauge invariant observables just for info.
  MesPlq(xml_out, "Observables", u);
  xml_out.flush();  
  pop(xml_out);

#if 0
  // Now need to create a SimpleFermState< LatticeFermion,  
  //                                       multi1d<LatticeColorMatrix>, 
  //                                       multi1d<LatticeColorMatrix> >
  multi1d<int> boundary(Nd);
  boundary[0] = 1;
  boundary[1] = 1;
  boundary[2] = 1;
  boundary[3] = -1;
  

#endif

  XMLReader xml_in("./t_temp_prec.ini.xml");
  
  std::string fermact, prec_fermact, ilu_fermact;
  read(xml_in, "/tempPrec/UnprecFermAct/FermAct", fermact);
  read(xml_in, "/tempPrec/PrecFermAct/FermAct", prec_fermact);
  read(xml_in, "/tempPrec/ILUPrecFermAct/FermAct", ilu_fermact);
  // Typedefs to save typing
  typedef LatticeFermion               T;
  typedef multi1d<LatticeColorMatrix>  P;
  typedef multi1d<LatticeColorMatrix>  Q;
  

  // Create unprec and Prec actions
  Handle< FermionAction<T,P,Q> >
    S_unprec(  TheFermionActionFactory::Instance().createObject(fermact,
								xml_in,
								"/tempPrec/UnprecFermAct")     );
  
  UnprecWilsonFermAct& S_un_wils = dynamic_cast< UnprecWilsonFermAct& >(*S_unprec);
  
  Handle< FermionAction<T,P,Q> >
    S_prec(    TheFermionActionFactory::Instance().createObject(prec_fermact,
								xml_in,
								"/tempPrec/PrecFermAct") );

  UnprecSpaceCentralPrecTimeWilsonFermAct& S_unprec_s_cprec_t = 
    dynamic_cast< UnprecSpaceCentralPrecTimeWilsonFermAct& >(*S_prec);


  Handle< FermionAction<T,P,Q> >
    S_iluprec_handle(    TheFermionActionFactory::Instance().createObject(ilu_fermact,
									   xml_in,
									   "/tempPrec/ILUPrecFermAct") );

  ILUPrecSpaceCentralPrecTimeWilsonTypeFermAct<T,P,Q>& S_iluprec=
    dynamic_cast<  ILUPrecSpaceCentralPrecTimeWilsonTypeFermAct<T,P,Q>& >(*S_iluprec_handle);

  
  
  // Now create a FermState
  Handle< FermState<T, P, Q> >  fs( S_un_wils.createState(u) );


  Handle< UnprecSpaceCentralPrecTimeLinearOperator<T,P,Q> > D_temp_prec_handle( S_unprec_s_cprec_t.linOp( fs ) );
  UnprecSCprecTWilsonLinOp& D_temp_prec = dynamic_cast<UnprecSCprecTWilsonLinOp&>(*D_temp_prec_handle);

  Handle< UnprecLinearOperator<T,P,Q> > D_w_handle( S_un_wils.linOp(fs) );
  UnprecWilsonLinOp& D_w = dynamic_cast< UnprecWilsonLinOp& >(*D_w_handle);

  Handle< ILUPrecSpaceCentralPrecTimeLinearOperator<T,P,Q> > D_tprec2_handle(  S_iluprec.linOp(fs) );
  ILUPrecSCprecTWilsonLinOp& D_temp_prec2 = dynamic_cast< ILUPrecSCprecTWilsonLinOp& >(*D_tprec2_handle);

  LatticeFermion chi;
  LatticeFermion psi1, psi2, tmp1, tmp2;
  LatticeFermion diff;

  psi1 = zero;
  psi2 = zero;
  
  gaussian(chi);
  D_temp_prec.getFermBC().modifyF(chi);
  
  // gamma_5 hermiticity:   psi1 = \gamma_5  C_R^{-1} \gamma_5 \chi
  tmp1 = Gamma(15)*chi;
  D_temp_prec.invCRightLinOp(tmp2, tmp1, PLUS);
  psi1 = Gamma(15)*tmp2;

  // psi2 = ( C_L^{-1} )^\dagger  \chi
  D_temp_prec.invCLeftLinOp(psi2, chi, MINUS);
  diff = psi2 - psi1;
  QDPIO::cout << " Gamma5_1  = " << sqrt( norm2(diff) / norm2(psi1) ) << endl;
  

  // Other way 
  // psi1 = \gamma_5  C_L^{-1} \gamma_5 \chi
  tmp1 = Gamma(15)*chi; 
  D_temp_prec.invCLeftLinOp(tmp2, tmp1, PLUS);
  psi1 = Gamma(15)*tmp2;

  // psi2 = ( C_R^{-1} )^\dagger  \chi
  D_temp_prec.invCRightLinOp(psi2, chi, MINUS);
  diff = psi2 - psi1;
  QDPIO::cout << " Gamma5_2  = " << sqrt( norm2(diff) / norm2(psi1) ) << endl;
  
  // Check  D_w = CL^{-1} CR^{-1} + D_s 
  // D_w
  D_w(psi1, chi, PLUS);

  // CL^{-1} CR^{-1} + D_s)
  D_temp_prec.invCRightLinOp(tmp1, chi, PLUS);
  D_temp_prec.invCLeftLinOp(tmp2, tmp1, PLUS);
  D_temp_prec.spaceLinOp(psi2, chi, PLUS);
  psi2+=tmp2;

  diff = psi2 - psi1;
  QDPIO::cout << " D_1 + = " << sqrt( norm2(diff)/norm2(psi1) ) << endl;

  
  // Check  D_w^\dag = CR^{-dagger} CL^{-dagger}  + D_s^\dag  (MINUS)
  // D_w 
  D_w(psi1, chi, MINUS);

  // CL^{-1} CR^{-1} + D_s  (MINUS)
  D_temp_prec.invCLeftLinOp(tmp1, chi, MINUS);  
  D_temp_prec.invCRightLinOp(tmp2, tmp1, MINUS);
  D_temp_prec.spaceLinOp(psi2, chi, MINUS);
  psi2+=tmp2;

  diff = psi2 - psi1;
  QDPIO::cout << " D_1 - = " << sqrt( norm2(diff)/norm2(psi1) ) << endl;

#if 0
  multi1d<LatticeColorMatrix> u_2(1);
  multi1d<LatticeColorMatrix> u_inv(1);
  
  // Hotstart
  HotSt(u_2);
  PScalar<PColorMatrix< RComplex<REAL>, 3> > munit;
  munit = Real(1).elem();
 

  // This used to test the invert 3by3 - I made that private now. 
  // Tested through context

  munit = Real(1).elem(); // Invert
  for(int site=all.start(); site <= all.end(); site++) { 
    u_2[0].elem(site) += munit;
    D_temp_prec.invert3by3( u_inv[0].elem(site), u_2[0].elem(site) );
  }

  LatticeColorMatrix identityP = u_2[0]*u_inv[0];
  LatticeColorMatrix ident;
  for(int site=all.start(); site <= all.end(); site++) { 
    identityP.elem(site) -= munit;
    ident.elem(site) = munit;
  }

  QDPIO::cout << "|| I  - u_2*inv ||/ || 1 || = " << sqrt( norm2( identityP ) / norm2(ident)  ) << endl;
#endif


#if 0
  // Used to test the actual T0 and InvT0 ops - these are private now.
  LatticeHalfFermion chi_half;
  LatticeHalfFermion psi_half;
  LatticeHalfFermion chi2_half;
  LatticeHalfFermion diff_half;

  gaussian(chi_half);
  D_temp_prec.TOp(psi_half, chi_half, PLUS);
  D_temp_prec.invTOp(chi2_half, psi_half, PLUS);

  diff_half = chi2_half - chi_half;
  QDPIO::cout << " || diff || / || chi_half || = " << sqrt(norm2(diff_half) / norm2(chi_half)) << endl;

  
  D_temp_prec.TOp(psi_half, chi_half, MINUS);
  D_temp_prec.invTOp(chi2_half, psi_half, MINUS);

  diff_half = chi2_half - chi_half;
  QDPIO::cout << " || diff || / || chi_half || = " << sqrt(norm2(diff_half) / norm2(chi_half)) << endl;
#endif

  // Test LeftInverse is inverse of Left (Both orderings, PLUS and MINUS)
  gaussian(chi);
  D_temp_prec.getFermBC().modifyF(chi);

  D_temp_prec.cLeftLinOp(psi1, chi, PLUS);
  D_temp_prec.invCLeftLinOp(tmp1, psi1, PLUS);

  diff = tmp1 - chi;
  QDPIO::cout << " CL CL_INV + = " << sqrt(norm2(diff)/norm2(chi)) << endl;


  D_temp_prec.cLeftLinOp(psi1, chi, MINUS);
  D_temp_prec.invCLeftLinOp(tmp1, psi1, MINUS);

  diff = tmp1 - chi;
  QDPIO::cout << " CL CL_INV - = " << sqrt(norm2(diff)/norm2(chi)) << endl;

  D_temp_prec.invCLeftLinOp(psi1, chi, PLUS);
  D_temp_prec.cLeftLinOp(tmp1, psi1, PLUS);

  diff = tmp1 - chi;
  QDPIO::cout << " CL_INV CL + = " << sqrt(norm2(diff)/norm2(chi)) << endl;

  D_temp_prec.invCLeftLinOp(psi1, chi, MINUS);
  D_temp_prec.cLeftLinOp(tmp1, psi1, MINUS);

  diff = tmp1 - chi;
  QDPIO::cout << " CL_INV CL - = " << sqrt(norm2(diff)/norm2(chi)) << endl;


  // Test C_R is inverse of C_R_Inverse both orders, PLUS & Minus
  gaussian(chi);
  D_temp_prec.getFermBC().modifyF(chi);
  D_temp_prec.cRightLinOp(psi1, chi, PLUS);
  D_temp_prec.invCRightLinOp(tmp1, psi1, PLUS);

  diff = tmp1 - chi;
  QDPIO::cout << " CR CR_INV + = " << sqrt(norm2(diff)/norm2(chi)) << endl;

  D_temp_prec.cRightLinOp(psi1, chi, MINUS);
  D_temp_prec.invCRightLinOp(tmp1, psi1, MINUS);

  diff = tmp1 - chi;
  QDPIO::cout << " CR CR_INV - = " << sqrt(norm2(diff)/norm2(chi)) << endl;

  D_temp_prec.invCRightLinOp(psi1, chi, PLUS);
  D_temp_prec.cRightLinOp(tmp1, psi1, PLUS);

  diff = tmp1 - chi;
  QDPIO::cout << " CR_INV CR + = " << sqrt(norm2(diff)/norm2(chi)) << endl;

  D_temp_prec.invCRightLinOp(psi1, chi, MINUS);
  D_temp_prec.cRightLinOp(tmp1, psi1, MINUS);

  diff = tmp1 - chi;
  QDPIO::cout << " CR_INV CR - = " << sqrt(norm2(diff)/norm2(chi)) << endl;


  // Now test against the normal D_op =  C_L^{-1}( 1 + C^L D_s C^R ) C_R^{-1}
  D_w(psi1, chi, PLUS);
  D_temp_prec.unprecLinOp(psi2, chi, PLUS);
  diff = psi2 - psi1;
  QDPIO::cout << " D_2 + = " << sqrt( norm2(diff) / norm2(psi1) ) << endl;
  
  D_w(psi1, chi, MINUS);
  D_temp_prec.unprecLinOp(psi2, chi, MINUS);
  diff = psi2 - psi1;
  QDPIO::cout << " D_2 - = " << sqrt( norm2(diff) / norm2(psi1) ) << endl;


  // ==============================================================================================
  // ILU Prec Op
  // ===========================================================================================


  QDPIO::cout << "ILU Prec Op: " << endl;
  //  ILUPrecSCprecTWilsonLinOp D_temp_prec2(fs, Mass, aniso );

  tmp1 = Gamma(15)*chi;
  D_temp_prec2.invCRightLinOp(tmp2, tmp1, PLUS);
  psi1 = Gamma(15)*tmp2;

  // psi2 = ( C_L^{-1} )^\dagger  \chi
  D_temp_prec2.invCLeftLinOp(psi2, chi, MINUS);
  diff = psi2 - psi1;
  QDPIO::cout << " Gamma5_1  = " << sqrt( norm2(diff) / norm2(psi1) ) << endl;
  

  // Other way 
  // psi1 = \gamma_5  C_L^{-1} \gamma_5 \chi
  tmp1 = Gamma(15)*chi; 
  D_temp_prec2.invCLeftLinOp(tmp2, tmp1, PLUS);
  psi1 = Gamma(15)*tmp2;

  // psi2 = ( C_R^{-1} )^\dagger  \chi
  D_temp_prec2.invCRightLinOp(psi2, chi, MINUS);
  diff = psi2 - psi1;
  QDPIO::cout << " Gamma5_2  = " << sqrt( norm2(diff) / norm2(psi1) ) << endl;

  // Test LeftInverse is inverse of Left (Both orderings, PLUS and MINUS)
  gaussian(chi);
  D_temp_prec2.getFermBC().modifyF(chi);
  D_temp_prec2.cLeftLinOp(psi1, chi, PLUS);
  D_temp_prec2.invCLeftLinOp(tmp1, psi1, PLUS);

  diff = tmp1 - chi;
  QDPIO::cout << " CL CL_INV + = " << sqrt(norm2(diff)/norm2(chi)) << endl;


  D_temp_prec2.cLeftLinOp(psi1, chi, MINUS);
  D_temp_prec2.invCLeftLinOp(tmp1, psi1, MINUS);

  diff = tmp1 - chi;
  QDPIO::cout << " CL CL_INV - = " << sqrt(norm2(diff)/norm2(chi)) << endl;

  D_temp_prec2.invCLeftLinOp(psi1, chi, PLUS);
  D_temp_prec2.cLeftLinOp(tmp1, psi1, PLUS);

  diff = tmp1 - chi;
  QDPIO::cout << " CL_INV CL + = " << sqrt(norm2(diff)/norm2(chi)) << endl;

  D_temp_prec2.invCLeftLinOp(psi1, chi, MINUS);
  D_temp_prec2.cLeftLinOp(tmp1, psi1, MINUS);

  diff = tmp1 - chi;
  QDPIO::cout << " CL_INV CL - = " << sqrt(norm2(diff)/norm2(chi)) << endl;


  // Test C_R is inverse of C_R_Inverse both orders, PLUS & Minus
  gaussian(chi);
  D_temp_prec2.getFermBC().modifyF(chi);
  D_temp_prec2.cRightLinOp(psi1, chi, PLUS);
  D_temp_prec2.invCRightLinOp(tmp1, psi1, PLUS);

  diff = tmp1 - chi;
  QDPIO::cout << " CR CR_INV + = " << sqrt(norm2(diff)/norm2(chi)) << endl;

  D_temp_prec2.cRightLinOp(psi1, chi, MINUS);
  D_temp_prec2.invCRightLinOp(tmp1, psi1, MINUS);

  diff = tmp1 - chi;
  QDPIO::cout << " CR CR_INV - = " << sqrt(norm2(diff)/norm2(chi)) << endl;

  D_temp_prec2.invCRightLinOp(psi1, chi, PLUS);
  D_temp_prec2.cRightLinOp(tmp1, psi1, PLUS);

  diff = tmp1 - chi;
  QDPIO::cout << " CR_INV CR + = " << sqrt(norm2(diff)/norm2(chi)) << endl;

  D_temp_prec2.invCRightLinOp(psi1, chi, MINUS);
  D_temp_prec2.cRightLinOp(tmp1, psi1, MINUS);

  diff = tmp1 - chi;
  QDPIO::cout << " CR_INV CR - = " << sqrt(norm2(diff)/norm2(chi)) << endl;


  // Now test against the normal  D_op =  C_L^{-1} unprecOp() C_R^{-1}
  D_w(psi1, chi, PLUS);
  D_temp_prec2.unprecLinOp(psi2, chi, PLUS);
  diff = psi2 - psi1;
  QDPIO::cout << " D_2 + = " << sqrt( norm2(diff) / norm2(psi1) ) << endl;
  
  D_w(psi1, chi, MINUS);
  D_temp_prec2.unprecLinOp(psi2, chi, MINUS);
  diff = psi2 - psi1;
  QDPIO::cout << " D_2 - = " << sqrt( norm2(diff) / norm2(psi1) ) << endl;

  // ==============================================================================================
  // ILU Prec Clover Op
  // ===========================================================================================

  CloverFermActParams p;
  p.Mass=0.01;
  p.clovCoeffR = 0.91;
  p.clovCoeffT = 1.07;
  p.u0 = 1.0;

  p.anisoParam.anisoP = true;
  p.anisoParam.t_dir = 3;
  p.anisoParam.xi_0 = 2.464;
  p.anisoParam.nu = 0.95;
  

  QDPIO::cout << "ILU Prec Clover Op: " << endl;
  ILUPrecSCprecTCloverLinOp D_temp_prec_clover(fs, p );

  UnprecCloverLinOp D_clov(fs, p);
  
  tmp1 = Gamma(15)*chi;
  D_temp_prec_clover.invCRightLinOp(tmp2, tmp1, PLUS);
  psi1 = Gamma(15)*tmp2;

  // psi2 = ( C_L^{-1} )^\dagger  \chi
  D_temp_prec_clover.invCLeftLinOp(psi2, chi, MINUS);
  diff = psi2 - psi1;
  QDPIO::cout << " Gamma5_1  = " << sqrt( norm2(diff) / norm2(psi1) ) << endl;
  

  // Other way 
  // psi1 = \gamma_5  C_L^{-1} \gamma_5 \chi
  tmp1 = Gamma(15)*chi; 
  D_temp_prec_clover.invCLeftLinOp(tmp2, tmp1, PLUS);
  psi1 = Gamma(15)*tmp2;

  // psi2 = ( C_R^{-1} )^\dagger  \chi
  D_temp_prec_clover.invCRightLinOp(psi2, chi, MINUS);
  diff = psi2 - psi1;
  QDPIO::cout << " Gamma5_2  = " << sqrt( norm2(diff) / norm2(psi1) ) << endl;

  // Test LeftInverse is inverse of Left (Both orderings, PLUS and MINUS)
  gaussian(chi);
  D_temp_prec_clover.getFermBC().modifyF(chi);
  D_temp_prec_clover.cLeftLinOp(psi1, chi, PLUS);
  D_temp_prec_clover.invCLeftLinOp(tmp1, psi1, PLUS);

  diff = tmp1 - chi;
  QDPIO::cout << " CL CL_INV + = " << sqrt(norm2(diff)/norm2(chi)) << endl;


  D_temp_prec_clover.cLeftLinOp(psi1, chi, MINUS);
  D_temp_prec_clover.invCLeftLinOp(tmp1, psi1, MINUS);

  diff = tmp1 - chi;
  QDPIO::cout << " CL CL_INV - = " << sqrt(norm2(diff)/norm2(chi)) << endl;

  D_temp_prec_clover.invCLeftLinOp(psi1, chi, PLUS);
  D_temp_prec_clover.cLeftLinOp(tmp1, psi1, PLUS);

  diff = tmp1 - chi;
  QDPIO::cout << " CL_INV CL + = " << sqrt(norm2(diff)/norm2(chi)) << endl;

  D_temp_prec_clover.invCLeftLinOp(psi1, chi, MINUS);
  D_temp_prec_clover.cLeftLinOp(tmp1, psi1, MINUS);

  diff = tmp1 - chi;
  QDPIO::cout << " CL_INV CL - = " << sqrt(norm2(diff)/norm2(chi)) << endl;


  // Test C_R is inverse of C_R_Inverse both orders, PLUS & Minus
  gaussian(chi);
  D_temp_prec_clover.getFermBC().modifyF(chi);
  D_temp_prec_clover.cRightLinOp(psi1, chi, PLUS);
  D_temp_prec_clover.invCRightLinOp(tmp1, psi1, PLUS);

  diff = tmp1 - chi;
  QDPIO::cout << " CR CR_INV + = " << sqrt(norm2(diff)/norm2(chi)) << endl;

  D_temp_prec_clover.cRightLinOp(psi1, chi, MINUS);
  D_temp_prec_clover.invCRightLinOp(tmp1, psi1, MINUS);

  diff = tmp1 - chi;
  QDPIO::cout << " CR CR_INV - = " << sqrt(norm2(diff)/norm2(chi)) << endl;

  D_temp_prec_clover.invCRightLinOp(psi1, chi, PLUS);
  D_temp_prec_clover.cRightLinOp(tmp1, psi1, PLUS);

  diff = tmp1 - chi;
  QDPIO::cout << " CR_INV CR + = " << sqrt(norm2(diff)/norm2(chi)) << endl;

  D_temp_prec_clover.invCRightLinOp(psi1, chi, MINUS);
  D_temp_prec_clover.cRightLinOp(tmp1, psi1, MINUS);

  diff = tmp1 - chi;
  QDPIO::cout << " CR_INV CR - = " << sqrt(norm2(diff)/norm2(chi)) << endl;


  // Now test against the normal  D_op =  C_L^{-1} unprecOp() C_R^{-1}
  D_clov(psi1, chi, PLUS);
  D_temp_prec_clover.unprecLinOp(psi2, chi, PLUS);
  diff = psi2 - psi1;
  QDPIO::cout << " D_2 + = " << sqrt( norm2(diff) / norm2(psi1) ) << endl;
  
  D_clov(psi1, chi, MINUS);
  D_temp_prec_clover.unprecLinOp(psi2, chi, MINUS);
  diff = psi2 - psi1;
  QDPIO::cout << " D_2 - = " << sqrt( norm2(diff) / norm2(psi1) ) << endl;
  
#if 0
  // Schur Style preconditioning
  EO3DPrecSCprecTWilsonLinOp D_schur_tprec(fs, Mass, aniso);

  QDPIO::cout << "Schur Style preconditioning tests " << endl;
  QDPIO::cout << "================================= " << endl;



  for(int cb=0; cb < 2; cb++) { 
    gaussian(chi,rb3[cb]);
    
    tmp1[rb3[cb]] = Gamma(15)*chi;
    D_schur_tprec.invCRightLinOp(tmp2, tmp1, PLUS,cb);
    psi1[rb3[cb]] = Gamma(15)*tmp2;
    
    // psi2 = ( C_L^{-1} )^\dagger  \chi
    D_schur_tprec.invCLeftLinOp(psi2, chi, MINUS,cb);
    diff[rb3[cb]] = psi2 - psi1;
    QDPIO::cout << "cb="<<cb<<" Gamma5_1  = " << sqrt( norm2(diff,rb3[cb]) / norm2(psi1,rb3[cb]) ) << endl;


    // Other way 
    // psi1 = \gamma_5  C_L^{-1} \gamma_5 \chi
    tmp1[rb3[cb]] = Gamma(15)*chi; 
    D_schur_tprec.invCLeftLinOp(tmp2, tmp1, PLUS, cb);
    psi1[rb3[cb]] = Gamma(15)*tmp2;

    // psi2 = ( C_R^{-1} )^\dagger  \chi
    D_schur_tprec.invCRightLinOp(psi2, chi, MINUS, cb);
    diff[rb3[cb]] = psi2 - psi1;
    QDPIO::cout << "cb="<<cb<< " Gamma5_2  = " << sqrt( norm2(diff, rb3[cb]) / norm2(psi1,rb3[cb]) ) << endl;

    // Test LeftInverse is inverse of Left (Both orderings, PLUS and MINUS)
    gaussian(chi,rb3[cb]);
    D_schur_tprec.cLeftLinOp(psi1, chi, PLUS,cb);
    D_schur_tprec.invCLeftLinOp(tmp1, psi1, PLUS,cb);

    diff[rb3[cb]] = tmp1 - chi;
    QDPIO::cout << "cb="<<cb<<" CL CL_INV + = " << sqrt(norm2(diff,rb3[cb])/norm2(chi,rb3[cb])) << endl;


    D_schur_tprec.cLeftLinOp(psi1, chi, MINUS,cb);
    D_schur_tprec.invCLeftLinOp(tmp1, psi1, MINUS,cb);

    diff[rb3[cb]] = tmp1 - chi;
    QDPIO::cout << "cb="<<cb<<" CL CL_INV - = " << sqrt(norm2(diff,rb3[cb])/norm2(chi,rb3[cb])) << endl;

    D_schur_tprec.invCLeftLinOp(psi1, chi, PLUS,cb);
    D_schur_tprec.cLeftLinOp(tmp1, psi1, PLUS,cb);

    diff[rb3[cb]] = tmp1 - chi;
    QDPIO::cout << "cb="<<cb<<" CL_INV CL + = " << sqrt(norm2(diff,rb3[cb])/norm2(chi,rb3[cb])) << endl;

    D_schur_tprec.invCLeftLinOp(psi1, chi, MINUS,cb);
    D_schur_tprec.cLeftLinOp(tmp1, psi1, MINUS,cb);

    diff[rb3[cb]] = tmp1 - chi;
    QDPIO::cout << "cb="<<cb<<" CL_INV CL - = " << sqrt(norm2(diff,rb3[cb])/norm2(chi,rb3[cb])) << endl;

    // Test C_R is inverse of C_R_Inverse both orders, PLUS & Minus
    gaussian(chi,rb3[cb]);
    D_schur_tprec.cRightLinOp(psi1, chi, PLUS,cb);
    D_schur_tprec.invCRightLinOp(tmp1, psi1, PLUS,cb);
    
    diff[rb3[cb]] = tmp1 - chi;
    QDPIO::cout << "cb="<<cb<< " CR CR_INV + = " << sqrt(norm2(diff,rb3[cb])/norm2(chi,rb3[cb])) << endl;
    
    D_schur_tprec.cRightLinOp(psi1, chi, MINUS,cb);
    D_schur_tprec.invCRightLinOp(tmp1, psi1, MINUS,cb);
    
    diff[rb3[cb]] = tmp1 - chi;
    QDPIO::cout << "cb="<<cb<< " CR CR_INV - = " << sqrt(norm2(diff,rb3[cb])/norm2(chi,rb3[cb])) << endl;
    
    D_schur_tprec.invCRightLinOp(psi1, chi, PLUS,cb);
    D_schur_tprec.cRightLinOp(tmp1, psi1, PLUS,cb);

    diff[rb3[cb]] = tmp1 - chi;
    QDPIO::cout << "cb="<<cb<<" CR_INV CR + = " << sqrt(norm2(diff,rb3[cb])/norm2(chi,rb3[cb])) << endl;
    
    D_schur_tprec.invCRightLinOp(psi1, chi, MINUS,cb);
    D_schur_tprec.cRightLinOp(tmp1, psi1, MINUS,cb);
    
    diff[rb3[cb]] = tmp1 - chi;
    QDPIO::cout << "cb="<<cb<<" CR_INV CR - = " << sqrt(norm2(diff,rb3[cb])/norm2(chi,rb3[cb])) << endl;
  }

  // D_w is unpreconditioned
  //  D_w1, chi, PLUS);
  gaussian(chi);
  D_schur_tprec.cRightLinOp(tmp1, chi, PLUS, 0);
  D_schur_tprec.cRightLinOp(tmp1, chi, PLUS, 1);

  D_w(tmp2, tmp1, PLUS);

  D_schur_tprec.cLeftLinOp(psi1, tmp2, PLUS, 0);
  D_schur_tprec.cLeftLinOp(psi1, tmp2, PLUS, 1);

  // --------------

  D_schur_tprec.evenEvenLinOp(psi2, chi, PLUS);
  D_schur_tprec.evenOddLinOp(tmp1, chi, PLUS);
  psi2[rb3[0]] += tmp1;

  D_schur_tprec.oddOddLinOp(psi2, chi, PLUS);
  D_schur_tprec.oddEvenLinOp(tmp1, chi, PLUS);
  psi2[rb3[1]] += tmp1;

 
  for(int cb=0; cb < 2; cb++) { 
    diff[rb3[cb]] = psi2 - psi1;
    QDPIO::cout << "cb="<<cb<<" D+ = " << sqrt(norm2(diff,rb3[cb])/norm2(chi,rb3[cb])) << endl;
  }

  gaussian(chi);
  D_w(psi1, chi, PLUS);
  D_schur_tprec.unprecLinOp(psi2,chi,PLUS);
  for(int cb=0; cb < 2; cb++) { 
    diff[rb3[cb]] = psi2 - psi1;
    QDPIO::cout << "cb="<<cb<<" D+ = " << sqrt(norm2(diff,rb3[cb])/norm2(chi,rb3[cb])) << endl;
  }


  gaussian(chi);
  D_schur_tprec.cLeftLinOp(tmp1, chi, MINUS, 0);
  D_schur_tprec.cLeftLinOp(tmp1, chi, MINUS, 1);

  D_w(tmp2, tmp1, MINUS);

  D_schur_tprec.cRightLinOp(psi1, tmp2, MINUS, 0);
  D_schur_tprec.cRightLinOp(psi1, tmp2, MINUS, 1);

  // --------------


  D_schur_tprec.evenEvenLinOp(psi2, chi, MINUS);
  D_schur_tprec.evenOddLinOp(tmp1, chi, MINUS);
  psi2[rb3[0]] += tmp1;

  D_schur_tprec.oddOddLinOp(psi2, chi, MINUS);
  D_schur_tprec.oddEvenLinOp(tmp1, chi, MINUS);
  psi2[rb3[1]] += tmp1;

 
  for(int cb=0; cb < 2; cb++) { 
    diff[rb3[cb]] = psi2 - psi1;
    QDPIO::cout << "cb="<<cb<<" D- = " << sqrt(norm2(diff,rb3[cb])/norm2(chi,rb3[cb])) << endl;
  }

  gaussian(chi);
  D_w(psi1, chi, MINUS);
  D_schur_tprec.unprecLinOp(psi2,chi,MINUS);
  for(int cb=0; cb < 2; cb++) { 
    diff[rb3[cb]] = psi2 - psi1;
    QDPIO::cout << "cb="<<cb<<" D- = " << sqrt(norm2(diff,rb3[cb])/norm2(chi,rb3[cb])) << endl;
  }



  // -----------------------------------------------------------------
  EO3DPrecSCprecTCloverLinOp D_schur_clov_tprec(fs,p);
  gaussian(chi);
  D_clov(psi1, chi, PLUS);
  D_schur_clov_tprec.unprecLinOp(psi2,chi,PLUS);
  for(int cb=0; cb < 2; cb++) { 
    diff[rb3[cb]] = psi2 - psi1;
    QDPIO::cout << "cb="<<cb<<" D+ = " << sqrt(norm2(diff,rb3[cb])/norm2(chi,rb3[cb])) << endl;
  }

  gaussian(chi);
  D_clov(psi1, chi, MINUS);
  D_schur_clov_tprec.unprecLinOp(psi2,chi,MINUS);
  for(int cb=0; cb < 2; cb++) { 
    diff[rb3[cb]] = psi2 - psi1;
    QDPIO::cout << "cb="<<cb<<" D- = " << sqrt(norm2(diff,rb3[cb])/norm2(chi,rb3[cb])) << endl;
  }
#endif
  
  // Time to bolt
  Chroma::finalize();

#endif

  exit(0);
}
