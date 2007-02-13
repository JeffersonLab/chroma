// $Id: t_temp_prec.cc,v 3.1 2007-02-13 22:25:04 bjoo Exp $
/*! \file
 *  \brief Test 4d fermion actions
 */

#include <iostream>
#include <cstdio>

#include "chroma.h"

#include "central_tprec_linop.h"
#include "actions/ferm/linop/unprec_s_cprec_t_wilson_linop_w.h"
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



int main(int argc, char **argv)
{
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

  // Now need to create a SimpleFermState< LatticeFermion,  
  //                                       multi1d<LatticeColorMatrix>, 
  //                                       multi1d<LatticeColorMatrix> >
  multi1d<int> boundary(Nd);
  boundary[0] = 1;
  boundary[1] = 1;
  boundary[2] = 1;
  boundary[3] = -1;

  Handle< FermState< LatticeFermion, 
                     multi1d<LatticeColorMatrix>, 
                     multi1d<LatticeColorMatrix> > > fs( new SimpleFermState< LatticeFermion, 
					                                      multi1d<LatticeColorMatrix>, 
					                                      multi1d<LatticeColorMatrix> >(boundary, u) );

  AnisoParam_t aniso;
  aniso.anisoP = true;
  aniso.t_dir = 3;
  aniso.xi_0 = 5;
  aniso.nu =0.8;

  Real Mass = 0.02;

  // Create the reference UnprecWilson LinOp
  UnprecWilsonLinOp D_w( fs, Mass, aniso );
  UnprecSCprecTWilsonLinOp D_temp_prec(fs, Mass, aniso );


  LatticeFermion chi;
  LatticeFermion psi1, psi2, tmp1, tmp2;
  LatticeFermion diff;

  psi1 = zero;
  psi2 = zero;
  
  gaussian(chi);
  
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
  
  
  // Time to bolt
  Chroma::finalize();

  exit(0);
}
