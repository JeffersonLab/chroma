// $Id: asqtad_linop_s.cc,v 1.1 2003-12-10 12:38:14 bjoo Exp $
/*! \file
 *  \brief Unpreconditioned Wilson linear operator
 */
// NEW $Id: asqtad_linop.cc 2003/11/13 steve

#include "chromabase.h"
//#include "actions/ferm/linop/asqtad_linop_s.h"
#include "asqtad_linop_s.h"

//! Creation routine
/*! \ingroup fermact
 *
 * \param _u_fat         fat7 links                        (Read)
 * \param _u_triple    triple links                        (Read)
 * \param _Mass        fermion mass   	                   (Read)
 */
void AsqtadLinOp::create(const multi1d<LatticeColorMatrix>& _u_fat, const multi1d<LatticeColorMatrix>& _u_triple, const Real& _Mass)
{
  Mass = _Mass;
  u_fat = _u_fat;
  u_triple = _u_triple;
  D.create(u_fat,u_triple);
  XMLFileWriter xml_out("output4.xml");
  push(xml_out, "more_tests");
  Write(xml_out, u_fat);
  Write(xml_out, u_triple);
  pop(xml_out);

//    CoeffWilsr_s = (AnisoP) ? Wilsr_s / xiF_0 : 1;
}


//! Apply Asqtad staggered fermion linear operator
/*!
 * \ingroup linop
 *
 * The operator acts on the checkerboarded lattice
 *
 * \param psi 	  Pseudofermion field     	       (Read)
 * \param isign   Flag ( PLUS | MINUS )   	       (Read)
 * \param cb      Checkerboard of OUTPUT VECTOR        (Read)
 */
void AsqtadLinOp::operator() (LatticeFermion& chi, const LatticeFermion& psi, enum PlusMinus isign) const
{
  LatticeFermion tmp;
  int cb;

  START_CODE("AsqtadLinOp");

  //
  //  Chi     =  2m Psi     +  D'      Psi
  //     E(O)          E(O)      EO(OE)   O(E)


  // THIS IS VERY NASTY, NEED TO THINK MORE ABOUT IT!!!
  
  if(isign == MINUS) cb=0;
    else cb=1;

  D.apply(tmp, psi, isign, cb);

  XMLFileWriter xml_out("output1.xml");
  
  push(xml_out, "Dslash");
  Write(xml_out, psi);
  Write(xml_out, tmp);
//  pop(xml_out);


  chi[rb[cb]] = 2*Mass*psi + tmp;
  
  Write(xml_out, chi);
  pop(xml_out);
  END_CODE("AsqtadLinOp");

//  return chi;
}
