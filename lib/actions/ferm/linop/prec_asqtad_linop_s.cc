// $Id: prec_asqtad_linop_s.cc,v 1.1 2003-12-10 14:27:08 bjoo Exp $
/*! \file
 *  \brief Unpreconditioned Wilson linear operator
 */
// NEW $Id: asqtad_linop.cc 2003/11/13 steve

#include "chromabase.h"
#include "actions/ferm/linop/asqtad_linop_s.h"


//! Apply Asqtad staggered fermion linear operator
/*!
 * \ingroup linop
 *
 * The operator acts on the checkerboarded lattice
 *
 * \param psi 	  Pseudofermion field     	       (Read)
 * \param isign   Flag ( PLUS | MINUS )   	       (Read)

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
  AsqtadDslash D;

  for(cb = 0; cb < 2: cb++) { 
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
