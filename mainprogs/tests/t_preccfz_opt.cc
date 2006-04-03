// $Id: t_preccfz_opt.cc,v 3.0 2006-04-03 04:59:15 edwards Exp $

#include <iostream>
#include <cstdio>

#include "chroma.h"
#include "actions/ferm/fermacts/zolotarev.h"

using namespace Chroma;


int main(int argc, char **argv)
{
  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);

  // Setup the layout
  const int foo[] = {2,2,2,2};
  multi1d<int> nrow(Nd);
  nrow = foo;  // Use only Nd elements
  Layout::setLattSize(nrow);
  Layout::create();

  XMLFileWriter xml("t_preccfz_opt.xml");
  push(xml, "t_preccfz_opt");

  multi1d<LatticeColorMatrix> u(Nd);
  for(int m=0; m < u.size(); ++m)
    gaussian(u[m]);

  // Create a FermBC with only periodic BC. Note the handle is on an abstract type.
  Handle<FermBC< multi1d<LatticeFermion> > >  fbc_a(new PeriodicFermBC< multi1d<LatticeFermion> >);

  // DWDslash class can be optimised
  
  EvenOddPrecOvlapContFrac5DFermActParams p;
  p.Mass =Real(0.06);
  p.RatPolyDeg=6;
  p.approximation_type=COEFF_TYPE_ZOLOTAREV;
  p.OverMass = Real(1.4);
  p.ApproxMin = 0.66;
  p.ApproxMax = 6.4635;

  EvenOddPrecOvlapContFrac5DFermActArray S_pdwf(fbc_a,p);

  // We have two cases.
  Handle<const ConnectState> state(S_pdwf.createState(u));
  
  multi1d<Real> alpha;
  multi1d<Real> beta;
  Real scale_factor;
      
  const OverlapConnectState& ov_state = 
    dynamic_cast<const OverlapConnectState&>(*state);

  S_pdwf.init(scale_factor, alpha, beta, ov_state);
      
  // Compare linops
  int N5 = S_pdwf.size();
  {
    bool isLastZeroP = ( p.RatPolyDeg % 2 == 0 ) ? true : false;

    QDPEvenOddPrecOvlapContFrac5DLinOpArray D_qdp(state,
						  p.Mass,
						  p.OverMass,
						  N5,
						  scale_factor,
						  alpha,
						  beta,
						  isLastZeroP);
    
    OptEvenOddPrecOvlapContFrac5DLinOpArray D_opt(state,
						  p.Mass,
						  p.OverMass,
						  N5,
						  scale_factor,
						  alpha,
						  beta,
						  isLastZeroP);
    
    multi1d<LatticeFermion>  chi5a(N5), chi5b(N5), psi5(N5), tmp1(N5);
    for(int m=0; m < N5; ++m)
    {
      gaussian(psi5[m]);
    }
    chi5a = chi5b = zero;

    D_qdp(chi5a, psi5, PLUS);
    D_opt(chi5b, psi5, PLUS);

//    push(xml,"Plus");
//    write(xml,"chi5a",chi5a);
//    write(xml,"chi5b",chi5b);
//    pop(xml);

    tmp1 = chi5a;
    tmp1 -= chi5b;

    QDPIO::cout << "Test eo-prec and opt eo-prec CFZ linop PLUS" << endl
		<< "|qCFZ|^2 = " << norm2(chi5a) << endl
		<< "|oCFZ|^2 = " << norm2(chi5b) << endl
		<< "|qCFZ - oCFZ|^2 = " << norm2(tmp1) << endl;

    D_qdp(chi5a, psi5, MINUS);
    D_opt(chi5b, psi5, MINUS);

//    push(xml,"Minus");
//    write(xml,"chi5a",chi5a);
//    write(xml,"chi5b",chi5b);
//    pop(xml);

    tmp1 = chi5a;
    tmp1 -= chi5b;

    QDPIO::cout << "Test eo-prec and opt eo-prec CFZ linop MINUS" << endl
		<< "|qCFZ|^2 = " << norm2(chi5a) << endl
		<< "|oCFZ|^2 = " << norm2(chi5b) << endl
		<< "|qCFZ - oCFZ|^2 = " << norm2(tmp1) << endl;
  }

  QDPIO::cout << "\n\n\nDone" << endl;

  pop(xml);

  // Time to bolt
  Chroma::finalize();

  exit(0);
}
