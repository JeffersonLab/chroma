// $Id: t_ovlap_bj.cc,v 1.2 2003-12-15 17:52:51 bjoo Exp $

#include <iostream>
#include <cstdio>

#include <stdlib.h>
#include <sys/time.h>

#include "chroma.h"
#include "actions/ferm/linop/lovlapms_w.h"
#include "actions/ferm/fermacts/zolotarev_state.h"
#include "actions/ferm/fermacts/zolotarev4d_fermact_bj_w.h"
#include "actions/ferm/linop/lovlapms_w.h"
using namespace QDP;

int main(int argc, char **argv)
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  // Setup the layout
  const int foo[] = {4,4,4,4};
  multi1d<int> nrow(Nd);
  nrow = foo;  // Use only Nd elements
  Layout::setLattSize(nrow);
  Layout::create();

  XMLFileWriter xml_out("t_ovlap.xml");

  //! Test out dslash
  multi1d<LatticeColorMatrix> u(Nd);
  for(int m=0; m < u.size(); ++m) {
    gaussian(u[m]);
    reunit(u[m]);
  }

  LatticeFermion psi, chi;

  random(psi);
  chi = zero;

  //! Create a linear operator
  QDPIO::cout << "Constructing Zolotarev4Dbj" << endl;


  //! Wilsoniums
  Real WilsonMass = -1.5;
  const UnprecWilsonFermAct D_w(WilsonMass);


  Real m_q = 0.1;
  XMLBufferWriter my_writer;

  //! N order Zolo approx, with wilson action.
  Zolotarev4DFermActBj   D(D_w, 
			   m_q,
			   22, 
			   1.0e-7,
			   5000,
			   my_writer);


  // Specify the approximation
  ZolotarevConnectState<LatticeFermion> connect_state(u, 0.05);

  // Make me a linop (this callls the initialise function)
  const LinearOperator<LatticeFermion>*  D_op = D.linOp((ConnectState &)connect_state);

  Double n2 = norm2(psi);
  psi /= n2;

  (*D_op)(chi,psi,(enum PlusMinus) PLUS);
  

  delete D_op;
  // Time to bolt
  QDP_finalize();

  exit(0);
}
