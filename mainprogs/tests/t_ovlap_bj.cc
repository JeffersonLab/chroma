// $Id: t_ovlap_bj.cc,v 1.6 2003-12-17 14:54:34 bjoo Exp $

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
  const int foo[] = {4,4,4,8};
  multi1d<int> nrow(Nd);
  nrow = foo;  // Use only Nd elements
  Layout::setLattSize(nrow);
  Layout::create();

  XMLFileWriter xml_out("t_ovlap.xml");

  //! Test out dslash
  multi1d<LatticeColorMatrix> u(Nd);
  
  //! Cold (order) start
  // for(int j = 0; j < Nd; j++) { 
  //   u(j) = Real(1);
  // }

  // Hot (disordered) start
  for(int j=0; j < Nd; j++) { 
    random(u(j));
    reunit(u(j));
  }

  //! Wilsoniums
  Real WilsonMass = -1.5;
  const UnprecWilsonFermAct S_w(WilsonMass);

  Real m_q = 0.0;
  XMLBufferWriter my_writer;

  //! N order Zolo approx, with wilson action.
  Zolotarev4DFermActBj   D(S_w, 
			   m_q,
			   28, 
			   1.0e-7,
			   5000,
			   my_writer);


  // Specify the approximation
  ZolotarevConnectState<LatticeFermion> connect_state(u, 0.001, 2.5);

  // Make me a linop (this callls the initialise function)
  const LinearOperatorProxy<LatticeFermion> D_op(D.linOp((ConnectState &)connect_state));

  LatticeFermion psi;
  gaussian(psi);
  Double n2 = norm2(psi);
  psi /= n2;

  LatticeFermion s1, s2, s3, tmp2;
  s1 = s2 = s3 = tmp2 = zero;

  D_op(s1,psi,PLUS);
  D_op(s2,psi,MINUS);
  D_op(tmp2, psi, PLUS);
  D_op(s3, tmp2, MINUS);

  s3 *= 2;
  s3 -= s1;
  s3 -= s2;

  cout << "Circle Norm: " << sqrt(norm2(s3)) << endl;

  // Time to bolt
  QDP_finalize();

  exit(0);
}
