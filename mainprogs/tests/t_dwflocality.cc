// $Id: t_dwflocality.cc,v 1.2 2004-03-17 15:33:03 kostas Exp $

#include <iostream>
#include <cstdio>

#include <stdlib.h>
#include <sys/time.h>

#include "chroma.h"


typedef multi1d<LatticeFermion>  MLF;

using namespace QDP;

int main(int argc, char **argv)
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  // Setup the layout
  //const int foo[] = {4,4,4,8};
  //const int boo[] = {1,1,1,0};

  
  multi1d<int> nrow(Nd); 
  multi1d<int> boundary(Nd);
  QDPIO::cout << "Enter lattice size" << endl;
  QDPIO::cin >> nrow;

  QDPIO::cout << "Enter BC" << endl;
  QDPIO::cin >> boundary;


  Layout::setLattSize(nrow);
  Layout::create();

  int Nt = nrow[3] ;
  XMLFileWriter xml("t_dwflocality.xml");
  push(xml,"DWF_locality");
  push(xml,"lattice");
  write(xml,"size",nrow);
  write(xml,"bc",boundary);
  pop(xml);

  int N5 ;
  Real WilsonMass ;
  Real m_q = 1.0;
  Real  RsdCG = 1.0e-7 ;
  int MaxCG = 5000 ;

  QDPIO::cout << "Enter Ls" << endl;
  QDPIO::cin >> N5;
  QDPIO::cout << "Enter WilsonMass" << endl;
  QDPIO::cin >> WilsonMass;
  push(xml,"ChiralParam");
  write(xml,"N5",N5);
  write(xml,"WilsonMass",WilsonMass);
  pop(xml);

  string cnf ;
  QDPIO::cout << "Enter SZIN gaugefield" << endl;
  QDPIO::cin >> cnf;

  push(xml,"Configuration");
  write(xml, "cfg", cnf);
  pop(xml);
  QDPIO::cout << "Read SZIN config from " << cnf << endl;


  multi1d<LatticeColorMatrix> u(Nd);    // Gauge field
  XMLReader gauge_xml;
  readSzin(gauge_xml, u, cnf);

  Double w_plaq, s_plaq, t_plaq, link;
  MesPlq(u, w_plaq, s_plaq, t_plaq, link);
  QDPIO::cout << " Initial plaqettes and link: " << w_plaq
              << " " << s_plaq << " " << t_plaq << " " << link << endl;
  push(xml,"gauge_observables");
  write(xml, "w_plaq", w_plaq);
  write(xml, "s_plaq", s_plaq);
  write(xml, "t_plaq", t_plaq);
  write(xml, "link", link);
  pop(xml);



  MLF psi(N5), chi(N5);

  // Create a FermBC with only periodic BC. Note the handle is on an abstract type.
  Handle<FermBC<MLF> >  fbc(new SimpleFermBC<MLF>(boundary));

  // DWDslash class can be optimised
  int n_count;

  InvType invType(CG_INVERTER) ;

  EvenOddPrecDWFermActArray S_pdwf(fbc, WilsonMass, m_q,N5);
  Handle<const ConnectState> state(S_pdwf.createState(u));

  LatticeFermion chi4;
  LatticeFermion res4;

  //! Create source
  QDPIO::cout << "Constructing source" << endl;
  multi1d<int> coor(4) ;
  coor[0]=coor[1]=coor[2]=0 ;
  coor[3]=Nt/2 ;
  chi4=zero ;
  srcfil(chi4, coor, 0,0);
  QDPIO::cout << "Done" << endl;
  
  QDPIO::cout << "source norm :" << norm2(chi4)<< endl;
  LatticeReal tt ;
  tt = localNorm2(chi4) ;
  
  chi = zero ;
  // Split the source to oposite walls according to chirality
  chi[0   ] = chiralProjectPlus(chi4) ;
  chi[N5-1] = chiralProjectMinus(chi4) ; 
  
QDPIO::cout << "5D source norm :" << norm2(chi)<< endl;

  S_pdwf.qpropT(psi, state, chi, invType, RsdCG, MaxCG, n_count);
  
  res4 = chiralProjectMinus(psi[0]) + chiralProjectPlus(psi[N5-1]) ;

  LatticeReal mag ;
  
  mag = sqrt(localNorm2(chi4 - 0.5*res4)) ;
  //mag = localNorm2(chi4) ;
  multi1d<int> c(4) ;
  {
    multi1d<Real> val(nrow[0]); 
    c = coor ;
    for (int x(0);x<nrow[0];x++){
      c[0] = x ;
      val[x] = peekSite(mag,c) ;
    }
    push(xml, "x_axis");
    write(xml,"val",val);
    pop(xml);
  }
  {
    multi1d<Real> val(nrow[1]); 
    c = coor ;
    for (int x(0);x<nrow[1];x++){
      c[1] = x ;
      val[x] = peekSite(mag,c) ;
      //QDPIO::cout << "source norm ("<<x<<"):" << peekSite(tt,c)<< endl;
    }
    push(xml, "y_axis");
    write(xml,"val",val);
    pop(xml);
  }
  {
    multi1d<Real> val(nrow[2]); 
    c = coor ;
    for (int x(0);x<nrow[2];x++){
      c[2] = x ;
      val[x] = peekSite(mag,c) ;
    }
    push(xml, "z_axis");
    write(xml,"val",val);
    pop(xml);
  }
  {
    multi1d<Real> val(nrow[3]); 
    c = coor ;
    for (int x(0);x<nrow[3];x++){
      c[3] = x ;
      val[x] = peekSite(mag,c) ;
    }
    push(xml, "t_axis");
    write(xml,"val",val);
    pop(xml);
  }

  pop(xml);

  // Time to bolt
  QDP_finalize();

  exit(0);
}
