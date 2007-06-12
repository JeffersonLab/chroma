// $Id: t_neflinop.cc,v 3.1 2007-06-12 13:59:14 edwards Exp $

#include <iostream>
#include <cstdio>

#include "chroma.h"

#include "qdp_util.h"

using namespace Chroma;


//! Read a SZIN fermion. This is a simple memory dump reader.
/*!
 * \ingroup io
 *
 * \param q          lattice fermion ( Modify )
 * \param file       path ( Read )
 */    

void readSzinFerm(multi1d<LatticeFermion>& q, const string& file)
{
  BinaryFileReader cfg_in(file);

  //
  // Read propagator field
  //
  multi1d<int> lattsize_cb = Layout::lattSize();
  lattsize_cb[0] /= 2;  // checkerboard in the x-direction in szin

  // Read prop
  for(int s=0; s < q.size(); ++s)
  {
    for(int cb=0; cb < 2; ++cb)
    {
      for(int sitecb=0; sitecb < Layout::vol()/2; ++sitecb)
      {
	multi1d<int> coord = crtesn(sitecb, lattsize_cb);

	// construct the checkerboard offset
	int sum = 0;
	for(int m=1; m<Nd; m++)
	  sum += coord[m];

	// The true lattice x-coord
	coord[0] = 2*coord[0] + ((sum + cb) & 1);

	read(cfg_in, q[s], coord); 	// Read in a site propagator
      }
    }
  }

  cfg_in.close();
}



int main(int argc, char **argv)
{
  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);

  // Setup the layout
  const int foo[] = {4,4,4,4};
  multi1d<int> nrow(Nd);
  nrow = foo;  // Use only Nd elements
  Layout::setLattSize(nrow);
  Layout::create();

  XMLFileWriter xml("t_neflinop.xml");
  push(xml, "t_neflinop");

  //! Test out dslash
  multi1d<LatticeColorMatrix> u(Nd);
#if 1
  for(int m=0; m < u.size(); ++m)
    gaussian(u[m]);
#else
  XMLReader gauge_xml;
  readSzin(gauge_xml, u, "small.cfg");
#endif

  // Create the BC objects
  const int bnd[] = {1,1,1,1};
  multi1d<int> boundary(Nd);
  boundary = bnd;

  // Create a FermBC
  Handle<FermBC<multi1d<LatticeFermion> > >  fbc_a(new SimpleFermBC<multi1d<LatticeFermion> >(boundary));

  Real WilsonMass = 1.5;
  Real m_q = 0.1;
  int  N5  = 8;

  multi1d<LatticeFermion> psi(N5), chi(N5),saveUprec(N5), savePrec(N5),tt(N5);

  for(int m=0; m < N5; ++m) 
    random(psi[m]);

  for(int m=0; m < N5; ++m)
    random(chi[m]);

  
  //#define SHAMIR
#define BORICI
 {
  //Shamir DWF case
#ifdef SHAMIR
  UnprecDWFermActArray S_DWF(fbc_a, WilsonMass, m_q, N5);
  UnprecNEFFermActArray S_NEF(fbc_a, WilsonMass, 1.0, 0, m_q, N5);
#endif
  //Borici DWF case
#ifdef BORICI
  UnprecOvDWFermActArray S_DWF(fbc_a, WilsonMass, m_q, N5);
  UnprecNEFFermActArray S_NEF(fbc_a, WilsonMass, 1.0, 1.0, m_q, N5);
#endif

  Handle<const ConnectState> stateDWF(S_DWF.createState(u));
  Handle<const LinearOperator< multi1d<LatticeFermion> > > Adwf(S_DWF.linOp(stateDWF));

  Handle<const ConnectState> stateNEF(S_NEF.createState(u));
  Handle<const LinearOperator< multi1d<LatticeFermion> > > Anef(S_NEF.linOp(stateNEF));

  multi1d<LatticeFermion> tmp1(N5);
  (*Adwf)(tmp1, psi, PLUS);
  multi1d<LatticeFermion> tmp2(N5);
  (*Anef)(tmp2, psi, PLUS);

  saveUprec = tmp2 ;

  multi1d<LatticeFermion> diff(N5);
  for(int m=0; m < N5; ++m)
    diff[m] = tmp2[m]-tmp1[m] ;



  (*Anef)(tmp1, psi, PLUS);
  DComplex nn1 = innerProduct(chi[0], tmp1[0]);
  for(int m=1; m < N5; ++m)
    nn1 += innerProduct(chi[m], tmp1[m]);
  
  (*Anef)(tmp2, chi, MINUS);
  DComplex nn2 = innerProduct(tmp2[0], psi[0]);
  for(int m=1; m < N5; ++m)
    nn2 += innerProduct(tmp2[m], psi[m]);



  push(xml,"innerprods");
  write(xml, "norm_psi" , Real(norm2(psi )));
  write(xml, "norm_chi" , Real(norm2(chi )));
  write(xml, "norm_tmp1", Real(norm2(tmp1)));
  write(xml, "norm_tmp2", Real(norm2(tmp2)));
  write(xml, "norm_diff", Real(norm2(diff)));
  write(xml, "nn1", nn1);
  write(xml, "nn2", nn2);

  pop(xml);
}

  {
  //Shamir case
#ifdef SHAMIR
  EvenOddPrecDWFermActArray S_DWF(fbc_a, WilsonMass, m_q, N5);
  EvenOddPrecNEFFermActArray S_NEF(fbc_a, WilsonMass, 1.0, 0, m_q, N5);
#endif
  //Borici DWF case
#ifdef BORICI
  EvenOddPrecOvDWFermActArray S_DWF(fbc_a, WilsonMass, m_q, N5);
  EvenOddPrecNEFFermActArray S_NEF(fbc_a, WilsonMass, 1.0, 1.0, m_q, N5);
#endif

  Handle<const ConnectState> stateDWF(S_DWF.createState(u));
  Handle<const EvenOddPrecLinearOperator< multi1d<LatticeFermion> > > Adwf(S_DWF.linOp(stateDWF));

  Handle<const ConnectState> stateNEF(S_NEF.createState(u));
  Handle<const EvenOddPrecLinearOperator< multi1d<LatticeFermion> > > Anef(S_NEF.linOp(stateNEF));

  multi1d<LatticeFermion> tmp1(N5);
  (*Adwf)(tmp1, psi, PLUS);
  multi1d<LatticeFermion> tmp2(N5);
  (*Anef)(tmp2, psi, PLUS);

  (*Anef).unprecLinOp(savePrec,psi,PLUS);

  for(int m=0; m < N5; ++m){
    tmp2[m][rb[0]] = tmp1[m] ;
  }

  multi1d<LatticeFermion> diff(N5);
  for(int m=0; m < N5; ++m)
    diff[m] = tmp2[m]-tmp1[m] ;



  (*Anef)(tmp1, psi, PLUS);
  DComplex nn1 = innerProduct(chi[0], tmp1[0]);
  for(int m=1; m < N5; ++m)
    nn1 += innerProduct(chi[m], tmp1[m]);
  
  (*Anef)(tmp2, chi, MINUS);
  DComplex nn2 = innerProduct(tmp2[0], psi[0]);
  for(int m=1; m < N5; ++m)
    nn2 += innerProduct(tmp2[m], psi[m]);



  push(xml,"prec_innerprods");
  write(xml, "norm_psi" , Real(norm2(psi )));
  write(xml, "norm_chi" , Real(norm2(chi )));
  //write(xml, "norm_tmp1", Real(norm2(tmp1)));
  //write(xml, "norm_tmp2", Real(norm2(tmp2)));
  write(xml, "norm_diff", Real(norm2(diff)));
  write(xml, "nn1", nn1);
  write(xml, "nn2", nn2);

  pop(xml);
  }
 for(int m=0; m < N5; ++m)
   tt[m] = savePrec[m]- saveUprec[m] ;

 push(xml,"prec_minus_unprec_innerprods");
 write(xml, "norm_Prec_minus_Uprec", Real(norm2(tt)));
 pop(xml);

 //TEST THE individual pieces of the operator
 {
  //Shamir case
#ifdef SHAMIR
  EvenOddPrecDWFermActArray S_DWF(fbc_a, WilsonMass, m_q, N5);
  EvenOddPrecNEFFermActArray S_NEF(fbc_a, WilsonMass, 1.0, 0, m_q, N5);
#endif
  //Borici DWF case
#ifdef BORICI
  EvenOddPrecOvDWFermActArray S_DWF(fbc_a, WilsonMass, m_q, N5);
  EvenOddPrecNEFFermActArray S_NEF(fbc_a, WilsonMass, 1.0, 1.0, m_q, N5);
#endif

  Handle<const ConnectState> stateDWF(S_DWF.createState(u));
  Handle<const EvenOddPrecLinearOperator< multi1d<LatticeFermion> > > Adwf(S_DWF.linOp(stateDWF));

  Handle<const ConnectState> stateNEF(S_NEF.createState(u));
  Handle<const EvenOddPrecLinearOperator< multi1d<LatticeFermion> > > Anef(S_NEF.linOp(stateNEF));

  multi1d<LatticeFermion> tmp1(N5);
  multi1d<LatticeFermion> tmp2(N5),diff_ee(N5), diff_oo(N5),diff_eo(N5) ; 
  multi1d<LatticeFermion>  diff_oe(N5),diff_inv_ee(N5),diff_inv_oo(N5) ;

  (*Adwf).evenEvenLinOp(tmp1, psi, PLUS);
  (*Anef).evenEvenLinOp(tmp2, psi, PLUS);
  
  for(int m=0; m < N5; ++m)
    diff_ee[m] = tmp2[m]-tmp1[m] ;

  (*Adwf).oddOddLinOp(tmp1, psi, PLUS);
  (*Anef).oddOddLinOp(tmp2, psi, PLUS);
  
  for(int m=0; m < N5; ++m)
    diff_oo[m] = tmp2[m]-tmp1[m] ;

  (*Adwf).evenEvenInvLinOp(tmp1, psi, PLUS);
  (*Anef).evenEvenInvLinOp(tmp2, psi, PLUS);
  
  for(int m=0; m < N5; ++m)
    diff_inv_ee[m] = tmp2[m]-tmp1[m] ;

  (*Adwf).evenOddLinOp(tmp1, psi, PLUS);
  (*Anef).evenOddLinOp(tmp2, psi, PLUS);
  
  for(int m=0; m < N5; ++m)
    diff_eo[m] = tmp2[m]-tmp1[m] ;

  (*Adwf).oddEvenLinOp(tmp1, psi, PLUS);
  (*Anef).oddEvenLinOp(tmp2, psi, PLUS);
  
  for(int m=0; m < N5; ++m)
    diff_oe[m] = tmp2[m]-tmp1[m] ;

  push(xml,"prec_pieces_PLUS");
  write(xml, "norm_ee" , Real(norm2(diff_ee)));
  write(xml, "norm_oo" , Real(norm2(diff_oo)));
  write(xml, "norm_eo" , Real(norm2(diff_eo)));
  write(xml, "norm_oe" , Real(norm2(diff_oe)));

  write(xml, "norm_inv_ee" , Real(norm2(diff_inv_ee)));

  pop(xml);

  (*Adwf).evenEvenLinOp(tmp1, psi, MINUS);
  (*Anef).evenEvenLinOp(tmp2, psi, MINUS);
  
  for(int m=0; m < N5; ++m)
    diff_ee[m] = tmp2[m]-tmp1[m] ;

  (*Adwf).oddOddLinOp(tmp1, psi, MINUS);
  (*Anef).oddOddLinOp(tmp2, psi, MINUS);
  
  for(int m=0; m < N5; ++m)
    diff_oo[m] = tmp2[m]-tmp1[m] ;

  (*Adwf).evenEvenInvLinOp(tmp1, psi, MINUS);
  (*Anef).evenEvenInvLinOp(tmp2, psi, MINUS);
  
  for(int m=0; m < N5; ++m)
    diff_inv_ee[m] = tmp2[m]-tmp1[m] ;

  (*Adwf).evenOddLinOp(tmp1, psi, MINUS);
  (*Anef).evenOddLinOp(tmp2, psi, MINUS);
  
  for(int m=0; m < N5; ++m)
    diff_eo[m] = tmp2[m]-tmp1[m] ;

  (*Adwf).oddEvenLinOp(tmp1, psi, MINUS);
  (*Anef).oddEvenLinOp(tmp2, psi, MINUS);
  
  for(int m=0; m < N5; ++m)
    diff_oe[m] = tmp2[m]-tmp1[m] ;

  push(xml,"prec_pieces_MINUS");
  write(xml, "norm_ee" , Real(norm2(diff_ee)));
  write(xml, "norm_oo" , Real(norm2(diff_oo)));
  write(xml, "norm_eo" , Real(norm2(diff_eo)));
  write(xml, "norm_oe" , Real(norm2(diff_oe)));

  write(xml, "norm_inv_ee" , Real(norm2(diff_inv_ee)));

  pop(xml);

 }

 pop(xml);

 // Time to bolt
 Chroma::finalize();

 exit(0);
}
