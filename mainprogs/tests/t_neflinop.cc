// $Id: t_neflinop.cc,v 1.1 2004-08-23 21:07:47 kostas Exp $

#include <iostream>
#include <cstdio>

#include "chroma.h"

#include "qdp_util.h"

using namespace QDP;


//! Read a SZIN fermion. This is a simple memory dump reader.
/*!
 * \ingroup io
 *
 * \param q          lattice fermion ( Modify )
 * \param file       path ( Read )
 */    

void readSzinFerm(multi1d<LatticeFermion>& q, const string& file)
{
  BinaryReader cfg_in(file);

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
  QDP_initialize(&argc, &argv);

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
  UnprecDWFermActArray S_f(fbc_a, WilsonMass, m_q, N5);


  Handle<const ConnectState> state(S_f.createState(u));
  Handle<const LinearOperator< multi1d<LatticeFermion> > > A(S_f.linOp(state));

  multi1d<LatticeFermion> psi(N5), chi(N5);

  for(int m=0; m < N5; ++m) 
    random(psi[m]);

  for(int m=0; m < N5; ++m)
    random(chi[m]);

  multi1d<LatticeFermion> tmp1(N5);
  (*A)(tmp1, psi, PLUS);
  DComplex nn1 = innerProduct(chi[0], tmp1[0]);
  for(int m=1; m < N5; ++m)
    nn1 += innerProduct(chi[m], tmp1[m]);

  multi1d<LatticeFermion> tmp2(N5);
  (*A)(tmp2, chi, MINUS);
  DComplex nn2 = innerProduct(tmp2[0], psi[0]);
  for(int m=1; m < N5; ++m)
    nn2 += innerProduct(tmp2[m], psi[m]);

  push(xml,"innerprods");
  write(xml, "norm_psi", Real(norm2(psi)));
  write(xml, "norm_chi", Real(norm2(chi)));
  write(xml, "nn1", nn1);
  write(xml, "nn2", nn2);
  pop(xml);

 
  pop(xml);

  // Time to bolt
  QDP_finalize();

  exit(0);
}
