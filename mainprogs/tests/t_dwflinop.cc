// $Id: t_dwflinop.cc,v 3.1 2007-06-12 13:59:14 edwards Exp $

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

  XMLFileWriter xml("t_dwflinop.xml");
  push(xml, "t_dwflinop");

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

  // Typedefs to save typing
  typedef LatticeFermion               T;
  typedef multi1d<LatticeColorMatrix>  P;
  typedef multi1d<LatticeColorMatrix>  Q;

  // Create a FermBC
  Handle< CreateFermState<T,P,Q> >  cfs(new CreateSimpleFermState<T,P,Q>(boundary));

  Real WilsonMass = 1.5;
  Real m_q = 0.1;
#if 1
  int  N5  = 8;
  UnprecDWFermActArray S_f(cfs, WilsonMass, m_q, N5);
#else
  int  N5  = 9;
  UnprecOvExtFermActArray S_f(cfs, WilsonMass, m_q, N5);
#endif

  Handle< FermState<T,P,Q> > state(S_f.createState(u));
  Handle< LinearOperatorArray<T> > A(S_f.linOp(state));

  multi1d<LatticeFermion> psi(N5), chi(N5);

#if 1
  for(int m=0; m < N5; ++m)
    random(psi[m]);

  for(int m=0; m < N5; ++m)
    random(chi[m]);
#else
  readSzinFerm(psi, "t_invert.psi0");
  readSzinFerm(chi, "t_invert.chi0");
#endif

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

 
#if 0
  if (N5 != Ls)
    QDP_error_exit("N5 != Ls");

  // Create a FermBC
  UnprecDWFermAct S_f_dwf(cfs, WilsonMass, m_q);
  LatticeDWFermion psi5, chi5, tmp1a;

  for(int m=0; m < N5; ++m)
  {
    pokeDW(psi5, psi[m], m);
    pokeDW(chi5, chi[m], m);
    pokeDW(tmp1a, tmp1[m], m);
  }

  const LinearOperatorProxy<LatticeDWFermion> B(S_f_dwf.linOp(u));

  LatticeDWFermion tmp5a;
  B(tmp5a, psi5, PLUS);
  DComplex nd1 = innerProduct(chi5, tmp5a);

  LatticeDWFermion tmp5b;
  B(tmp5b, chi5, MINUS);
  DComplex nd2 = innerProduct(tmp5b, psi5);

  push(xml,"innerprods_dwf");
  write(xml, "norm_psi5", Real(norm2(psi5)));
  write(xml, "norm_chi5", Real(norm2(chi5)));
  write(xml, "nd1", nd1);
  write(xml, "nd2", nd2);
  write(xml, "norm_tmp1_tmp1a", Real(norm2(tmp5a-tmp1a)));
  pop(xml);

#endif

#if 0

  // Test inverter
  QDPIO::cout << "|psi|^2 = " << norm2(psi) << endl;
  QDPIO::cout << "|chi|^2 = " << norm2(chi) << endl;

  chi = A(psi, PLUS);

  QDPIO::cout << "|chi|^2 = " << norm2(chi) << endl;

  Real RsdCG = 1.0e-5;
  int  MaxCG = 1000;
  int  n_count;
  InvCG2(*A, chi, psi, RsdCG, MaxCG, n_count);   // subtelty here: deref A to get object pointed by handle
 

  multi1d<LatticeFermion> psi2(N5), chi2(N5);

  readSzinFerm(psi2, "t_invert.psi1");
  readSzinFerm(chi2, "t_invert.chi1");

  for(int m=0; m < N5; ++m)
  {
    tmp1[m] = psi2[m] - psi[m];
    tmp2[m] = chi2[m] - chi[m];
  }

  QDPIO::cout << "|psi2-psi|^2 = " << norm2(tmp1) << endl;
  QDPIO::cout << "|chi2-chi|^2 = " << norm2(tmp2) << endl;


  LatticePropagator quark_propagator, quark_prop_source;
  XMLBufferWriter xml_buf;
  int ncg_had;

  quark_prop_source = 1;
  quark_propagator  = zero;

  quarkProp4(quark_propagator, xml_buf, quark_prop_source,
	     S_f, u,  CG_INVERTER, RsdCG, MaxCG, ncg_had);

  LatticePropagator q2;
  XMLReader prop_xml;
  readSzinQprop(prop_xml, q2, "szin.prop");

  QDPIO::cout << "|quark_prop|^2 = " << norm2(quark_propagator) << endl;
  QDPIO::cout << "|q2|^2 = "         << norm2(q2) << endl;
  QDPIO::cout << "|q2-quark_prop|^2 = " << norm2(q2-quark_propagator) << endl;
#endif

  pop(xml);

  // Time to bolt
  Chroma::finalize();

  exit(0);
}
