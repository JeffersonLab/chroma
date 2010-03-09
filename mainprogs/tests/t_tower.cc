#include <chroma.h>
#include <tower.h>

using namespace QDP;
using namespace Chroma;

int main(int argc, char **argv)
{
  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);
  typedef LatticeFermion T;
  typedef multi1d<LatticeColorMatrix> P;
  typedef multi1d<LatticeColorMatrix> Q;
  typedef LatticeColorMatrix U;

  // Setup the layout
  const int foo[] = {2,2,2,2};
  multi1d<int> nrow(Nd);
  nrow = foo;  // Use only Nd elements
  Layout::setLattSize(nrow);
  Layout::create();

  Tower<T> psi(4);
  Tower<T> chi(4);
  Q u(Nd);
  P p(Nd);

  for(int mu=0; mu < Nd; mu++) {
    gaussian(u[mu]);
    gaussian(p[mu]);
    reunit(u[mu]);
  }

  Handle<FermState<T,P,Q> > fs(new PeriodicFermState<T,P,Q>(u));

  WilsonDslash dslash(fs);


  chi = zero;
  psi = zero;
  gaussian(psi[0]);

  dslash.applyTower(chi,psi,p,PLUS,1);

  T X, Y;
  Tower<T> Xt(2);
  Tower<T> Yt(2);
  gaussian(X);
  gaussian(Y);
  Xt[0] = X;
  Yt[0] = Y;

  Q ds_u(Nd);
  ds_u = zero;
  dslash.deriv(ds_u, X, Y, PLUS);

  TowerArray<U> ds_t(2);
  dslash.deriv(ds_t, Xt, Yt, PLUS);

  for(int mu=0; mu < Nd; mu++) { 
    U diff = zero;
    diff = ds_u[mu]-ds_t[mu][0];
    QDPIO::cout << "mu = " << mu << " norm2(diff) = " << norm2(diff) << endl;
  }


  // Time to bolt
  Chroma::finalize();

  exit(0);
}
