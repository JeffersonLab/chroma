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

  Tower<T> psi(3);
  Tower<T> chi(3);
  Q u(Nd);
  P p(Nd);

  for(int mu=0; mu < Nd; mu++) {
    gaussian(u[mu]);
    gaussian(p[mu]);
    reunit(u[mu]);
  }

  Handle<FermState<T,P,Q> > fs(new PeriodicFermState<T,P,Q>(u));
  WilsonDslash dslash(fs);
  // --- Dslash Apply test 
  { 
 
    for(int cb=1; cb >=0; cb--) { 
      
      multi1d<T> psi_f(3);
      T chi_f; gaussian(chi_f);
      
      for(int i=0; i < 3; i++) {
	psi_f[i][rb[cb]] = zero;
      }
      Q pu(Nd);

      for(int mu=0; mu < Nd; mu++) { 
	pu[mu] = u[mu];
      }
    
      for(int i=0; i < 3; i++) { 
	Handle<FermState<T,P,Q> > fs_tmp(new PeriodicFermState<T,P,Q>(pu));
	WilsonDslash dslash_tmp(fs_tmp);
	dslash_tmp.apply(psi_f[i], chi_f, PLUS,cb);
	U tmp;
	for(int mu=0; mu < Nd;mu++) {
	  tmp = pu[mu];
	  pu[mu] = -p[mu]*tmp;
	}
      }
      // So now psi_f[i] should hold        [D(P^2U) chi, D(-PU)chi, D(U) chi]
      multi1d<T> result(3);
      result=zero;
      result[0][rb[cb]] = psi_f[0];
      result[1][rb[cb]] = psi_f[0] +   psi_f[1];
      result[2][rb[cb]] = psi_f[0] + Real(2)*psi_f[1] + psi_f[2];
  
      
      psi=zero; chi=zero;
      chi[0]= chi_f;
      chi[1]= chi_f;
      chi[2]= chi_f;

      dslash.applyTower(psi,chi,p, PLUS, cb);
      for(int i=0; i < 3; i++) { 
	T diff = zero;
	diff[rb[cb]] = psi[i]-result[i];
	QDPIO::cout << "level=" << i << " cb=" << cb << " diff =" << norm2(diff,rb[cb]) << endl; 
      }
    }
  }

  

 
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
  ds_t =zero;
  
  dslash.deriv(ds_t, Xt, Yt, PLUS);
  QDPIO::cout << "Dslash" << endl;
  for(int mu=0; mu < Nd; mu++) { 
    U diff = zero;
    diff = ds_u[mu]-ds_t[mu][0];
    QDPIO::cout << "mu = " << mu << " norm2(diff) = " << norm2(diff) << endl;
  }

  Real Mass=0.1;

  UnprecWilsonLinOp W(fs, Mass);
  ds_u = zero;
  ds_t = zero;
  W.deriv(ds_u, X, Y, PLUS);
  W.deriv(ds_t, Xt, Yt, PLUS);

  QDPIO::cout << "Unprec Wilson LinOp:" << endl;
  for(int mu=0; mu < Nd; mu++) { 
    U diff = zero;
    diff = ds_u[mu]-ds_t[mu][0];
    QDPIO::cout << "mu = " << mu << " norm2(diff) = " << norm2(diff) << endl;
  }

  ds_u = zero;
  ds_t = zero;
  W.deriv(ds_u, X, Y, MINUS);
  W.deriv(ds_t, Xt, Yt, MINUS);
  QDPIO::cout << "Unprec Wilson LinOp:" << endl;
  for(int mu=0; mu < Nd; mu++) { 
    U diff = zero;
    diff = ds_u[mu]-ds_t[mu][0];
    QDPIO::cout << "mu = " << mu << " norm2(diff) = " << norm2(diff) << endl;
  }

  // Time to bolt
  Chroma::finalize();

  exit(0);
}
