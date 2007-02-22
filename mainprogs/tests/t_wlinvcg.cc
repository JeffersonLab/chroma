// $Id: t_wlinvcg.cc,v 3.1 2007-02-22 21:11:50 bjoo Exp $

#include <iostream>
#include <cstdio>

#include "chroma.h"
#include "primitives.h" // GTF: for PLUS



using namespace Chroma;


#include "chromabase.h"
#include "actions/ferm/invert/invcg2.h"

//! Conjugate-Gradient (CGNE) algorithm for a generic Linear Operator
/*! \ingroup invert
 * This subroutine uses the Conjugate Gradient (CG) algorithm to find
 * the solution of the set of linear equations
 *
 *   	    Chi  =  A . Psi
 *
 * where       A = M^dag . M
 *
 * Algorithm:

 *  Psi[0]  :=  initial guess;    	       Linear interpolation (argument)
 *  r[0]    :=  Chi - M^dag . M . Psi[0] ;     Initial residual
 *  p[1]    :=  r[0] ;	       	       	       Initial direction
 *  IF |r[0]| <= RsdCG |Chi| THEN RETURN;      Converged?
 *  FOR k FROM 1 TO MaxCG DO    	       CG iterations
 *      a[k] := |r[k-1]|**2 / <Mp[k],Mp[k]> ;
 *      Psi[k] += a[k] p[k] ;   	       New solution vector
 *      r[k] -= a[k] M^dag . M . p[k] ;        New residual
 *      IF |r[k]| <= RsdCG |Chi| THEN RETURN;  Converged?
 *      b[k+1] := |r[k]|**2 / |r[k-1]|**2 ;
 *      p[k+1] := r[k] + b[k+1] p[k];          New direction
 *
 * Arguments:
 *
 *  \param M       Linear Operator    	       (Read)
 *  \param chi     Source	               (Read)
 *  \param psi     Solution    	    	       (Modify)
 *  \param RsdCG   CG residual accuracy        (Read)
 *  \param MaxCG   Maximum CG iterations       (Read)
 *  \param n_count Number of CG iteration      (Write)
 *
 * Local Variables:
 *
 *  p   	       Direction vector
 *  r   	       Residual vector
 *  cp  	       | r[k] |**2
 *  c   	       | r[k-1] |**2
 *  k   	       CG iteration counter
 *  a   	       a[k]
 *  b   	       b[k+1]
 *  d   	       < p[k], A.p[k] >
 *  Mp  	       Temporary for  M.p
 *
 * Subroutines:
 *                             +               
 *  A       Apply matrix M or M  to vector
 *
 * Operations:
 *
 *  2 A + 2 Nc Ns + N_Count ( 2 A + 10 Nc Ns )
 */

void WlInvCG2(const LinearOperator& M,
	      const LatticeFermion& chi,
	      LatticeFermion& psi,
	      const Real& RsdCG, 
	      int MaxCG, 
	      int& n_count)
{
  const Subset& s = M.subset();

  LatticeFermion mp;
  Real a;
  Real b;
  Double c;
  Double d;
  
  Real rsd_sq = (RsdCG * RsdCG) * Real(norm2(chi,s));

  //                                            +
  //  r[0]  :=  Chi - A . Psi[0]    where  A = M  . M
    
  //                      +
  //  r  :=  [ Chi  -  M(u)  . M(u) . psi ]
  LatticeFermion r;
  r[s] = chi - M(M(psi, PLUS), MINUS);

  //  p[1]  :=  r[0]
  LatticeFermion p;
  p[s] = r;
  
  //  Cp = |r[0]|^2
  Double cp = norm2(r, s);   	       	   /* 2 Nc Ns  flops */

  QDPIO::cout << "WlInvCG: k = 0  cp = " << cp << "  rsd_sq = " << rsd_sq << endl;

  //  IF |r[0]| <= RsdCG |Chi| THEN RETURN;
  if ( toBool(cp  <=  rsd_sq) )
  {
    n_count = 0;
    return;
  }

  //
  //  FOR k FROM 1 TO MaxCG DO
  //
  for(int k = 1; k <= MaxCG; ++k)
  {
    //  c  =  | r[k-1] |**2
    c = cp;

    //  a[k] := | r[k-1] |**2 / < p[k], Ap[k] > ;
    //      	       	       	       	       	  +
    //  First compute  d  =  < p, A.p >  =  < p, M . M . p >  =  < M.p, M.p >
    //  Mp = M(u) * p
    mp[s] = M(p, PLUS);

    //  d = | mp | ** 2
    d = norm2(mp, s);	/* 2 Nc Ns  flops */

    a = Real(c)/Real(d);

    //  Psi[k] += a[k] p[k]
    psi[s] += a * p;	/* 2 Nc Ns  flops */

    //  r[k] -= a[k] A . p[k] ;
    //      	       +            +
    //  r  =  r  -  M(u)  . Mp  =  M  . M . p  =  A . p
    r[s] -= a * M(mp, MINUS);

    //  IF |r[k]| <= RsdCG |Chi| THEN RETURN;

    //  cp  =  | r[k] |**2
    cp = norm2(r, s);	                /* 2 Nc Ns  flops */

    QDPIO::cout << "WlInvCG: k = " << k << "  cp = " << cp << endl;

    if ( toBool(cp  <=  rsd_sq) )
    {
      n_count = k;
      return;
    }

    //  b[k+1] := |r[k]|**2 / |r[k-1]|**2
    b = Real(cp) / Real(c);

    //  p[k+1] := r[k] + b[k+1] p[k]
    p[s] = r + b*p;	/* Nc Ns  flops */
  }
  n_count = MaxCG;
  QDP_error_exit("too many CG iterations: count = %d", n_count);
}









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

  XMLFileWriter xml("t_lwldslash.xml");
  push(xml,"t_lwldslash");

  //! Test out dslash
  multi1d<LatticeColorMatrix> u(Nd);
  for(int m=0; m < u.size(); ++m)
    gaussian(u[m]);

  LatticeFermion psi, chi;
  random(psi);
  chi = zero;


  // Construct szin-like gauge field
  multi3d<ColorMatrix> u_tmp(Nd,2,Layout::sitesOnNode()/2);
  for(int m=0; m < u.size(); ++m)
  {
    multi1d<ColorMatrix> u_tt(Layout::sitesOnNode());
    QDP_extract(u_tt, u[m], all);

    // Use this guy
    for(int i=0; i < u_tt.size(); ++i)
    {
      int ii = i >> 1;
      int cb = i / (Layout::sitesOnNode()/2);
      u_tmp[m][cb][ii] = transpose(u_tt[ii]);
    }
  }



  pack_gauge_field(u_tmp);

  //! Create a linear operator
  WilsonDslash D(u);
  chi = D(psi, PLUS, 0);

  write(xml,"Nd", Nd);
  write(xml,"Nc", Nc);
  write(xml,"Ns", Ns);
  write(xml,"nrow", nrow);
  write(xml,"psi", psi);
  write(xml,"chi", chi);

  //! Create and try a more sophisticated operator
  Real Kappa = 0.1;
  PreconditionedWilson  M(u,Kappa);
  LatticeFermion eta;
  eta = M(psi, PLUS);

  write(xml,"eta", eta);

  pop(xml);

  // Time to bolt
  Chroma::finalize();

  exit(0);
}
