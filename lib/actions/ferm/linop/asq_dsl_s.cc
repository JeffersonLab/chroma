/*  $Id: asq_dsl_s.cc,v 1.1 2003-12-10 16:20:59 bjoo Exp $  */

#include "chromabase.h"
//#include "actions/ferm/fermacts/linop/asq_dsl_s.h
#include "asq_dsl_s.h"

using namespace QDP;

/* This routine is specific to staggered fermions! */

/* Asq_Dsl */

/* Description: */

/* This routine applies the "asq" or "asqtad" operator D' to Psi,  */
/* putting the result in Chi. */

/*	       Nd-1 */
/*	       --- */
/*	       \                     F */                     
/*   chi(x)  :=  >  isign eta  (x) [U  (x) psi(x+mu) */
/*	       /             mu	     mu */
/*	       --- */
/*	       mu=0 */

/*			+ c_3 U  (x) U  (x+mu) U  (x+2mu) psi(x+3mu) ] */
/*                             mu     mu        mu */

/*	             Nd-1 */
/*	             --- */
/*	             \                      +F */
/*                -    >  isign eta  (x)  [U  (x-mu) psi(x-mu) */
/*	             /             mu	    mu */
/*	             --- */
/*	             mu=0 */

/*                             +      +          + */
/*			+ c_3 U  (x) U  (x-2mu) U  (x-3mu) psi(x-3mu) ] */
/*                             mu     mu         mu */
/* Note the KS phase factors are already included in the U's! */

/* Arguments: */

/*  U_fat     'fat link' Gauge field				(Read) */
/*  U_triple    'triple-links' UUU				(Read) */
/*  Psi	      Pseudofermion field				(Read) */
/*  Chi	      Pseudofermion field				(Write) */
/*		      + */
/*  ISign      D' or D'  ( +1 | -1 ) respectively		(Read) */
/*  CB	      Checkerboard of OUTPUT vector			(Read) */

/* NOTE: the coefficient c_3 is included in u_triple! */

//
//  recoded by mcneile
//

//! Creation routine
void QDPStaggeredDslash::create(const multi1d<LatticeColorMatrix>& _u_fat, const multi1d<LatticeColorMatrix>& _u_triple)
{ 
  u_fat = _u_fat;
  u_triple = _u_triple;
  
//  XMLFileWriter xml_out("output2.xml");
//  push(xml_out, "more_tests");
//  Write(xml_out, u_fat);
//  Write(xml_out,u_triple);
//  pop(xml_out);

}

void QDPStaggeredDslash::apply (LatticeFermion& chi, const LatticeFermion& psi, enum PlusMinus isign, int cb) const
{
  START_CODE("asqdslash")

//  multi1d<LatticeColorMatrix> u_fat(Nd);
//  multi1d<LatticeColorMatrix> u_triple(Nd);

  XMLFileWriter xml_out("output2.xml");
  push(xml_out, "more_tests");
  Write(xml_out, u_fat);
  Write(xml_out,u_triple);
  pop(xml_out);
//  LatticeFermion chi;

//  need convention on isign
//

//void asq_dsl(
//	     multi1d<LatticeColorMatrix> & u_fat,
//	     multi1d<LatticeColorMatrix> & u_triple,
//	     LatticeFermion & psi,
//	     LatticeFermion & chi,
//	     int isign, int cb)


  LatticeFermion tmp_0;
  LatticeFermion tmp_1;
  LatticeFermion tmp_2;

  int mu;

// get the fat7 and triple links

//  Fat7_Links(u, u_fat, u0);         BAD IDEA!!!
//  Triple_Links(u, u_triple, u0);

  /* Forward one-hop and three-hop neigbhors */
  /* Note the KS phase factors are already included in the U's! */
  mu = 0;
  tmp_0 = shift(psi, FORWARD, mu);
  chi[rb[cb]] = u_fat[mu] * tmp_0;
  //  NEIGHBOUR(tmp_0, tmp_1, FORWARD, mu, 1);
  tmp_1 = shift(tmp_0, FORWARD, mu);
  tmp_2 = shift(tmp_1, FORWARD, mu);
  chi[rb[cb]] += u_triple[mu] * tmp_2;

  for(mu = 1;mu  <= ( Nd - 1); ++mu )
  {
    tmp_0 = shift(psi, FORWARD, mu);
    chi[rb[cb]] += u_fat[mu] * tmp_0;
    //    NEIGHBOUR(tmp_0, tmp_1, FORWARD, mu, 1);
    tmp_1 = shift(tmp_0, FORWARD, mu);
    tmp_2 = shift(tmp_1, FORWARD, mu);
    chi[rb[cb]] += u_triple[mu] * tmp_2;
  }

  /* Backward one-hop and three-hop neigbhors */
  /* Note the KS phase factors are already included in the U's! */
  for(mu = 0;mu  <= ( Nd - 1); ++mu )
  {
    chi[rb[cb]] -= shift(adj(u_fat[mu]), BACKWARD, mu) * shift(psi, 
BACKWARD, mu);
    tmp_0 = shift(adj(u_triple[mu]),  BACKWARD, mu) * shift(psi, BACKWARD, 
mu);
    // NEIGHBOUR(tmp_0, tmp_1, BACKWARD, mu, 1);
    tmp_1 = shift(tmp_0, BACKWARD, mu);
    tmp_2 = shift(tmp_1, BACKWARD, mu);

    chi[rb[cb]] -= tmp_2;
  }


  // this requires more thought
  const int MINUS = -1 ;
     
  if(isign == MINUS)
    chi = -chi;

}



