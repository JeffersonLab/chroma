/*  $Id: asq_dsl_s.h,v 1.1 2003-12-10 12:38:13 bjoo Exp $  */

/* This routine is specific to staggered fermions! */

#ifndef __asqdslash_h__
#define __asqdslash_h__

#include "linearop.h"

using namespace QDP;

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

class QDPStaggeredDslash : public DslashLinearOperator<LatticeFermion>
{  
public:
  //! Empty constructor. Must use create later
  QDPStaggeredDslash() {}
 
  //! Full constructor
  QDPStaggeredDslash(const multi1d<LatticeColorMatrix>& _u_fat, const multi1d<LatticeColorMatrix>& _u_triple) 
  {create(_u_fat,_u_triple);}
 
  //! Creation routine  
  void create(const multi1d<LatticeColorMatrix>& _u_fat, const multi1d<LatticeColorMatrix>& _u_triple);
 
  //! No real need for cleanup here
  ~QDPStaggeredDslash() {}

/* Arguments: */

/*  u_fat     Fat7 links                	  		(Read) */
/*  u_triple  triple links					(Read) */
/*  Psi	      Pseudofermion field - Source		        (Read) */
/*		      + */
/*  ISign      D' or D'  ( +1 | -1 ) respectively		(Read) */
/*  CB	      Checkerboard of OUTPUT vector			(Read) */

  void apply (LatticeFermion& chi, const LatticeFermion& psi, enum PlusMinus isign, int cb) const;
  
  //! Subset is all here
  const OrderedSubset& subset() const {return all;}
    
private:
  multi1d<LatticeColorMatrix> u_fat;
  multi1d<LatticeColorMatrix> u_triple;
};

#endif


