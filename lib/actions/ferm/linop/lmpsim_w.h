// -*- C++ -*-
// $Id: lmpsim_w.h,v 1.6 2003-11-09 22:35:19 edwards Exp $

// #pragma ident "Id"

/*! \file
 *  \brief Preconditioned Wilson linear operator
 */

#ifndef __wilson_prec_w_h__
#define __wilson_prec_w_h__

#include "linearop.h"
#include "dslash_w.h"

using namespace QDP;

//! Preconditioned Wilson-Dirac operator
/*!
 * \ingroup linop
 *
 * This routine is specific to Wilson fermions!
 *
 *                                                    ~      ~+
 * This subroutine applies the preconditioned matrix  M  or  M   the vector
 * Psi
 *
 *      	       	   {   ~
 *      	       	   {   M(U) . Psi      	       if  ISign = PLUS
 *      	   Chi  =  {
 *      	       	   {   ~   +
 *      	       	   {   M(U)  . Psi     	       if  ISign = MINUS
 *
 * Algorithm:
 *
 * The kernel for Wilson fermions is
 *
 *      M  =  1 - k D'
 *
 * Note, we also allow for a parity breaking (1 + i H gamma_5 ) term,
 * if desired, in which case the preconditioning changes somewhat!
 *
 * where we use the notation  D'  for "D slash". In a "even-odd" or "red-black"
 * basis  D'  may be written as
 *
 *      [       1       	   -k D'       ]
 *      [        E,E    	        E,O    ]
 *      [       	       	       	       ]
 *      [   -k D'       	       1       ]
 *      [        O,E    	        O,O    ]
 *
 * The preconditioning consists of defining the lower and upper triangular
 * parts of this matrix,  L  and  U, and using the matrix
 *
 *      	       	       	   [   	   1   	         	   0   	       	   ]
 *      ~      -1    -1 	   [   	    E,E	         	       	       	   ]
 *      M  =  L   M U   	=  [   	       	       	            2    	   ]
 *      	       	       	   [   	   0   	       	    1    - k  D'    D' 	   ]
 *      	       	       	   [   	       	             O,O        O,E  E,O   ]
 *
 * instead of  M [Ref. DeGrand & Rossi, UCSD-HET90-01]. This means that we only
 * need the pseudofermion fields on the odd sites: the number of Flops required
 * is therefore unchanged, but the convergence properties are improved.
 *
 * Arguments:
 *
 *  Psi 	Pseudofermion field     	       (Read)
 *  Kappa       Hopping parameter       	       (Read)
 *  Chi 	Pseudofermion field     	       (Write)
 *  ISign       Flag ( PLUS | MINUS )   	       (Read)
 *
 * Local variables:
 *
 *  Tmp 	       Temporary (lives on even sites)
 *
 * Operations:
 *
 *  2 ( Nc Ns + DSlash )  flops 
 */

class PreconditionedWilson : public LinearOperator<LatticeFermion>
{
public:
  //! Full constructor
  PreconditionedWilson(const multi1d<LatticeColorMatrix>& _u, const Real& _Kappa)
    {create(_u,_Kappa);}

  //! Destructor is automatic
  ~PreconditionedWilson() {}

  //! Only defined on the odd subset
  const OrderedSubset& subset() const {return rb[1];}

  //! Creation routine
  void create(const multi1d<LatticeColorMatrix>& _u, const Real& _Kappa);

  //! Apply the operator onto a source vector
  LatticeFermion operator() (const LatticeFermion& psi, enum PlusMinus isign) const;

private:
  Real Kappa;
  multi1d<LatticeColorMatrix> u;
  WilsonDslash D;
};

#endif
