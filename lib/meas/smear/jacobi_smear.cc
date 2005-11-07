// $Id: jacobi_smear.cc,v 2.1 2005-11-07 21:18:59 edwards Exp $

#error "NOT FULLY CONVERTED - NOT SURE THIS IS REALLY NEEDED"

/* Do a covariant JACOBI smearing of a color vector field 
 * N iterations of the smearing procedure are performed
 * 
 * Resulting in J_{N}, using the recurrence:
 *
 *  J_{N}(x) = s_0 + kappa_sm D_{jacobi} J_{n-1}
 *
 *  where J_{0} = s_0
 *
 */

/* Arguments: */

/*  u -- gauge field ( Read ) */
/*  chi -- color vector field ( Modify ) */
/*  kappa_sm  -- The Jacobi smearing parameter kappa */
/*  n_jacobi -- The number of Jacobi applications */
/*  j_decay  -- direction of decay ( Read ) */

include(types.mh)

SUBROUTINE(jacobi_smear, chi, u, kappa_sm, n_jacobi, j_decay)

multi1d<LatticeColorMatrix> u(Nd);
LatticeColorVector chi;
Real kappa_sm;
int n_jacobi;
int j_decay;
{ /* Local variables */
  include(COMMON_DECLARATIONS)

  LatticeColorVector chi_zero;
  LatticeColorVector tmp;
  
  int i;
  int cb;
  Double factor;
  Real rfactor;
  START_CODE();


    
  for(cb = 0; cb < Nsubl; cb++) { 
    chi_zero[cb] = chi[cb];
  }

  /* Do n_jacobi iterations */
  for(i=0; i < n_jacobi; i++) { 

    /* tmp = D_{jacobi} chi_{n} */
    jacobi_term (u, chi, tmp, j_decay);

    /* chi_{n+1} = chi_{0} + kappa_sm D_{jacobi} chi_{n}  */
    for(cb = 0; cb < Nsubl; cb++) { 

      /* chi_{n+1} = chi_{0} */
      chi[cb] = chi_zero[cb];

      /* chi_{n+1} += kappa_sm D chi_{n}
                    = kappa_sm tmp
      */
      chi[cb] += tmp[cb] * kappa_sm;
    }

  }

    
  factor=TO_DOUBLE(kappa_sm) - TO_DOUBLE(JacobiKappaCrit);
  factor*= - JacobiAlphaFactor * n_jacobi;
  factor = exp(factor);
  rfactor = FLOAT(factor);
  for(cb = 0; cb < Nsubl; cb++) {
    chi[cb] = chi[cb] * rfactor;
  }
  
  push(xml_out,"jacobi_smearing");
write(xml_out, "n_jacobi", n_jacobi);
pop(xml_out);
  END_CODE();
}
