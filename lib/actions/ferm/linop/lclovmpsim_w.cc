// $Id: lclovmpsim_w.cc,v 1.1 2005-02-17 02:52:36 edwards Exp $

#error "NOT FULLY CONVERTED"

/* This routine is specific to Wilson fermions! */

/* LCLOVMPSIM */
/*      	       	       	       	       	       	     ~  */
/* This subroutine applies the preconditioned matrix  M  or  M   the vector */
/* Psi, */

/*                  {   ~ */
/*                  {   M(U) . Psi              if  ISign = PLUS */
/*          Chi  =  { */
/*                  {   ~   + */
/*                  {   M(U)  . Psi             if  ISign = MINUS */

/* Algorithm: */

/* The kernel for Wilson fermions with clover action is */

/*      M  =  A - k D' */

/* where we use the notation  D'  for "D slash". In a "even-odd" or "red-black" */
/* basis  D'  may be written as */

/*      [      A       	   -k D'       ] */
/*      [        E,E    	        E,O    ] */
/*      [       	       	       	       ] */
/*      [   -k D'       	      A        ] */
/*      [        O,E    	        O,O    ] */

/* The preconditioning consists of left multiplying by */


/*      [       1      	       0       ] */
/*      [        E,E    	        E,O    ] */
/*      [       	       	       	       ] */
/*      [  k D'    A^(-1)       1       ] */
/*      [      O,E   E,E	        O,O    ] */


/* Resulting in a new  M */

/*      [      A       	            -k D'                  ] */
/*      [        E,E    	                 E,O               ] */
/*      [       	       	       	                           ] */
/*      [      0       	   A     - k^2 D'    A^(-1)  D'    ] */
/*      [        O,E    	     O,O         O,E   E,E     E,O */


/* This routine computes */

/*      ~  */
/*      M  =  A(o,o) - k^2 * D'(o,e) . A^-1(e,e) . D'(e,o) */


/* Arguments: */

/*  U   	          Gauge field        	               (Read) */
/*  clov           'clover term' off-diag part          (Read) */
/*  diag_clov      'clover term' diag part              (Read) */
/*  invclov        inverse 'clover term' off-diag part  (Read) */
/*  diag_invclov   inverse 'clover term' diag part      (Read) */
/*  Psi 	          Pseudofermion field     	       (Read) */
/*  Kappa_ds       Dslash hopping parameter             (Read) */
/*  Chi 	          Pseudofermion field     	       (Write) */
/*  ISign          Flag ( PLUS | MINUS )                (Read) */

/* Local variables: */

/*  Tmp 	       Temporary (lives on even sites) */

/* Operations: */

/*  2 ( Nc Ns + DSlash + InvClvMs )  flops */
include(types.mh)

SUBROUTINE(lclovmpsim, A, psi, chi, Ncb, isign)

LINEAR_OPERATOR(A);
LatticeFermion psi;
LatticeFermion chi;
int Ncb;
int isign;
{ /* Local variables */
  include(COMMON_DECLARATIONS)

  LATTICE_TRIANG(clov);
  LATTICE_TRIANG(invclov);
  /* multi1d<LatticeColorMatrix> u(Nd); */
  Real Kappa_ds;

  LatticeFermion tmp1;
  LatticeFermion tmp2;
  Real dummy;
  int temp;

  LINEAR_OPERATOR(D);
  PROTOTYPE(`D', `DESC', `DESC', `DESC', `VAL', `VAL', `VAL')`'

  START_CODE();

  /* Extract arguments */

  UNPACK_LINEAR_OPERATOR(A, D, clov, invclov, Kappa_ds);

  /*  Chi   =  D'    A^(-1)     D'	   Psi  */
  /*     O      O,E        E,E   E,O          O */
      D (D, psi, tmp1, 1, isign, 1);
  ldumul (invclov, tmp1, tmp2);
  D (D, tmp2, tmp1, 1, isign, 0);
  
  /*      	       	          2 */
  /*  Chi   =  A    Psi  -  k  Chi  */
  /*     O      O,O    O          O */
  ldumul (clov, psi, chi);
  dummy = - (Kappa_ds*Kappa_ds);
  chi += tmp1 * dummy;       	       /* Nc Ns  flops */
    
  END_CODE();
}
