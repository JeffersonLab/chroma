// $Id: eig_spec_w.cc,v 1.1 2004-01-04 21:56:04 edwards Exp $
/*! \file
 *  \brief Compute low lying eigenvalues of the hermitian H = gamma_5 A 
 */

#error "CONVERSION NOT COMPLETE: NEED TO MAKE APPROPRIATE SUBCLASSING"

#include "chromabase.h"
#include "meas/eig/ritz.h"

using namespace QDP;

//! Compute low lying eigenvalues of the hermitian H = gamma_5 A
/*!
 * \ingroup eig
 *
 *  Compute low lying eigenvalues of the hermitian H = gamma_5 A
 *  using the Ritz functional minimization routine,
 *  if desired in the Kalkreuter-Simma algorithm

 *  A		The operator A as a lin op		(Read)
 *  B		The "operator square"		(Read)
 *  psi		Eigenvectors			(Modify)
 *  N_eig	Eigenvalue number		(Read)
 *  Ncb		Number of sublattices		(Read)
 *  lambda_H	Eigenvalues			(Write)
 *  valid_eig	flag for valid eigenvalues	(Write)
 *  IOper	Int flag for relation between A and B	(Read)
 *  N_min	Minimal CG iterations		(Read)
 *  N_max	Maximal CG iterations		(Read)
 *  Kalk_Sim	Use Kalkreuter-Simma criterion	(Read)
 *  n_renorm	Renormalize every n_renorm iter.	(Read)
 *  N_KS_max	Max number of Kalkreuter-Simma iterations	(Read)
 *  RsdR_a	(absolute) residue		(Read)
 *  RsdR_r	(relative) residue		(Read)
 *  Rsdlam	relative accuracy of lambda	(Read)
 *  Cv_fact	"Convergence factor" required	(Read)
 *  NCG_tot	Total number of CG iter		(Write)
 *  n_KS	total number of Kalkreuter-Simma iterations	(Write)
 *  n_valid	number of valid eigenvalues	(Write)
 *  n_zero	number of zero eigenvalues	(Write)
 *  ProjApsiP	flag for projecting A.psi	(Read)
 *  SU2AdjP	flag for SU(2) adjoint repr	(Read)
 *  OverlapP	flag for overlap fermions	(Read) 
 */


void Eig_Spec(A, B, psi, N_eig, Ncb, lambda_H, valid_eig, IOper,
	      N_min, N_max, Kalk_Sim, n_renorm, N_KS_max,
	      RsdR_a, RsdR_r, Rsdlam, Cv_fact,
	      NCG_tot, n_KS, n_valid, n_zero, ProjApsiP, SU2AdjP, OverlapP)
  LINEAR_OPERATOR(A);
LINEAR_OPERATOR(B);
multi2d<LatticeFermion> psi(Ncb, N_eig);
multi1d<Real> lambda_H(N_eig);
Real RsdR_a;
Real RsdR_r;
Real Rsdlam;
Real Cv_fact;
multi1d<int> valid_eig(N_eig);
int N_eig;
int Ncb;
int IOper;
int N_min;
int N_max;
int Kalk_Sim;
int n_renorm;
int N_KS_max;
int NCG_tot;
int n_KS;
int n_valid;
int n_zero;
bool ProjApsiP;
bool SU2AdjP;
bool OverlapP;
{
  START_CODE("Eig_spec");
  
  LINEAR_OPERATOR(C);
  LINEAR_OPERATOR(M);

  multi1d<LatticeFermion> psi_t(Ncb);
  multi1d<LatticeFermion> tmp(Ncb);
  LatticeFermion tmp3;
  LatticeSpinMatrix fe_tmp;
  multi1d<Complex> off_diag(N_eig*(N_eig-1)/2);
  Complex eigen;
  Double ddummy;
  Double ddumm2;
  multi1d<Real> lambda(N_eig);
  multi1d<Real> lambda_old(N_eig);
  multi1d<Real> lambda_cck(N_eig);
  Real dummy;
  Real lambda_t;
  Real lambda_t2;
  Real chiral;
  Real del_lamb;
  Real acc;
  Real rsdl_sq;
  Real rsdl_zero;
  int i;
  int n;
  int j;
  int ij;
  int cb;
  int cc;
  int sc;
  int toggle;
  int n_ritz;
  int n_jacob;
  bool OvlapChiralP;
  int ichiral;
  
  int G5 = Ns*Ns-1;

  if (N_eig == 1)
    Kalk_Sim = NO;

  if (! Kalk_Sim)
    N_KS_max = 0;

  bool NonConv;

  /* Determine machine accuracy */
  FILL(acc, FUZZ);
  acc /= Real(5);
  acc *= acc;
  rsdl_sq = Rsdlam * Rsdlam;
  rsdl_sq = max(rsdl_sq, acc);
  rsdl_zero = Real(10) * rsdl_sq;

  if (N_eig > 1)
    
    NCG_tot = 0;
  n_zero = 0;

  if (OverlapP && ((FermAct == OVERLAP_POLE) || (FermAct == ZOLOTAREV_4D)) && Kalk_Sim == YES)
  {
    UNPACK_LINEAR_OPERATOR(A, M);

    ichiral = PLUS;
    CONSTRUCT_LINEAR_OPERATOR(C, lovddag, M, ichiral);
    OvlapChiralP = YES;

    /* Hardwire to some overlap specific values, since accuracy */
    /* of H^2 and H are very similar here. */
    rsdl_sq = Real(5.0e-7);
    rsdl_zero = Real(5.0e-6);
  }
  else
  {
    C = B;
    OvlapChiralP = NO;
  }

  if (Kalk_Sim == YES)
  {
    lambda_old = 1;
    off_diag = 0;
    del_lamb = 1;

    if (OvlapChiralP == YES)
    {
      /* First try one positive and one negative chiral vector */
      n = 1;
      dummy = Rsdlam;
      Rsdlam *= Real(100);
      del_lamb = sqrt(rsdl_sq);
      del_lamb = max(del_lamb, Rsdlam);
      psi_t = Gamma(G5)*psi[0] + psi[0];
      Ritz(C, psi_t, n, lambda_t, del_lamb, Rsdlam, Rsdlam, rsdl_sq, 
	   n_renorm, Ncb, n_ritz, NO, N_min, 0, Cv_fact, ProjApsiP);
      NCG_tot += n_ritz;

      FREE_LINEAR_OPERATOR(C);
      ichiral = MINUS;
      CONSTRUCT_LINEAR_OPERATOR(C, lovddag, M, ichiral);

      tmp = Gamma(G5)*psi[1] - psi[1];
      Ritz(C, tmp, n, lambda_t2, del_lamb, Rsdlam, Rsdlam, rsdl_sq, 
	   n_renorm, Ncb, n_ritz, NO, N_min, 0, Cv_fact, ProjApsiP);
      NCG_tot += n_ritz;

      push(xml_out,"Trial_Ritzes");
      Write(xml_out, NCG_tot);
      Write(xml_out, n_ritz);
      Write(xml_out, lambda_t);
      Write(xml_out, lambda_t2);
      pop(xml_out);
      Rsdlam = dummy;

      if (lambda_t2 < lambda_t)
      {
	QDPIO::cout << "Do Kalkreuter-Simma for negative chirality" << endl;
	/* Do negative chirality Kalkreuter-Simma */
	lambda[0] = lambda_t2;
	lambda[1] = lambda_t;
	lambda_t2 = lambda_t;
	/* Keep tmp, and construct negative chirality vector from psi_t */
	psi[0] = tmp;
	A(tmp, psi_t, PLUS);
	psi_t = Gamma(G5)*tmp;
	tmp -= psi_t;
	psi[1] = tmp;
	
	/* Now make the other random vectors negative chiral */
	for(i = 2; i < N_eig; i++)
	  psi_t = Gamma(G5)*psi[i] - psi[i];
      }
      else
      {
	QDPIO::cout << "Do Kalkreuter-Simma for positive chirality" << endl;
	/* Do positive chirality Kalkreuter-Simma */
	lambda[0] = lambda_t;
	lambda[1] = lambda_t2;
	/* Keep psi_t, and construct positive chirality vector from tmp */
	psi[0] = psi_t;
	A(psi_t, tmp, PLUS);
	tmp = Gamma(G5)*psi_t + psi_t;
	psi[1] = tmp;
	
	/* Construct positive chirality operator, again */
	FREE_LINEAR_OPERATOR(C);
	ichiral = PLUS;
	CONSTRUCT_LINEAR_OPERATOR(C, lovddag, M, ichiral);

	/* Now make the other random vectors positive chiral */
	for(i = 2; i < N_eig; i++)
	{
	  psi_t = Gamma(G5) * psi[i];
	  psi[i] += psi_t;
	}
      }
    }

    for(n_KS = 1; n_KS <= N_KS_max && del_lamb > Rsdlam; n_KS++)
    {
      ij = 0;
      toggle = false;
      for(n = 1, i = 0; n <= N_eig; n++, i++)
      {
	if (! SU2AdjP || i%2 == 0)
        {
	  if (OverlapP && toggle && ! OvlapChiralP)
	  {
	    /* New trial vector is (gamma_5*prev_vector) */
	    psi_t = Gamma(G5)*psi[i-1];
	    psi[i] = psi_t;
	    /* Note orthogonalization and normalization are done in Ritz */
	  }

	  Ritz(C, psi, n, lambda_t, RsdR_a, RsdR_r, Rsdlam, rsdl_sq, 
	       n_renorm, Ncb, n_ritz, Kalk_Sim, N_min, N_max, Cv_fact, ProjApsiP);
	  NCG_tot += n_ritz;

	  lambda[i] = lambda_t;
	  psi_t = psi;

	  if (OverlapP && ! OvlapChiralP)
	  {
	    /* Test if mode looks sufficiently chiral */
	    	    	    	    	    
	    tmp = Gamma(G5)*psi_t;
	    chiral = fabs(innerProductReal(psi_t, tmp));

	    if (i == n_zero && chiral > 0.98)
	      n_zero++;

	    if (! toggle && chiral < 0.98 && n_KS < 5)
	    {
	      toggle = true;
	    }
	    else if(toggle)
	      toggle = false;

	    QDPIO::cout << " chiral-test: i=" << k << ", chiral=" << chiral
			<< ", toggle=" << toggle << endl;
	  }
	}
	else
	{
	  if (SU2AdjP)
	  {
	    /* Adjoint SU(2): 2nd of the pair of degenerate eigenvalues */
	    lambda[i] = lambda[i-1];
	    n_ritz = 0;

	    for(cb = 0; cb < Ncb; cb++)
	    {
	      tmp3 = psi[i-1][cb];
	      fe_tmp = CAST(tmp3);
	      for(cc = 0; cc < Nc; cc++)
		for(sc = 0; sc < Ns; sc++)
		  fe_tmp[sc][cc] = adj(fe_tmp[sc][cc]);
	      tmp3 = CAST(fe_tmp);
	      tmp[cb] = tmp3;
	    }
	    psi_t = Gamma(5) * tmp;
	    psi[i] = psi_t;
	  }
	  else
	  {
	    QDP_error_exit("Internal error - SU2AdjP and OverlapP not set\n");
	  }
	}

					
	C(tmp, psi_t, PLUS);

	for(j = 0; j < i; j++)
	{
	  off_diag[ij] = innerProduct(psi[j], tmp);
	  ij++;
	}
      }

      SN_Jacob(psi, N_eig, lambda, off_diag, RsdR_a, 50, Ncb, n_jacob);
      push(xml_out,"MatElem_after");
      Write(xml_out, n_jacob);
      Write(xml_out, lambda);
      pop(xml_out);
      if (OverlapP && ! OvlapChiralP)
	QDPIO::cout << " at n_KS " << n_ks << ": n_zero = " << n_zero << endl;

      lambda_old[0] -= lambda[0];
      lambda_t = fabs(lambda_old[0]);
      if (lambda_t > rsdl_sq)
      {
	del_lamb = fabs(lambda_t / lambda[0]);
      }
      else
	del_lamb = 0;
      lambda_old[0] = lambda[0];
      for(i = 1; i < N_eig; i++)
      {
	lambda_old[i] -= lambda[i];
	lambda_t = fabs(lambda_old[i]);
	if (lambda_t > rsdl_sq)
	{
	  lambda_t /= lambda[i];
	  lambda_t = fabs(lambda_t);
	  del_lamb = max(del_lamb, lambda_t);
	}
	lambda_old[i] = lambda[i];
      }
    }

    if (n_KS >= N_KS_max)
    {
      NonConv = true;
      push(xml_out,"NONConversion_warning");
      Write(xml_out, n_KS);
      Write(xml_out, N_KS_max);
      pop(xml_out);
    }

    if (OvlapChiralP)
    {
                              
      del_lamb = sqrt(Rsdlam);
      n_valid = 0;
      valid_eig = 0;

      if (ichiral == PLUS)
      {
	/* How many positive chirality zero eigenvalues? */
	for(i = 0; i < N_eig; i++)
	{
	  psi_t = psi[i];

	  A(tmp, psi_t, PLUS);
	  lambda_t = innerProductReal(psi_t, tmp);

	  psi_t = Gamma(G5) * tmp;
	  tmp -= psi_t;
	  ddummy = norm2(tmp);
	  ddumm2 = 0.5*sqrt(1 - ddummy);
	  if (lambda[i] < Real(0.5))
	  {
	    lambda_cck[i] = 0.5 - ddumm2;
	  }
	  else
	  {
	    lambda_cck[i] = 0.5 + ddumm2;
	  }

	  if (lambda_cck[i] < rsdl_zero)
	  {
	    n_zero++;
	    lambda_H[i] = lambda_t;
	    lambda_t = fabs(lambda_t);
	    if (lambda_t < Rsdlam)
	    {
	      valid_eig[i] = 1;
	      n_valid++;
	    }
	  }
	  else
	  {
	    lambda_H[i] = sqrt(lambda_t);
	    lambda_t -= lambda_cck[i];
	    lambda_t = fabs(lambda_t);
	    lambda_t = lambda_t / lambda[i];
	    if (lambda_t < del_lamb)
	    {
	      valid_eig[i] = 1;
	      n_valid++;
	    }

	    if (i == n_zero && lambda_t2 < lambda[i])
	    {
	      lambda_t = lambda[i];
	      push(xml_out,"Zero_identification_warning");
	      Write(xml_out, i);
	      Write(xml_out, lambda_t);
	      Write(xml_out, lambda_t2);
	      pop(xml_out);
	      QDPIO::cout << "Warning: lambda_t2 = " << lambda_t2
			  << " < lambda(" << i << ") = " << lambda_t << endl;
	    }

	    dummy = Real(1) / sqrt(ddummy);
	    /* Pack both chirality eigenvectors into same Dirac spinor */
	    for(cb = 0; cb < Ncb; cb++)
	      psi[i][cb] += tmp[cb] * dummy;
	  }
	}

	push(xml_out,"Positive_Chirality");
	Write(xml_out, n_zero);
	Write(xml_out, n_KS);
	pop(xml_out);
      }
      else	/* ichiral == MINUS */
      {
	/* How many negative chirality zero eigenvalues? */
	for(i = 0; i < N_eig; i++)
	{
	  psi_t = psi[i];
	  A(tmp, psi_t, PLUS);
	  lambda_t = innerProductReal(psi_t, tmp);

	  psi_t = Gamma(G5)*tmp;
	  tmp += psi_t;
	  ddummy = norm2(tmp);
	  ddumm2 = 0.5*sqrt(1 - ddummy);
	  if (lambda[i] < Real(0.5))
	  {
	    lambda_cck[i] = 0.5 - ddumm2;
	  }
	  else
	  {
	    lambda_cck[i] = 0.5 + ddumm2;
	  }

	  if (lambda_cck[i] < rsdl_zero)
	  {
	    n_zero++;
	    lambda_H[i] = -lambda_t;
	    lambda_t = fabs(lambda_t);
	    if (lambda_t < Rsdlam)
	    {
	      valid_eig[i] = 1;
	      n_valid++;
	    }
	  }
	  else
	  {
	    lambda_t = - lambda_t;
	    lambda_H[i] = sqrt(lambda_t);
	    lambda_t -= lambda_cck[i];
	    lambda_t = fabs(lambda_t);
	    lambda_t = lambda_t / lambda[i];
	    if (lambda_t < del_lamb)
	    {
	      valid_eig[i] = 1;
	      n_valid++;
	    }

	    if (i == n_zero && lambda_t2 < lambda[i])
	    {
	      lambda_t = lambda[i];
	      push(xml_out,"Zero_identification_warning");
	      Write(xml_out, i);
	      Write(xml_out, lambda_t);
	      Write(xml_out, lambda_t2);
	      pop(xml_out);
	      FPRINTF(trm_out,"Warning: lambda_t2 = %g < lambda(%d) = %g\n",
		      lambda_t2, i, lambda_t);
	    }

	    dummy = Real(1) / sqrt(ddummy);
	    /* Pack both chirality eigenvectors into same Dirac spinor */
	    for(cb = 0; cb < Ncb; cb++)
	      psi[i][cb] += tmp[cb] * dummy;
	  }
	}

	push(xml_out,"Negative_Chirality");
	Write(xml_out, n_zero);
	Write(xml_out, n_KS);
	pop(xml_out);
	n_zero = - n_zero;
      }

                              
      push(xml_out,"Final_Chiral_Overlap");
      Write(xml_out, n_zero);
      Write(xml_out, n_valid);
      Write(xml_out, lambda_H);
      Write(xml_out, lambda);
      Write(xml_out, lambda_cck);
      Write(xml_out, valid_eig);
      pop(xml_out);

    }		/* End of OvlapChiralP == YES */
    else if (IOper < 2)
    {
      /* Prepare to eigenvalues of the operator, not only of the "square" */
                              
      ij = 0;
      for(i = 0; i < N_eig; i++)
      {
	psi_t = psi[i];

	A(tmp, psi_t, PLUS);
	lambda_H[i] = innerProductReal(psi_t, tmp);
	lambda_old[i] = fabs(lambda_H[i]);

	for(j = 0; j < i; j++)
	{
	  off_diag[ij] = innerProduct(psi[j], tmp);
	  ij++;
	}
      }

    }
  }
  else	/* Kalk_Sim == NO */
  {
    if (N_eig > 1)
      off_diag = 0;

    ij = 0;
    for(n = 1, i = 0; n <= N_eig; n++, i++)
    {
      if ((! SU2AdjP && (OverlapP == NO || FermAct == TRUNC_OVERLAP)) ||
	  i <= n_zero || (i+n_zero)%2 == 0)
      {
	Ritz (B, psi, n, lambda_t, RsdR_a, RsdR_r, Rsdlam, rsdl_sq, 
	      n_renorm, Ncb, n_ritz, NO, 0, 0, Cv_fact, ProjApsiP);
	NCG_tot += n_ritz;
	lambda[i] = lambda_t;

	if (OverlapP == YES && i <= n_zero)
	{
	  /* Test if mode looks sufficiently chiral */
	  tmp = Gamma(G5) * psi;
	  chiral = fabs(innerProductReal(psi[i], tmp));
	  QDPIO::cout << " chiral-test i=" << i << "  " << chiral << endl;
	  if (i == n_zero && chiral > 0.999)
	    n_zero++;

	}
      }
      else
      {
	if (SU2AdjP)
	{
	  /* Adjoint SU(2): 2nd of the pair of degenerate eigenvalues */
	  lambda[i] = lambda[i-1];
	  n_ritz = 0;

	  for(cb = 0; cb < Ncb; cb++)
	  {
	    tmp3 = psi[i-1][cb];
	    fe_tmp = CAST(tmp3);
	    for(cc = 0; cc < Nc; cc++)
	      for(sc = 0; sc < Ns; sc++)
		fe_tmp[sc][cc] = adj(fe_tmp[sc][cc]);
	    tmp3 = CAST(fe_tmp);
	    tmp[cb] = tmp3;
	  }
	  psi_t = Gamma(5) * tmp;
	  psi[i] = psi_t;
	}
	else if (OverlapP == YES && FermAct != TRUNC_OVERLAP)
	{
	  lambda[i] = lambda[i-1];
	  n_ritz = 0;

	  /* New vector is (gamma_5*prev_vector) orthog wrt to previous vects */
	  psi_t = Gamma(G5) * psi[i-1];

	  GramSchm(psi_t, 1, psi, i, Ncb);

	  /* Normalize */
	  dummy = Real(1) / sqrt(norm2(psi_t));
	  for(cb = 0; cb < Ncb; ++cb)
	  {
	    psi_t[cb] *= dummy;
	    psi[i][cb] = psi_t[cb];
	  }
	}
	else
	{
	  QDP_error_exit("Internal error - SU2AdjP and OverlapP not set\n");
	}
      }

      if (IOper < 2)
      {
	/* Prepare to eigenvalues of the operator, not only of the "square" */
	psi_t = psi[i];

	A(tmp, psi_t, PLUS);
	lambda_H[i] = innerProduct(psi_t, tmp);
	lambda_old[i] = fabs(lambda_H[i]);

	for(j = 0; j < i; j++)
	{
	  off_diag[ij] = innerProduct(psi[j], tmp);
	  ij++;
	}

      }
    }

    if (OverlapP)
      QDPIO::cout << " n_zero = " << n_zero << endl;
  }

  if (! OvlapChiralP)
  {
    if (IOper < 2)
    {
      /* Do final Jacoby for ev's of A when B=A^2 */
      n = N_eig;

      dummy = sqrt(Rsdlam);

      if (N_eig > 1)
      {
	SN_Jacob (psi, n, lambda_H, off_diag, RsdR_a, 50, Ncb, n_jacob);

	/* Label eigenvalues okay, or not */
	n_valid = 0;
	for(n = 0; n < N_eig; n++)
	{
	  valid_eig[n] = 0;
	  lambda_t = lambda_H[n] * lambda_H[n];
	  if( lambda_t < rsdl_zero )
	  {
	    for(j = n_valid; j < N_eig; j++)
	      if( lambda[j] < rsdl_zero )
	      {
		valid_eig[n] = 1;
		n_valid++;
		break;
	      }
	  }
	  else
	  {
	    for(j = n_valid; j < N_eig; j++)
	    {
	      del_lamb = lambda[j];
	      del_lamb -= lambda_t;
	      del_lamb = fabs(del_lamb);
	      del_lamb = del_lamb / lambda[j];
	      if( del_lamb < dummy )
	      {
		valid_eig[n] = 1;
		n_valid++;
		break;
	      }
	    }
	  }
	}

	push(xml_out,"Final_jacobi");
	Write(xml_out, n_jacob);
	Write(xml_out, n_valid);
	Write(xml_out, lambda_old);
	Write(xml_out, lambda_H);
	Write(xml_out, valid_eig);
	pop(xml_out);

      }
      else	/* N_eig = 1 */
      {
	del_lamb = lambda_H[0] * lambda_H[0];
	del_lamb -= lambda[0];
	del_lamb = fabs(del_lamb);
	del_lamb = del_lamb / lambda[0];
	if( del_lamb < dummy )
	{
	  valid_eig[0] = 1;
	  n_valid = 1;
	}
	else
	{
	  valid_eig[0] = 0;
	  n_valid = 0;
	}
	push(xml_out,"Final_result");
	Write(xml_out, n_valid);
	Write(xml_out, lambda);
	Write(xml_out, lambda_H);
	Write(xml_out, valid_eig);
	pop(xml_out);
      }
    }
    else
    {
      /* A and B are identical */
      n_valid = N_eig;
      ij = -1;
      FILL(valid_eig, ij);
      lambda_H = lambda;
      push(xml_out,"Final_result_sq");
      Write(xml_out, lambda);
      pop(xml_out);

      if (N_eig > 1)
	}
  }
  else	/* OvlapChiralP == YES */
  {
    FREE_LINEAR_OPERATOR(C);
  }

  END_CODE("Eig_spec");
}
