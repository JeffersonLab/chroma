/* $Id: ovdsduf_w.cc,v 3.0 2006-04-03 04:59:10 edwards Exp $ ($Date: 2006-04-03 04:59:10 $) */

/* This routine is specific to Wilson fermions! */

/* ovdsduf -- computes the derivative of the overlap pole fermionic action */
/* respect the link field */
/* u -- gauge field ( Read ) */

/*         |  dS      dS_f */
/* ds_u -- | ----   + -----   ( Modify ) */
/*         |  dU       dU */

/* psi -- [1./(M_dag*M)]*chi_  ( read ) */

void OvDsDuf(multi1d<LatticeColorMatrix>& ds_u,
	     const multi1d<LatticeColorMatrix>& u,
	     const LatticeFermion& psi)
{
  LINEAR_OPERATOR(A);
  LINEAR_OPERATOR(B);
  LINEAR_OPERATOR(MdagM);
  LINEAR_OPERATOR(M);
  int numroot;
  int NEig;
  multi1d<LatticeColorMatrix> utmp_p(Nd);
  LatticeColorMatrix utmp_t;
  multi1d<LatticeFermion> ftmp_1(Ncb); 
  multi1d<LatticeFermion> ftmp_2(Ncb); 
  multi1d<LatticeFermion> ftmp_3(Ncb); 
  multi2d<LatticeFermion> psin(Ncb, numroot);
  LatticeComplex lctmp1;
  LatticeComplex lctmp2;
  LatticeReal lrtmp1;
  LatticeReal lrtmp2;
  multi1d<Complex> ccoeff(NEig);
  Complex cconsts;
  DComplex dxp;
  multi1d<Real> resP(numroot);
  multi1d<Real> rootQ(numroot);
  Real dummyg;
  Real dummy1;
  Real dummy2;
  Double dummyd;
  multi2d<LatticeFermion> EigVec(Ncb, NEig);
  multi1d<Real> EigValFunc(NEig);
  Real constP;
  Real RsdCG;
  multi1d<Real> RsdCGV(numroot);
  Real m_q;
  Double psi_norm;
  int n_count;

  int s;
  int mu;
  int cb;
  int isz;
  int G5;
  int ichiral;
  int i;
  
  START_CODE("subroutine");;

            
  if ( FermAct != OVERLAP_POLE )
    QDP_error_exit("Wrong fermion action", FermAct);

  if ( KappaMC != KappaMD )
  {
    /* for the overlap we use for molecular dynamics (Kappa==KappaMD) the 
       Trunc overlap operator, but with Kappa=KappaMC.  */
    ConsLinOp (A, u, KappaMC, TRUNC_OVERLAP);
  }
  else
  {		
    ConsLinOp (A, u, KappaMD, FermAct);
  }

  UNPACK_LINEAR_OPERATOR(A, MdagM, M, m_q, 
			 numroot, constP, resP, rootQ, 
			 EigVec, EigValFunc, NEig, RsdCG);
    
  ftmp_1 = psi;

  if (NEig > 0)
  {
    /* Project eigenvectors out of ftmp_1 */
                
    for(i = 0; i < NEig; ++i)
    {
      lctmp1 = trace(adj[EigVec[i][0]] * ftmp_1[0]);
      lctmp2 = trace(adj[EigVec[i][1]] * ftmp_1[1]);
      lctmp1 += lctmp2;

      dxp = sum(lctmp1);
      ccoeff[i] = FLOAT(dxp);

      ftmp_1[0] -= EigVec[i][0] * ccoeff[i];
      ftmp_1[1] -= EigVec[i][1] * ccoeff[i];
    }

  }

  /* Getting the norm of (subtracted) psi for MinvCGm */
  psi_norm = norm2(ftmp_1[0]);
  psi_norm += norm2(ftmp_1[1]);
  psi_norm = sqrt(psi_norm); 

  FILL(RsdCGV, RsdCG);
  /* smallest mass, subject to change...*/
  isz = numroot - 1;
  /* calc all psi_n via multimass solver */

  MInvCGm (MdagM, ftmp_1, psin, psi_norm, rootQ, numroot, isz, RsdCGV, Ncb, n_count, EigVec, NEig);
  FPRINTF(trm_out,"In MInvCGm used CG steps: %d\n",n_count);

  if (NEig > 0)
  {
    /* Now add in contributions from eigenvectors */
    /* Note: EigValFunc does not contain the eigenvalues we need
       so we compute them quickly here. */
            
    for(i = 0; i < NEig; ++i)
    {
      for(cb = 0; cb < Ncb; ++cb )
	ftmp_1[cb] = EigVec[i][cb];
      MdagM (MdagM, ftmp_1, ftmp_2, Ncb, PLUS);

      lrtmp1 = real(trace(adj[ftmp_1[0]] * ftmp_2[0]));
      lrtmp2 = real(trace(adj[ftmp_1[1]] * ftmp_2[1]));
      lrtmp1 += lrtmp2;

      dummyd = sum(lrtmp1);
      dummy1 = FLOAT(dummyd);

      for(s = 0; s < numroot; ++s)
      {
	dummy2 = dummy1 + rootQ[s];
	/* cconsts = ccoeff[i] / dummy2; */
	dummy2 = TO_REAL(1) / dummy2;
	cconsts = ccoeff[i] * dummy2;
	psin[s][0] += EigVec[i][0] * cconsts;
	psin[s][1] += EigVec[i][1] * cconsts;
      }
    }

  }

  /* up to scalar factors                                 */
  /* dsdu+= psi^+ dH_wdU psi                              */
  /* calc CG(H_w+b_n,psi)=psi_n) all in ones with minvcgm */
  /* for n=0,N-1 do                          		  */           
  /*     dsdu+= psi_n^+ dH_w/dU psi_n      		  */
  /*     dsdu+= psi_n^+ H_w dH_w/dU H_w psi_n             */

  G5 = Ns*Ns-1;
  utmp_p = 0;

  /* Global factor for the action */
  /* including kappa of the hopping terms of the H_W operator */
  dummyg=TO_REAL(-0.5)*(TO_REAL(1)-m_q*m_q)
    *TO_REAL(0.5)/(TO_REAL(Nd)-OverMass);

  /* Calculate the chirality for the global sign */
  ischiral (psi, ichiral, Ncb);
  if (ichiral == 0)
  {
    QDP_error_exit("Psi must be chiral", ichiral);
  }
  else if(ichiral == -1)
    dummyg = - dummyg;

  dummy1=dummyg*constP;
  /* This copy is done only to have the same code as in the next section */
  for(cb = 0;cb  < ( Ncb ); ++cb )
  {
    ftmp_3[cb] = psi[cb];
  }
  for(mu = 0;mu  < ( Nd); ++mu )
  {
    /* ftmp_1 =  (gamma(mu))*psi */
    for(cb = 0; cb < Ncb; ++cb)
      SPIN_PRODUCT(ftmp_3(cb),(INTEGER_LSHIFT_FUNCTION(1,mu)),ftmp_1(cb));
    /* ftmp_2 = -(gamma(mu))*psi */
    ftmp_2 = -ftmp_1;
    /* ftmp_2 = (1 - gamma(mu))*psi */
    ftmp_2 += ftmp_3;

    /* ftmp_1 = gamma(5)*(1 - gamma(mu))*psi */
    for(cb = 0; cb < Ncb; ++cb)
      SPIN_PRODUCT(ftmp_2(cb),G5,ftmp_1(cb));
    
    for(cb = 0; cb < Ncb; ++cb)
    {
      /* utmp_t= diractrace[gamma(5)*(1-gamma(mu))psi(x+mu)psi(x)^dag ] */
      utmp_t = shift(ftmp_1[1-cb], cb, FORWARD, mu) * adj(ftmp_3[cb]);
      utmp_p[mu][cb] += utmp_t * dummy1;

      /* utmp_t= diractrace[gamma(5)*(1+gamma(mu))psi(x)psi(x+mu)^dag ] */
/* Not needed, since this will come automatically when doing the
   anti-hermitian projection! (UMH) */
    }
  }

  for(s = 0; s < numroot; ++s)
  {
    dummy1=  dummyg*resP[s]*rootQ[s];
    dummy2= -dummyg*resP[s];
    /* utmp_1+= psi_s dHdU psi_s */

    /* This copy is done only to have the same code as in the next section */
    for(cb = 0;cb  < ( Ncb ); ++cb )
    {
      ftmp_3[cb] = psin[s][cb];
    }
    for(mu = 0;mu  < ( Nd); ++mu )
    {
      /* ftmp_1 =  (gamma(mu))*psin */
      for(cb = 0; cb < Ncb; ++cb)
	SPIN_PRODUCT(ftmp_3(cb),(INTEGER_LSHIFT_FUNCTION(1,mu)),ftmp_1(cb));
      /* ftmp_2 = -(gamma(mu))*psin */
      ftmp_2 = -ftmp_1;
      /* ftmp_2 = (1 - gamma(mu))*psin */
      ftmp_2 += ftmp_3;

      /* ftmp_1 = gamma(5)*(1 - gamma(mu))*psin */
      for(cb = 0; cb < Ncb; ++cb)
	SPIN_PRODUCT(ftmp_2(cb),G5,ftmp_1(cb));
    
      for(cb = 0;cb  < ( Ncb ); ++cb )
      {
	/* utmp_t= diractrace[gamma(5)*(1-gamma(mu))psin(x+mu)psin(x)^dag ] */
	utmp_t = shift(ftmp_1[1-cb], cb, FORWARD, mu) * adj(ftmp_3[cb]);
	utmp_p[mu][cb] += utmp_t * dummy1;
      }
    }
    /* utmp_1+= psi_s H dHdU H psi_s */
    /* Maybe the opertator H_w is out anywhere ! */

    for(cb = 0;cb  < ( Ncb ); ++cb )
    {
      ftmp_3[cb] = psin[s][cb];
    }
    M (M, ftmp_3, ftmp_1, Ncb, PLUS);

    for(cb = 0; cb < Ncb; ++cb)
      SPIN_PRODUCT(ftmp_1(cb),G5,ftmp_3(cb));

    for(mu = 0;mu  < ( Nd); ++mu )
    {
      /* ftmp_1 =  (gamma(mu))*H*psin */
      for(cb = 0; cb < Ncb; ++cb)
	SPIN_PRODUCT(ftmp_3(cb),(INTEGER_LSHIFT_FUNCTION(1,mu)),ftmp_1(cb));
      /* ftmp_2 = -(gamma(mu))*H*psin */
      ftmp_2 = -ftmp_1;
      /* ftmp_2 = (1 - gamma(mu))*H*psin */
      ftmp_2 += ftmp_3;

      /* ftmp_1 = gamma(5)*(1 - gamma(mu))*H*psin */
      for(cb = 0; cb < Ncb; ++cb)
	SPIN_PRODUCT(ftmp_2(cb),G5,ftmp_1(cb));
    
      for(cb = 0;cb  < ( Ncb ); ++cb )
      {
	/* utmp_t= diractrace[gamma(5)*(1-gamma(mu))H*psin(x+mu)H*psin(x)^dag ] */
	utmp_t = shift(ftmp_1[1-cb], cb, FORWARD, mu) * adj(ftmp_3[cb]);
	utmp_p[mu][cb] += utmp_t * dummy2;
      }
    }
  }

  for(mu = 0;mu  < ( Nd); ++mu )
  {
    for(cb = 0;cb  < ( Ncb ); ++cb )
    { 
      ds_u[mu][cb] += u[mu][cb] * utmp_p[mu][cb];
    }
  }

  DestLinOp (A, FermAct);
            
  END_CODE("subroutine");;
}
