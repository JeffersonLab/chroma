// $Id: pade_trln_w.cc,v 3.0 2006-04-03 04:59:06 edwards Exp $

#error "NOT FULLY CONVERTED"


include(types.mh)


/* Parameters:
   
   u             -- The Gauge Field         (Read)
   chi           -- Z2 noise vector         (Read)
   Ncb           -- No of checkerboards     (Read)
   Nkappa        -- The no of kappa values  (Read)
   Npade         -- The order of the Pade approx (Read)
   Kappas        -- The array of Kappa values (Read)
   Pade_c        -- Array holding the denominator Pade coeffs  (Read)
   Pade_b0       -- b0 from the Pade Expansion (Read)
   Pade_b        -- Array holding numerator Pade coeffs (Read)
   RsdMR         -- Stopping Criterion for the MR algorithm    (Read)
   psibp         -- An array to hold the psibar-psi values     (Write)
   trlns         -- An array to hold the trace log values      (Write)
   n_count       -- The number of MR iterations taken          (n_count)
*/
SUBROUTINE(PadeTrLn, u, chi, Ncb, Nkappa, Kappas, Npade,
	   Pade_b0, Pade_b, Pade_c, Nimpr, Omega_impr,
	   RsdMR, psibp_unimp, trlns_unimp, psibp_imp, trlns_imp,
	   n_count)

     /* Gauge field */
     multi1d<LatticeColorMatrix> u(Nd);

     /* Z2 src */
     multi1d<LatticeFermion> chi(Ncb);
     
     /* No of checkerboards */
     int Ncb;
     
     /* No of Kappa values */
     int Nkappa;
    
     /* No of Pade Coeffs */
     int Npade;

     /* The array of Kappas */
     multi1d<Real> Kappas(Nkappa);

     /* Index of largest kappa (smallest mass) */
     /* The array of Pade_c */
     multi1d<Real> Pade_c(Npade);
     Real Pade_b0;
     multi1d<Real> Pade_b(Npade);

     /* Stopping criterion for the Multi Mass MR solver */
     Real RsdMR;
     
     /* Array for unimproved psibar psi estimates */
     multi1d<Complex> psibp_unimp(Nkappa);

     /* Array of unimproved Trace Log estimates */
     multi1d<Complex> trlns_unimp(Nkappa);

    /* Array for improved psibar psi estimates */
     multi1d<Complex> psibp_imp(Nkappa);

     /* Array of improved Trace Log estimates */
     multi1d<Complex> trlns_imp(Nkappa);

     /* No of solver iters taken */
     int n_count;
     
     /* Order of unbiased subtraction */
     int Nimpr;
     multi1d<Real> Omega_impr(Nimpr);

{
  include(COMMON_DECLARATIONS)

  multi2d<LatticeFermion> psi(Ncb, Nshift);   /* Array of solver solutions */
  
  multi1d<LatticeFermion> psi_test(Ncb); /* SOmething to check the residuals with */
  multi1d<LatticeFermion> r(Ncb);

  multi1d<LatticeFermion> Dchi(Ncb);

  Double chi_norm;               /* The norm of the source */
  Double r_norm;                 /* The norm of the residuals */

  /* The shifted and unshifted operators */

  /* A dslash operator */
  LINEAR_OPERATOR(D);

  LINEAR_OPERATOR(Mshifted);
  PROTOTYPE(`Mshifted', `DESC', `DESC', `DESC', `VAL', `VAL')

  LINEAR_OPERATOR(Munshifted);
  PROTOTYPE(`Munshifted', `DESC', `DESC', `DESC', `VAL', `VAL')

  /* The number of shifted solutions */
  int Nshift;             /* No of shifts for solver */
  multi1d<Real> Shifts(Nshift);        /* The array of shifts */

  multi1d<Real> RsdMR_array(Nshift);   /* An array version of RsdMR needed by solver */

  Real Shiftmass;             /* Useful for testing the solutions */

  Real MaxRsd;
  /* Various Counters */
  int i;
  int j;
  int l;
  int kappa;
  int cb;

   /* Two variables for tracing */
  LatticeComplex lc;             
  LatticeComplex lct; 
  DComplex dcsum;
  Complex csum;

  /* Temporaries */
  multi1d<Real> first_terms(Nkappa);
  multi1d<Complex> pade_terms(Npade);
		    
  
  /* Stuff for Trln */
  Real TrI;
  Real minusone;
  Real tmp;
  Real kappa_pow;
  Real pade_c_pow;


  /* Dslash Bilinears */
  multi1d<Complex> bilinear(Nimpr);

  /* Psibar Psi improvement terms */
  multi2d<Complex> psib_offterms(Nkappa, Nimpr);
  multi3d<Complex> trln_offterms(Nkappa, Npade, Nimpr);
  multi2d<Complex> sum_trln_offterms(Nkappa, Nimpr);
  
  /* Glue improvement : 'Plaquette' */
  Double w_plaq;
  Double s_plaq;
  Double t_plaq;
  Double link;
  Real f_plaq;
  Complex c_plaq;

  /* index of the smallest shift */
  int isz;
  Real minshift;
  
  START_CODE();
  
  minusone = - TO_REAL(1);

  /* First work out the shifts from the kappas and the sigmas */
  Nshift = Nkappa * ( Npade + 1 );

  /* Allocate the space */
      
  /* Work out the shifts */
  j=0;
  
  /* First just the 1/kappa's */
  for(i = 0; i < Nkappa; i++) {
    Shifts[j] = TO_REAL(1) / Kappas[i]; 
    j++;
  }

  /* Then the (1 + shift)/kappa's */
  for(i = 0; i < Nkappa; i++) { 
    for(l = 0; l < Npade; l++) { 
      Shifts[j] = ( TO_REAL(1) + Pade_c[l] ) / Kappas[i];
      j++;
    }
  }
  
  /* Set up the solver */
    
  FILL(RsdMR_array, RsdMR);
  /* push(xml_out,"pade_trln");
write(xml_out, "Nshift", Nshift);
write(xml_out, "Shifts", Shifts);
write(xml_out, "RsdMR", RsdMR);
write(xml_out, "MaxCG", MaxCG);
pop(xml_out); */

  /* Fill all Psi-s with 0 to start with */
  psi = 0;

  /* Compute the norm of chi */
  chi_norm=TO_REAL(0);

  for(cb=0; cb < Ncb; cb++) {
    chi_norm += norm2(chi[cb]);
  }
  chi_norm = sqrt(chi_norm);

  /* Flick on fermionic Bc's */
  phfctr (u, FORWARD);


    Shiftmass = 0;
  ConsDslash (D, u, Shiftmass, WILSON_DSLASH);

    CONSTRUCT_LINEAR_OPERATOR(Munshifted, ldunshifted, D);

  /* Work out what the smallest shift */
  minshift = Shifts[0];
  isz = 0;
  for(i=1; i < Nshift; i++) { 
    if(  Shifts[i] < minshift ) { 
      minshift = Shifts[i];
      isz = i;
    }
  }

  /* Call the solver */
  MInvMRm (Munshifted, chi, psi, chi_norm, Shifts, Nshift, isz, RsdMR_array, Ncb, n_count);

  /* push(xml_out,"Pade_MR_Solver");
write(xml_out, "n_count", n_count);
pop(xml_out); */

  FREE_LINEAR_OPERATOR(Munshifted);
  
        
  for(i=0; i < Nshift; i++) { 
    
    for(cb = 0; cb < Ncb; cb++) { 
      psi_test[cb] = psi[i][cb];
    }
   

    Shiftmass = Shifts[i];
    CONSTRUCT_LINEAR_OPERATOR(Mshifted, ldshifted, D, Shiftmass);

    Mshifted (Mshifted, psi_test, r, Ncb, PLUS);

    FREE_LINEAR_OPERATOR(Mshifted);
    
    Shiftmass = - TO_REAL(1);
    for(cb = 0; cb < Ncb; cb++) { 
      r[cb] = r[cb] * Shiftmass;
      r[cb] = chi[cb] + r[cb];
    }


    r_norm = 0;
    for(cb = 0; cb < Ncb; cb++) { 
      r_norm += norm2(r[cb]);
    }

    r_norm = sqrt(r_norm);
    r_norm /= chi_norm;
    RsdMR_array[i] = FLOAT(r_norm); 
  }

  MaxRsd = RsdMR_array[0];
  for( i=1; i < Nshift; i++) { 
    if( RsdMR_array[i] > MaxRsd ) { 
      MaxRsd  = RsdMR_array[i];
    }
  }
  PRINTF("TrLnSolve: n_count = %d, RsdMR = %16.8e, MaxRsd = %16.8e\n", 
	 n_count, RsdMR, MaxRsd);

    /* Flip off Boundary Conditions */
  phfctr (u, BACKWARD);


  /* Now make sure we have all converged */
  /* push(xml_out,"pade_trln_rsd");
write(xml_out, "RsdMR_array", RsdMR_array);
pop(xml_out); */

  

  for(i=0; i < Nkappa; i++) { 
    tmp = TO_REAL(1) / Kappas[i];
    for(cb = 0; cb < Ncb; cb++) { 
      psi[i][cb] = psi[i][cb] * tmp;
    }
  }

  j=Nkappa; 
  for(i=0; i < Nkappa; i++) { 
    tmp = TO_REAL(1) / Kappas[i];
    for(l=0; l < Npade; l++) { 
      for(cb = 0; cb < Ncb; cb++) {
	psi[j][cb] = psi[j][cb] * tmp; 
      }
      j++;
    }
  }

  /* Work out dumb psibar psi estimates */
      

  for(i = 0; i < Nkappa; i++) { 

    lc = trace(adj[chi[0]] * psi[i][0]);
   
    for(cb = 1; cb < Ncb; cb++) { 
      lct = trace(adj[chi[cb]] * psi[i][cb]);
      lc += lct;
    }
 
    dcsum = sum(lc);
   
    psibp_unimp[i] = FLOAT(dcsum);
  
  }
  

  /* Now try to compute TrLn = b0 TrI - sum Pade_c(l) chi (M + c(l))^{-1} chi */

  /* First find the first term: Pade_b0 * Tr I, where Tr I = Nc * Nd * V */
    TrI = TO_REAL(Ns*Nc)*TO_REAL(2)*TO_REAL(vol_cb);
  trlns_unimp = 1;
  tmp = TrI * Pade_b0;
  FILL(first_terms, tmp);
  trlns_unimp = first_terms * trlns_unimp;
  
  /* Now get the remaining Npade terms */
  
  j = Nkappa;
  for(i = 0 ; i < Nkappa; i++) {
    for(l = 0; l < Npade; l++) {

      /* pade_diag_term(l) = < chi, M(kappa, pade_c) chi > */
      lc = trace(adj[chi[0]] * psi[j][0]);
      for(cb = 1; cb < Ncb; cb++) {
	lct = trace(adj[chi[cb]] * psi[j][cb]);
	lc += lct;
      }

      /* Sum up the traces */
      dcsum = sum(lc);

      /* Change to floating point */
      csum = FLOAT(dcsum);
      csum = Pade_b[l] * csum;
      pade_terms[l] = csum * minusone;
      trlns_unimp[i] += pade_terms[l];
      j++;
      
    }
  }

  /* Work out subtraction terms. First the estimators of the 
     traceless operators */

  /* 1 Work out plaquette for 4-th order term */
  MesPlq (u, w_plaq, s_plaq, t_plaq, link);
  
  f_plaq = FLOAT(w_plaq);

  
     
    bilinear = 0;

  
  
  /* Copy chi -> r */
  for(cb=0; cb < Ncb; cb++) {
    r[cb] = chi[cb];
  }
  
  
  /* Flip bc-s back on */
  phfctr (u, FORWARD);

  for(i=0; i < Nimpr; i++) { 
    
    /* Apply Dslash: Dchi = Dslash r */
    for(cb=0; cb < Ncb; cb++) {
      dslash (u, r[cb], Dchi[1-cb], PLUS, cb);
    }

    /* Now trace Dchi with chi to get the bilinear(term) chi D chi */
    lc = trace(adj[chi[0]] * Dchi[0]);
    for(cb = 1; cb < Ncb; cb++) {
      lct = trace(adj[chi[cb]] * Dchi[cb]);
      lc += lct;
    }
    
    /* Sum up the traces */
    dcsum = sum(lc);

    /* Change to floating point */
    bilinear[i] = FLOAT(dcsum);

    /* If i corresponds to an even powered term (recall i goes from 0) */
    if ( (i + 1) % 2 == 0 ) {

      /* Deal with the even power appropriately */
      switch( i + 1 ) { 
      case 2:
	/* Second order term is traceless */
	break;
      case 4:
	
	/* Here we will put the plaquette */
	c_plaq = 1;
	tmp = TO_REAL(64)*f_plaq*TO_REAL(3)*TO_REAL(6)*TO_REAL(Ncb)*TO_REAL(vol_cb);
	c_plaq = tmp * c_plaq;
	bilinear[i] += c_plaq;
	break;
      default:

	/* Anything higher than 4th order we will ignore for now */
	bilinear[i] = 0;
	break;
      }

    }

    /* Now copy Dchi back to R for the next iteration */
    for(cb = 0; cb < Ncb; cb++) { 
      r[cb] = Dchi[cb];
    }
  }

  phfctr (u, BACKWARD);

  /* We are done with Dchi -> NOw dealing with the bilinears only */
  
  for(kappa=0; kappa < Nkappa; kappa++) { 
    for(i = 0; i < Nimpr; i++) { 
      
      /* psib_offterm = kappa^{i+1} bilinear(i) */
      kappa_pow = pow((double)Kappas[kappa], (double)TO_REAL(i+1));
      
      psib_offterms[i][kappa] = kappa_pow * bilinear[i];

      for(l = 0; l < Npade; l++) { 
	
	/* Pade_offterm(order) = kappa^{i+1}/(1 + c(order))^{i+2} */
	tmp = TO_REAL(1) / (TO_REAL(1) + Pade_c[l]);
	pade_c_pow = pow((double)tmp, (double)TO_REAL(i+2));
	tmp = -Pade_b[l]*kappa_pow * pade_c_pow;
	trln_offterms[i][l][kappa] = tmp * bilinear[i];
      }
    }
  }

  
  sum_trln_offterms = 0;

  for(kappa=0; kappa < Nkappa; kappa++) { 
    for(i = 0; i < Nimpr; i++) { 
      for(l = 0; l < Npade; l++) { 
	sum_trln_offterms[i][kappa] += trln_offterms[i][l][kappa];
      }

      /* Multiply each term by the subtraction coefficient omega */
      sum_trln_offterms[i][kappa] = Omega_impr[i] * sum_trln_offterms[i][kappa];

      psib_offterms[i][kappa] = Omega_impr[i] * psib_offterms[i][kappa];
    }
  }

  /* push(xml_out,"PsiBarPsiOffterms");
write(xml_out, "psib_offterms", psib_offterms);
pop(xml_out);
     push(xml_out,"TrlnOffterms");
write(xml_out, "sum_trln_offterms", sum_trln_offterms);
pop(xml_out); */

  /* Now construct improved estimates */
  psibp_imp = psibp_unimp;
  trlns_imp = trlns_unimp;

  for(kappa=0; kappa < Nkappa; kappa++) { 

    for(i=0; i < Nimpr; i++) { 
      psibp_imp[kappa] -= psib_offterms[i][kappa];
      trlns_imp[kappa] -= sum_trln_offterms[i][kappa];
    }
  }

  /* Normalise psi bar psi */
  TrI = TO_REAL(1)/ TrI;

  for(kappa=0; kappa < Nkappa; kappa++) { 
    psibp_imp[kappa] = psibp_imp[kappa] * TrI;
    psibp_unimp[kappa] = psibp_unimp[kappa] * TrI;
  }


        
              DestDslash (D, WILSON_DSLASH);
       
      
  END_CODE();
}
