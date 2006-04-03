namespace Chroma {


/* $Id: pbg5p_w.cc,v 3.0 2006-04-03 04:59:04 edwards Exp $ */

/* PBG5P - Calculates noise estimator for the trace used in the
/*         following constructions */

/*   <psibar_psi(x)> = Tr D^{-1}  */
/*   <psibar psi(x) psibar psi(y)> = Tr D^{-2} - Tr D{-1}  */
/*   <psibar gamma_5 psi psibar gamma_5 psi> = = Tr H^{-2}    */

/* What is computed are the traces so iso(flavor) triplets and singlets can be
/* reconstructed */ 

/* This routine is specific to Wilson and really Overlap fermions! */

/* u           -- gauge field ( Read ) */
/* nhit        -- number of Gaussian hits ( Read ) */
/* Mass        -- array of quark mass values ( Read ) */
/* numMass     -- number of quark masses ( Read ) */
/* RsdCG       -- the CG accuracy used for multishift CG ( Read ) */

include(types.mh)

void Pbg5p(XMLWriter& xml_out,
	   const WilsonTypeFermAct< multi1d<LatticeFermion> >& S_f,
	   Handle<const ConnectState> state,
	   const multi1d<Real>& Mass, 
	   enum InvType invType,
	   const multi1d<Real>& RsdCG,
	   int MaxCG,
	   int nhit)
{
  START_CODE();

  int n_congrd;
  int n_congrd_tot;
  multi1d<Double> TrDinv(numMass);
  multi1d<Double> TrHinv_sq(numMass);
  multi1d<Double> TrDinv_avg(numMass);
  multi1d<Double> TrHinv_sq_avg(numMass);
  LatticeFermion eta;
  LatticeFermion eta_sq;
  LatticeFermion tmp;
  LatticeFermion b;
  LatticeFermion b_orig;
  LatticeReal lrtrace_aux;

  Double ddummy1;
  Double ddummy2;
  Double one;
  Double two;
  Double mass;
  int i;

  int G5 = Ns*Ns - 1;
  isign = PLUS;

        
  TrDinv_avg = 0;
  TrHinv_sq_avg = 0;
  n_congrd = 0;
  n_congrd_tot = 0;

      
  for(int ihit = 1; ihit <= nhit; ++ihit)
  {
    /* Fill with random gaussian noise such that < b_dagger * b > = 1 per d.o.f. */
    gaussian(b);

#if 1
    /* If using Schroedinger functional, zero out the boundaries */
    S_f.getFermBC().modifyF(b);
#endif

    /* Make the source chiral */
    LatticeFermion tmp = chiralProjectPlus(b);
    b = tmp;
    
    /* Prop. solution */
    /* eta = D^(-1) * b   */
    /* eta_sq = (D^dag*D)^(-1) * b   */
    /* The inverse of the un-preconditioned matrix is needed here. This */
    /* is given by Qprop! */
    /* Construct   */
    /*   Tr Dinv     = < b_dag * eta >    */
    /*   Tr Hinv_sq  = < b_dag * eta_sq >     */
    
    TrDinv = 0;
    TrHinv_sq = 0;

    n_congrd = 0;
    eta = 0;
    eta_sq = 0;
    b_orig = b;

    for(i=0; i < numMass; ++i)
    {
      gaussian(tmp);
      CHIRAL_PROJECT(tmp, isign, eta, REPLACE);
      CHIRAL_PROJECT(tmp, isign, eta_sq, REPLACE);

      Qprop (u, b, Mass[i], RsdCG[i], eta, n_congrd);
      b = b_orig;

      tmp = Gamma(G5)*eta;
      Qprop (u, tmp, Mass[i], RsdCG[i], eta_sq, n_congrd);
      tmp = Gamma(G5)*eta_sq;
      eta_sq = tmp;
      
      TrDinv[i] += real(innerProduct(b, eta));
      TrHinv_sq[i] += real(innerProduct(b, eta_sq));
    }

        
    /* Normalize and print */
    ddummy1 = Double(1) / Double(vol);
    for(i=0; i < numMass; ++i)
    {
      TrDinv[i] *= ddummy1;
      TrHinv_sq[i] *= ddummy1;
    }

    if (nhit > 1)
    {
      push(xml_out,"Pbpg5_trace_hit");
      write(xml_out, "ihit", ihit);
      write(xml_out, "nhit", nhit);
      write(xml_out, "TrDinv", TrDinv);
      write(xml_out, "TrHinv_sq", TrHinv_sq);
      write(xml_out, "n_congrd", n_congrd);
      pop(xml_out);
    }

    /* Add on to running average */
    TrDinv_avg += TrDinv;
    TrHinv_sq_avg += TrHinv_sq;
    n_congrd_tot += n_congrd;
  }

      
  /* Normalize */
  ddummy1 = Double(1) / Double(nhit);
  for(i=0; i < numMass; ++i)
  {
    TrDinv_avg[i] *= ddummy1;
    TrHinv_sq_avg[i] *= ddummy1;
  }

  push(xml_out,"Pbp5g_trace_avg");
  write(xml_out, "nhit", nhit);
  write(xml_out, "Mass", Mass);
  write(xml_out, "TrDinv_avg", TrDinv_avg);
  write(xml_out, "TrHinv_sq_avg", TrHinv_sq_avg);
  write(xml_out, "n_congrd_tot", n_congrd_tot);
  pop(xml_out);

       
  END_CODE();
}

}  // end namespace Chroma
