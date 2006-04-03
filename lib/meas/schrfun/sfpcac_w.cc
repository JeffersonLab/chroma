// $Id: sfpcac_w.cc,v 3.0 2006-04-03 04:59:04 edwards Exp $

#error "NOT FULLY CONVERTED"

/* This routine is specific to Wilson fermions! */

/* NOTE: this routine assumes the chiral basis for the gamma matrices, */
/* in particular the specific forms of gamma_0 (or gamma_4, which here is */
/* actually Gamma(8)) and of gamma_5! */

/* Compute correlation functions between axial current or pseudescalar */
/* density and boundary fields using Schroedinger BC. */

/* Also computed, on demand, are correlation functions between both */
/* boundaries with zero, one (vector current) and two (axial current or */
/* pseudoscalar density) insertions in the bulk. */

/* Compute quark propagators by inverting the Wilson operator using */
/* the kappa's in the array "Kappa".  */
/* The initial guess for the inverter is zero. */

/* The results are written to the namelist file. */

/* For further details see the comments in the dependent subroutines. */

/* u -- gauge field ( Read ) */
/* j_decay  -- direction along which the exponentrial is allowed to decay ( Read ) */
/* Kappa    -- array of kappa's to use in the inversion ( Read ) */
/* ClCoeff  -- array of clover coefficients to use in the inversion ( Read ) */
/* numKappa -- length of the Kappa array. ( Read ) */
/* RsdCG    -- CG accuracy ( Read ) */
/* ncg_had  -- total number of CG iterations for quark proagators ( Modify ) */
/* ZVfactP -- flag for doing Z_V measurements ( Read ) */
/* ZAfactP -- flag for doing Z_A measurements ( Read ) */
/* x0, y0 -- time slices with axial current insertions ( Read ) */

include(types.mh)

SUBROUTINE(SFpcac, u, j_decay, Kappa, ClCoeffR, ClCoeffT, numKappa, 
	   RsdCG, ncg_had, ZVfactP, ZAfactP, x0, y0)

multi1d<LatticeColorMatrix> u(Nd);
int j_decay;
multi1d<Real> Kappa(numKappa);
multi1d<Real> ClCoeffR(numKappa);
multi1d<Real> ClCoeffT(numKappa);
Real RsdCG;
int numKappa;
int ncg_had;
int ZVfactP;
int ZAfactP;
int x0;
int y0;

{ /* Local Variables */
  include(COMMON_DECLARATIONS)
  LatticeFermion chi;
  LatticeFermion psi;
  LatticePropagator quark_propagator;
  LatticePropagator quark_prop_f;
  LatticePropagator quark_prop_b;
  LatticePropagator tmp_prop;
  
  LatticeFermion tmp1;
  LatticeFermion tmp2;
  LatticeFermion tmp3;
  LatticeReal r_tmp1;
  LatticeReal r_tmp2;
  LatticeReal faa_tmp;
  LatticeReal fap_tmp;
  LatticeReal fpp_tmp;
  LatticeInteger t_coord;
  LatticeBoolean tmask;

  multi1d<Double> hrsum(length);
  Double ddummy;
  multi1d<Real> pseudo_prop(length);
  multi1d<Real> axial_prop(length);
  multi1d<Real> pseudo_prop_b(length);
  multi1d<Real> axial_prop_b(length);
  multi1d<Real> vector_corr(length);

  multi1d<DPropagator> slice_prop(length);
  Propagator kprop;
  
  Real Kappa_mes;
  Real ClovCoeff_mes;
  Real ClovCoeffR_mes;
  Real ClovCoeffT_mes;
  Real ftmp;
  Real norm;
  Real f_1;
  Real f_AA;
  Real f_AP_PA;
  Real f_PP;
  Real fd_AA;
  Real fd_AP_PA;
  Real fd_PP;
  Real dummy;
  multi1d<int> t_source(Nd);
  int length;
  int colour_source;
  int spin_source;
  int spin_source2;
  int G5;
  int jd;
  int n;
  int cb;
  int loop;
  int direction;
  int tmin;
  int tmax;
  int t;
  int t_eff;
  
  START_CODE();
  
  if ( SchrFun == NO ) 		/* Code is only for Schroedinger BC */
  {
    END_CODE();
    return;
  }
  
  if ( Ns != 4 )
  {
    END_CODE();
    return;
  }
  
  /* Because of implied structure of (1 +/- gamma_0) in the chiral basis: */
  if ( j_decay != 3 )
    QDP_error_exit("SFpcac requires j_decay=3", j_decay);

  if ( numKappa < 1 )
    QDP_error_exit("illegal value", numKappa);
  
  if ( ZAfactP == YES && x0 < y0 )
    QDP_error_exit("Z_A computation requires x0 > y0", x0, y0);

  /* Multiply in the fermion phases */
  phfctr (u, FORWARD);

  /* The dimension of the propagators */
  length = nrow[j_decay];
  
  G5 = Ns*Ns-1;
  jd = INTEGER_LSHIFT_FUNCTION(1,j_decay);

        
  if (ZAfactP == YES)
  {
          }

  if (ZVfactP == YES || ZAfactP == YES)
  {
              }

  norm = 2;
  for(n = 0; n < Nd; ++n)
    if (n != j_decay)
      norm *= nrow[n];

  norm = TO_REAL(1) / norm;


  /* Location of upper wall source */
  if (GlueImp == 0)
  {
    tmin = 0;
    tmax = cb_nrow[j_decay] - 2;
  }
  else
  {
    tmin = 1;
    tmax = cb_nrow[j_decay] - 3;
  }

  for(loop = 0; loop < numKappa; ++loop)
  {
    if (ZVfactP == YES || ZAfactP == YES)
    {
      quark_prop_f = 0;
      quark_prop_b = 0;
    }

    Kappa_mes = Kappa[loop];
    if (AnisoP == YES)
    {
      ClovCoeffR_mes = ClovCoeffR = ClCoeffR[loop];
      ClovCoeffT_mes = ClovCoeffT = ClCoeffT[loop];
    }
    else
    {
      ClovCoeff_mes = ClovCoeff = ClCoeffR[loop];
    }

    
    for(direction = -1; direction <= 1; direction+=2)
    {
      t_source = 0;
      if (direction == -1)
	t_source[j_decay] = tmax;
      else
	t_source[j_decay] = tmin;


      pseudo_prop = 0;
      axial_prop = 0;
      quark_propagator = 0;

      for(colour_source = 0; colour_source < Nc; ++colour_source)
      {
	for(spin_source = 0; spin_source < Ns/2; ++spin_source)
	{
	  spin_source2 = spin_source + 2;

	  /* Compute quark propagator "psi" using source "chi" with type specified */
	  /* by WALL_SOURCE, and colour and spin equal to colour_source and spin_source. */
	  	  	  	  	  maksrc (tmp1, colour_source, spin_source, t_source, j_decay, OPTION[WALL_SOURCE]);

	  if (direction == -1)
	    for(cb = 0; cb < 2; ++cb)
	    {
	      tmp2 = u[j_decay][cb] * tmp1[cb];
	      psi[cb] = tmp2;
	      psi[cb] -= Gamma[jd] * tmp2;
	      chi[cb] = psi[cb] * Kappa_mes;
	    }
	  else
	    for(cb = 0; cb < 2; ++cb)
	    {
	      tmp2 = adj(u[j_decay][cb]) * tmp1[cb];
	      psi[cb] = tmp2;
	      psi[cb] += Gamma[jd] * tmp2;
	      chi[1-cb] = shift(psi[cb], 1-cb, BACKWARD, j_decay) * Kappa_mes;
	    }

	  	  
	  psi = 0;
	  Qprop (u, chi, Kappa_mes, RsdCG, psi, ncg_had);

	  /* Store in quark_prpagator and in quark_prop_f/b */
	  if (direction == -1)
	  {
	    chi = -psi;
	    for(cb = 0; cb  < 2; ++cb )
	    {
	      /* top spin contribution */
	      transf (psi[cb], quark_propagator[cb], colour_source, spin_source, cb, FORWARD);
	      /* bottom spin contribution */
	      transf (chi[cb], quark_propagator[cb], colour_source, spin_source2, cb, FORWARD);
	    }

	    /* Store the first two source spin components, also in slots */
	    /* of the last LAST two! (Recall that in the chiral basis used */
	    /* the last two components are minus the first two, so we
	       actually store q_prop * gamma_5!) */
	    if (ZVfactP == YES || ZAfactP == YES)
	      for(cb = 0; cb  < 2; ++cb )
	      {
		transf (psi[cb], quark_prop_b[cb], colour_source, spin_source, cb, FORWARD);
		transf (psi[cb], quark_prop_b[cb], colour_source, spin_source2, cb, FORWARD);
	      }
	  }
	  else
	  {
	    for(cb = 0; cb  < 2; ++cb )
	    {
	      /* top spin contribution */
	      transf (psi[cb], quark_propagator[cb], colour_source, spin_source, cb, FORWARD);
	      /* bottom spin contribution */
	      transf (psi[cb], quark_propagator[cb], colour_source, spin_source2, cb, FORWARD);
	    }

	    if (ZVfactP == YES || ZAfactP == YES)
	    {
	      	      	      for(cb = 0; cb  < 2; ++cb )
	      {
		tmp2 = adj(u[j_decay][cb]) * psi[cb];
		tmp3 = Gamma[jd] * tmp2;
		tmp2 += tmp3;
		tmp3 = tmp2 * Kappa_mes;

		transf (tmp3, quark_prop_f[cb], colour_source, spin_source, cb, FORWARD);
		transf (tmp3, quark_prop_f[cb], colour_source, spin_source2, cb, FORWARD);
	      }
	      	      	    }
	  }

	  	  	}
      }

      /* Construct the pion axial current divergence and the pion correlator */
      SFcurcor (u, quark_propagator, pseudo_prop, axial_prop, j_decay);

      /* Time reverse backwards source */
      if (direction == -1)
      {
	if ( ZAfactP == YES )
	  for(t = 0; t < length; t++)
	  {
	    pseudo_prop_b[t] = pseudo_prop[t] * norm;
	    axial_prop_b[t] = axial_prop[t] * norm;
	  }

	for(t = 0; t < (length-1)/2 + 1; ++t)
	{
	  t_eff = length - t - 1;

	  ftmp = pseudo_prop[t];
	  pseudo_prop[t] = pseudo_prop[t_eff];
	  pseudo_prop[t_eff] = ftmp;

	  ftmp = axial_prop[t];
	  axial_prop[t] = -axial_prop[t_eff];
	  axial_prop[t_eff] = -ftmp;
	}
      }

      /* Normalize to compare to everybody else */
      for(t = 0; t < length; ++t)
      {
	pseudo_prop[t] = pseudo_prop[t] * norm;
	axial_prop[t] = axial_prop[t] * norm;
      }

      if (AnisoP == YES)
      {
	push(xml_out,"PCAC_measurements");
write(xml_out, "loop", loop);
write(xml_out, "direction", direction);
write(xml_out, "AnisoP", AnisoP);
write(xml_out, "Kappa_mes", Kappa_mes);
write(xml_out, "ClovCoeffR_mes", ClovCoeffR_mes);
write(xml_out, "ClovCoeffT_mes", ClovCoeffT_mes);
write(xml_out, "pseudo_prop", pseudo_prop);
write(xml_out, "axial_prop", axial_prop);
pop(xml_out);
      }
      else
      {
	push(xml_out,"PCAC_measurements");
write(xml_out, "loop", loop);
write(xml_out, "direction", direction);
write(xml_out, "AnisoP", AnisoP);
write(xml_out, "Kappa_mes", Kappa_mes);
write(xml_out, "ClovCoeff_mes", ClovCoeff_mes);
write(xml_out, "pseudo_prop", pseudo_prop);
write(xml_out, "axial_prop", axial_prop);
pop(xml_out);
      }
    }

    
    if (ZVfactP == YES || ZAfactP == YES)
    {
      /* Sum the forward propagator over time slices to get Kprop */
      ftmp = TO_REAL(2) * norm;
                  tmp_prop = quark_prop_f[0];
      tmp_prop += quark_prop_f[1];
      slice_prop = sumMulti(tmp_prop, timeslice);
      kprop = FLOAT(slice_prop[tmax]);
      kprop = kprop * ftmp;
      
      /* quark_prop_f is no longer needed, and can be re-used below */

      /* Construct f_1 */
      f_1 = real(trace(adj[kprop] * kprop));
      ftmp = WORD_VALUE(WORD_ftmp,HALF) / Kappa_mes;
      f_1 *= WORD_VALUE(WORD_ftmp,HALF)*ftmp*ftmp;

      /* Construct H'' = H' gamma_5 K, where H' is the propagator */
      /* from the upper boundary. The gamma_5 multiplication was done. */
      /* H' is no longer needed: quark_prop_b can be overwritten. */
      for(cb = 0; cb < 2; cb++)
      {
	tmp_prop = quark_prop_b[cb] * kprop;
	quark_prop_b[cb] = tmp_prop;
      }

          }

    if (ZVfactP == YES)
    {
      /* Construct f_V */
                              
      n = G5 ^ jd;
      r_tmp1 = 0;
      for(cb = 0; cb < 2; cb++)
      {
        tmp_prop = Gamma[n] * quark_propagator[cb];
        r_tmp2 = real(trace(adj[quark_prop_b[cb]] * tmp_prop));
	r_tmp1 += r_tmp2;
      }

      ftmp = norm / (TO_REAL(2)*Kappa_mes);
      hrsum = sumMulti(r_tmp1, timeslice);
      for(t = 0; t < length; t++)
      {
	dummy = FLOAT(hrsum[t])
	vector_corr[t] = dummy * ftmp;
      }

      if (AnisoP == YES)
      {
	push(xml_out,"ZV_measurements");
write(xml_out, "loop", loop);
write(xml_out, "Kappa_mes", Kappa_mes);
write(xml_out, "ClovCoeffR_mes", ClovCoeffR_mes);
write(xml_out, "ClovCoeffT_mes", ClovCoeffT_mes);
write(xml_out, "f_1", f_1);
write(xml_out, "vector_corr", vector_corr);
pop(xml_out);
      }
      else
      {
	push(xml_out,"ZV_measurements");
write(xml_out, "loop", loop);
write(xml_out, "Kappa_mes", Kappa_mes);
write(xml_out, "ClovCoeff_mes", ClovCoeff_mes);
write(xml_out, "f_1", f_1);
write(xml_out, "vector_corr", vector_corr);
pop(xml_out);
      }
	
                                  }

    if (ZAfactP == YES)
    {
      /* Construct the ingredients for f^I_AA */
                                                                  faa_tmp = 0;
      fap_tmp = 0;
      fpp_tmp = 0;
      t_coord = Layout::latticeCoordinate(j_decay);

      /* "right" A_0 insertion at x */
      tmask = t_coord == x0;
      quark_prop_f = 0;
      for(colour_source = 0; colour_source < Nc; ++colour_source)
	for(spin_source = 0; spin_source < Ns/2; ++spin_source)
	{
	  spin_source2 = spin_source + 2;

	  for(cb = 0; cb < 2; ++cb )
	    transf (chi[cb], quark_propagator[cb], colour_source, spin_source, cb, BACKWARD);

	  n = jd ^ G5;
	  for(cb = 0; cb  < 2; ++cb )
	    psi[cb] = Gamma[n] * chi[cb];
	  /* This gives multiplication with gamma_5 * gamma_0. */
	  /* Include the minus sign for multiplication with */
	  /* gamma_0 * gamma_5 below. */

	  chi = 0;
	  for(cb = 0; cb < 2; ++cb )
	    copymask(chi[cb], tmask, psi[cb], NEGATE);

	  psi = 0;
	  Qprop (u, chi, Kappa_mes, RsdCG, psi, ncg_had);

	  for(cb = 0; cb < 2; ++cb )
	  {
	    transf (psi[cb], quark_prop_f[cb], colour_source, spin_source, cb, FORWARD);
	    transf (psi[cb], quark_prop_f[cb], colour_source, spin_source2, cb, FORWARD);
	  }
	}

      /* "left" P insertion at y */
      r_tmp1 = 0;
      for(cb = 0; cb < 2; ++cb )
      {
	r_tmp2 = real(trace(adj[quark_prop_b[cb]] * quark_prop_f[cb]));
	r_tmp1 -= r_tmp2;
      }
      t = y0 + 1;
      tmask = t_coord == t;
      copymask(fap_tmp, tmask, r_tmp1, ADD);
      t = y0 - 1;
      tmask = t_coord == t;
      copymask(fap_tmp, tmask, r_tmp1, SUBTRACT);

      /* "left" A_0 insertion at y */
      r_tmp1 = 0;
      for(cb = 0; cb  < 2; ++cb )
      {
	tmp_prop = Gamma[jd] * quark_prop_f[cb];
	r_tmp2 = real(trace(adj[quark_prop_b[cb]] * tmp_prop));
	r_tmp1 += r_tmp2;
      }
      tmask = t_coord == y0;
      copymask(faa_tmp, tmask, r_tmp1, ADD);

      /* "right" A_0 insertion at y */
      tmask = t_coord == y0;
      n = jd ^ G5;
      quark_prop_f = 0;
      for(colour_source = 0; colour_source < Nc; ++colour_source)
	for(spin_source = 0; spin_source < Ns/2; ++spin_source)
	{
	  spin_source2 = spin_source + 2;

	  for(cb = 0; cb < 2; ++cb )
	    transf (chi[cb], quark_propagator[cb], colour_source, spin_source, cb, BACKWARD);

	  for(cb = 0; cb  < 2; ++cb )
	    psi[cb] = Gamma[n] * chi[cb];
	  /* This gives multiplication with gamma_5 * gamma_0. */
	  /* Include the minus sign for multiplication with */
	  /* gamma_0 * gamma_5 below. */

	  chi = 0;
	  for(cb = 0; cb < 2; ++cb )
	    copymask(chi[cb], tmask, psi[cb], NEGATE);

	  psi = 0;
	  Qprop (u, chi, Kappa_mes, RsdCG, psi, ncg_had);

	  for(cb = 0; cb < 2; ++cb )
	  {
	    transf (psi[cb], quark_prop_f[cb], colour_source, spin_source, cb, FORWARD);
	    transf (psi[cb], quark_prop_f[cb], colour_source, spin_source2, cb, FORWARD);
	  }
	}

      /* "left" P insertion at x */
      r_tmp1 = 0;
      for(cb = 0; cb < 2; ++cb )
      {
	r_tmp2 = real(trace(adj[quark_prop_b[cb]] * quark_prop_f[cb]));
	r_tmp1 += r_tmp2;
      }
      t = x0 + 1;
      tmask = t_coord == t;
      copymask(fap_tmp, tmask, r_tmp1, ADD);
      t = x0 - 1;
      tmask = t_coord == t;
      copymask(fap_tmp, tmask, r_tmp1, SUBTRACT);

      /* "left" A_0 insertion at x */
      r_tmp1 = 0;
      for(cb = 0; cb  < 2; ++cb )
      {
	tmp_prop = Gamma[jd] * quark_prop_f[cb];
	r_tmp2 = real(trace(adj[quark_prop_b[cb]] * tmp_prop));
	r_tmp1 -= r_tmp2;
      }
      tmask = t_coord == x0;
      copymask(faa_tmp, tmask, r_tmp1, ADD);

      /* "right" P insertion at x */
      quark_prop_f = 0;
      for(colour_source = 0; colour_source < Nc; ++colour_source)
	for(spin_source = 0; spin_source < Ns/2; ++spin_source)
	{
	  spin_source2 = spin_source + 2;

	  for(cb = 0; cb < 2; ++cb )
	    transf (chi[cb], quark_propagator[cb], colour_source, spin_source, cb, BACKWARD);

	  for(cb = 0; cb  < 2; ++cb )
	    psi[cb] = Gamma[G5] * chi[cb];

	  chi = 0;
	  t = x0 + 1;
	  tmask = t_coord == t;
	  for(cb = 0; cb < 2; ++cb )
	    copymask(chi[cb], tmask, psi[cb], REPLACE);
	  t = x0 - 1;
	  tmask = t_coord == t;
	  for(cb = 0; cb < 2; ++cb )
	    copymask(chi[cb], tmask, psi[cb], SUBTRACT);

	  psi = 0;
	  Qprop (u, chi, Kappa_mes, RsdCG, psi, ncg_had);

	  for(cb = 0; cb < 2; ++cb )
	  {
	    transf (psi[cb], quark_prop_f[cb], colour_source, spin_source, cb, FORWARD);
	    transf (psi[cb], quark_prop_f[cb], colour_source, spin_source2, cb, FORWARD);
	  }
	}

      /* "left" P insertion at y */
      r_tmp1 = 0;
      for(cb = 0; cb < 2; ++cb )
      {
	r_tmp2 = real(trace(adj[quark_prop_b[cb]] * quark_prop_f[cb]));
	r_tmp1 -= r_tmp2;
      }
      t = y0 + 1;
      tmask = t_coord == t;
      copymask(fpp_tmp, tmask, r_tmp1, ADD);
      t = y0 - 1;
      tmask = t_coord == t;
      copymask(fpp_tmp, tmask, r_tmp1, SUBTRACT);

      /* "left" A_0 insertion at y */
      r_tmp1 = 0;
      for(cb = 0; cb  < 2; ++cb )
      {
	tmp_prop = Gamma[jd] * quark_prop_f[cb];
	r_tmp2 = real(trace(adj[quark_prop_b[cb]] * tmp_prop));
	r_tmp1 += r_tmp2;
      }
      tmask = t_coord == y0;
      copymask(fap_tmp, tmask, r_tmp1, ADD);

      /* "right" P insertion at y */
      quark_prop_f = 0;
      for(colour_source = 0; colour_source < Nc; ++colour_source)
	for(spin_source = 0; spin_source < Ns/2; ++spin_source)
	{
	  spin_source2 = spin_source + 2;

	  for(cb = 0; cb < 2; ++cb )
	    transf (chi[cb], quark_propagator[cb], colour_source, spin_source, cb, BACKWARD);

	  for(cb = 0; cb  < 2; ++cb )
	    psi[cb] = Gamma[G5] * chi[cb];

	  chi = 0;
	  t = y0 + 1;
	  tmask = t_coord == t;
	  for(cb = 0; cb < 2; ++cb )
	    copymask(chi[cb], tmask, psi[cb], REPLACE);
	  t = y0 - 1;
	  tmask = t_coord == t;
	  for(cb = 0; cb < 2; ++cb )
	    copymask(chi[cb], tmask, psi[cb], SUBTRACT);

	  psi = 0;
	  Qprop (u, chi, Kappa_mes, RsdCG, psi, ncg_had);

	  for(cb = 0; cb < 2; ++cb )
	  {
	    transf (psi[cb], quark_prop_f[cb], colour_source, spin_source, cb, FORWARD);
	    transf (psi[cb], quark_prop_f[cb], colour_source, spin_source2, cb, FORWARD);
	  }
	}

      /* "left" P insertion at x */
      r_tmp1 = 0;
      for(cb = 0; cb < 2; ++cb )
      {
	r_tmp2 = real(trace(adj[quark_prop_b[cb]] * quark_prop_f[cb]));
	r_tmp1 += r_tmp2;
      }
      t = x0 + 1;
      tmask = t_coord == t;
      copymask(fpp_tmp, tmask, r_tmp1, ADD);
      t = x0 - 1;
      tmask = t_coord == t;
      copymask(fpp_tmp, tmask, r_tmp1, SUBTRACT);

      /* "left" A_0 insertion at x */
      r_tmp1 = 0;
      for(cb = 0; cb  < 2; ++cb )
      {
	tmp_prop = Gamma[jd] * quark_prop_f[cb];
	r_tmp2 = real(trace(adj[quark_prop_b[cb]] * tmp_prop));
	r_tmp1 -= r_tmp2;
      }
      tmask = t_coord == x0;
      copymask(fap_tmp, tmask, r_tmp1, ADD);

                                          
            ftmp = norm / TO_REAL(2);
      hrsum = sumMulti(faa_tmp, timeslice);
      dummy = FLOAT(hrsum[x0]);
      f_AA = dummy * norm;
      dummy = FLOAT(hrsum[y0]);
      f_AA += dummy * norm;

      hrsum = sumMulti(fap_tmp, timeslice);
      dummy = FLOAT(hrsum[x0]);
      f_AP_PA = dummy * ftmp;
      dummy = FLOAT(hrsum[y0]);
      f_AP_PA += dummy * ftmp;
      t = x0 + 1;
      dummy = FLOAT(hrsum[t]);
      f_AP_PA += dummy * ftmp;
      t = x0 - 1;
      dummy = FLOAT(hrsum[t]);
      f_AP_PA += dummy * ftmp;
      t = y0 + 1;
      dummy = FLOAT(hrsum[t]);
      f_AP_PA += dummy * ftmp;
      t = y0 - 1;
      dummy = FLOAT(hrsum[t]);
      f_AP_PA += dummy * ftmp;

      ftmp = norm / TO_REAL(4);
      hrsum = sumMulti(fpp_tmp, timeslice);
      t = x0 + 1;
      dummy = FLOAT(hrsum[t]);
      f_PP = dummy * ftmp;
      t = x0 - 1;
      dummy = FLOAT(hrsum[t]);
      f_PP += dummy * ftmp;
      t = y0 + 1;
      dummy = FLOAT(hrsum[t]);
      f_PP += dummy * ftmp;
      t = y0 - 1;
      dummy = FLOAT(hrsum[t]);
      f_PP += dummy * ftmp;

      fd_AA = axial_prop_b[y0] * axial_prop[x0];
      fd_AA -= axial_prop_b[x0] * axial_prop[y0];

      ftmp = 0.5;
      fd_AP_PA = pseudo_prop_b[y0+1] * axial_prop[x0];
      fd_AP_PA -= pseudo_prop_b[y0-1] * axial_prop[x0];
      fd_AP_PA += axial_prop_b[y0] * pseudo_prop[x0+1];
      fd_AP_PA -= axial_prop_b[y0] * pseudo_prop[x0-1];
      fd_AP_PA -= axial_prop_b[x0] * pseudo_prop[y0+1];
      fd_AP_PA += axial_prop_b[x0] * pseudo_prop[y0-1];
      fd_AP_PA -= pseudo_prop_b[x0+1] * axial_prop[y0];
      fd_AP_PA += pseudo_prop_b[x0-1] * axial_prop[y0];
      fd_AP_PA = fd_AP_PA * ftmp;

      ftmp *= ftmp;
      fd_PP = pseudo_prop_b[y0+1] * pseudo_prop[x0+1];
      fd_PP -= pseudo_prop_b[y0+1] * pseudo_prop[x0-1];
      fd_PP -= pseudo_prop_b[y0-1] * pseudo_prop[x0+1];
      fd_PP += pseudo_prop_b[y0-1] * pseudo_prop[x0-1];
      fd_PP -= pseudo_prop_b[x0+1] * pseudo_prop[y0+1];
      fd_PP += pseudo_prop_b[x0+1] * pseudo_prop[y0-1];
      fd_PP += pseudo_prop_b[x0-1] * pseudo_prop[y0+1];
      fd_PP -= pseudo_prop_b[x0-1] * pseudo_prop[y0-1];
      fd_PP = fd_PP * ftmp;

      if (AnisoP == YES)
      {
	push(xml_out,"ZA_measurements");
write(xml_out, "loop", loop);
write(xml_out, "Kappa_mes", Kappa_mes);
write(xml_out, "ClovCoeffR_mes", ClovCoeffR_mes);
write(xml_out, "ClovCoeffT_mes", ClovCoeffT_mes);
write(xml_out, "f_1", f_1);
write(xml_out, "f_AA", f_AA);
write(xml_out, "f_AP_PA", f_AP_PA);
write(xml_out, "f_PP", f_PP);
write(xml_out, "fd_AA", fd_AA);
write(xml_out, "fd_AP_PA", fd_AP_PA);
write(xml_out, "fd_PP", fd_PP);
pop(xml_out);
      }
      else
      {
	push(xml_out,"ZA_measurements");
write(xml_out, "loop", loop);
write(xml_out, "Kappa_mes", Kappa_mes);
write(xml_out, "ClovCoeff_mes", ClovCoeff_mes);
write(xml_out, "f_1", f_1);
write(xml_out, "f_AA", f_AA);
write(xml_out, "f_AP_PA", f_AP_PA);
write(xml_out, "f_PP", f_PP);
write(xml_out, "fd_AA", fd_AA);
write(xml_out, "fd_AP_PA", fd_AP_PA);
write(xml_out, "fd_PP", fd_PP);
pop(xml_out);
      }

                            }
  }
  
  if (ZVfactP == YES || ZAfactP == YES)
  {
              }
  
  if (ZAfactP == YES)
  {
          }
  
      
  /* Turn off the fermion phases */
  phfctr (u, BACKWARD);

  END_CODE();
}
