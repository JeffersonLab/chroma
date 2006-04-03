// $Id: sfcurcor_w.cc,v 3.0 2006-04-03 04:59:04 edwards Exp $

#error "NOT FULLY CONVERTED"

/* This routine is specific to Wilson fermions! */

/* Construct 'current correlators' and axial density used for the PCAC determination
 *  in the Schroedinger Functional
 */

/* u -- gauge field ( Read ) */
/* quark_propagator -- quark propagator ( Read ) */
/* pseudoscalar_propagator -- pion correlator ( Modify ) */
/* axial_current -- axial-current to pion_1 correlators ( Modify ) */
/* j_decay -- direction of the exponential decay ( Read ) */

/*         ____ */
/*         \ */
/* cc(t) =  >  < m(0, 0) c(t + L, x) > */
/*         /                     */
/*         ---- */
/*           x */
include(types.mh)

SUBROUTINE(SFcurcor, u, quark_propagator, pseudoscalar_propagator, 
	   axial_current, j_decay)

multi1d<LatticeColorMatrix> u(Nd);
LatticePropagator quark_propagator;
multi1d<Real> pseudoscalar_propagator(length);
multi1d<Real> axial_current(length);
int j_decay;

{ /* Local Variables. */
  include(COMMON_DECLARATIONS)
  LatticePropagator tmp_prop;
  LatticeReal psi_tmp;
  LatticeReal psi_sq;

  multi1d<Double> hsum(length);
  Real dummy;
  int length;
  int t;
  int t_eff;
  int jd;
  int cb;

  START_CODE();

  length = nrow[j_decay];

        
  /* First compute the pion correlator - the pseudoscalar */
  psi_sq = 0;
  for(cb = 0; cb < 2; ++cb)
  {
    psi_tmp = real(trace(adj[quark_propagator[cb]] * quark_propagator[cb]));
    psi_sq += psi_tmp;
  }

  /* Do a slice-wise sum. */
  hsum = sumMulti(psi_sq, timeslice);

  for(t = 0; t < length; ++t)
  {
    t_eff = INTEGER_MOD_FUNCTION(t + length, length);

    dummy = FLOAT(hsum[t]);
    pseudoscalar_propagator[t_eff] += dummy;
  }


  /* Construct the axial-current to pion correlator */
  jd = INTEGER_LSHIFT_FUNCTION(1,j_decay);
  psi_sq = 0;
  for(cb = 0; cb < 2; ++cb)
  {
    tmp_prop = Gamma[jd] * quark_propagator[cb];
    psi_tmp = real(trace(adj[quark_propagator[cb]] * tmp_prop));
    psi_sq -= psi_tmp;
  }

  /* Do a slice-wise sum. */
  hsum = sumMulti(psi_sq, timeslice);

  for(t = 0; t < length; ++t)
  {
    t_eff = INTEGER_MOD_FUNCTION(t + length, length);

    dummy = FLOAT(hsum[t]);
    axial_current[t_eff] += dummy;
  }

        
  END_CODE();
}
