// $Id: stagpbp_s.cc,v 3.0 2006-04-03 04:59:04 edwards Exp $
/*! \file
 *  \brief Calculates noise estimator for the staggered trace 
 */

#error "Converted but not tested"

#include "chromabase.h"
#include "meas/pbp/stagpbp_s.h"

namespace Chroma {

//! Calculates noise estimator for the staggered trace 

/*! Calculates noise estimator for the staggered trace used in the
 *           following constructions
 *   <psibar_psi(x)> = Tr D^{-1} 
 *   <psibar psi(x) psibar psi(y)> = Tr D^{-2}
 *   <psibar "gamma_5" psi psibar "gamma_5" psi> = Tr (D^dag D)^{-1}
 *                                 = Tr H^{-2}

 * What is computed are the traces so iso(flavor) triplets and singlets can be
 * reconstructed 

 * This routine is specific to staggered fermions!

 * u           -- gauge field ( Read )
 * nhit        -- number of Gaussian hits ( Read )
 * Mass        -- array of quark mass values ( Read )
 * numMass     -- number of quark masses ( Read )
 * RsdCG       -- the CG accuracy used for multishift CG ( Read ) 
 */

void StagPbp(XMLWriter& xml_out,
	     const StaggeredFermActBase& S_f,
	     Handle<const ConnectState> state,
	     const multi1d<Real>& Mass, 
	     enum InvType invType,
	     const multi1d<Real>& RsdCG,
	     int MaxCG, 
	     int nhit)
{
  START_CODE();

  const int numMass = Mass.size();

  int n_congrd;
  int n_congrd_tot;
  multi1d<Double> TrDinv(numMass);
  multi1d<Double> TrHinv_sq(numMass);
  multi1d<Double> TrDinv_sq(numMass);
  Double TrBdag_B;
  multi1d<Double> TrBdag_Eta(numMass);
  multi1d<Double> TrEtadag_Eta(numMass);
  Double TrBdag_B_avg;
  multi1d<Double> TrBdag_Eta_avg(numMass);
  multi1d<Double> TrEtadag_Eta_avg(numMass);
  multi1d<LatticeStaggeredFermion> eta(numMass);
  LatticeStaggeredFermion b;
  LatticeReal lrtrace_aux;

  Double ddummy1;
  Double two;
  Double mass;
  int ihit;
  int i;

  /* Code is specific to staggered fermions */
#if 0
  if (FermAct != STAGGERED)
    QDP_error_exit("only support staggered fermions", FermAct);
#endif

        
  TrBdag_B_avg = 0;
  TrBdag_Eta_avg = 0;
  TrEtadag_Eta_avg = 0;
  n_congrd_tot = 0;

    
  for(ihit = 1; ihit <= nhit; ++ihit)
  {
    /* Fill with random gaussian noise such that < b_dagger * b > = 1 per d.o.f. */
    gaussian(b);
    n_congrd = 0;

#if 0
    /* If using Schroedinger functional, zero out the boundaries */
    S_f.getFermBC().modifyF(b);
#endif

    /* Keep track of actual norm of noise */
    TrBdag_B = norm2(b);
       

    /* Prop. solution */
    /* eta = (D*D_dag)^(-1) * b   */
    MStQprop (u, b, Mass, numMass, RsdCG, eta, 4, n_congrd);


    /* Construct, using that b and eta live only on the even sublattice   */
    /*   Tr Bdag_B     = < b_dag * b >    */
    /*   Tr Bdag_Eta   = < b_dag * eta >     */
    /*   Tr Etadag_Eta = < eta_dag * eta >     */
    /* Use these to make */
    /*   Tr D^{-1} = m * < b_dag * eta > */
    /*   Tr (D^dag D)^{-1} = Tr H^{-2} = < b_dag * eta > */
    /*   Omega = Tr (D^dag D)^{-1}  + Tr D^{-2} */
    /*         = 2 * m^2 * < eta^dag * eta > */
    /*   Tr D^{-2} = Omega - Tr (D^dag D)^{-1} */
    
    TrBdag_Eta = 0;
    TrEtadag_Eta = 0;

    for(i=0; i < numMass; ++i)
    {
      lrtrace_aux = real(trace(adj[b] * eta[i]));
      TrBdag_Eta[i] += sum(lrtrace_aux);

      TrEtadag_Eta[i] += norm2(eta[i]);
    }

    

    /* Add on to running average */
    TrBdag_B_avg += TrBdag_B;
    TrBdag_Eta_avg += TrBdag_Eta;
    TrEtadag_Eta_avg += TrEtadag_Eta;
    n_congrd_tot += n_congrd;

    /* Normalize and print */
    ddummy1 = TO_DOUBLE(1) / TO_DOUBLE(vol_cb);
    TrBdag_B = TrBdag_B * ddummy1;
    for(i=0; i < numMass; ++i)
    {
      TrBdag_Eta[i] = TrBdag_Eta[i] * ddummy1;
      TrEtadag_Eta[i] = TrEtadag_Eta[i] * ddummy1;
    }

    push(xml_out,"Staggered_trace_hit");
    write(xml_out, "ihit", ihit);
    write(xml_out, "nhit", nhit);
    write(xml_out, "TrBdag_B", TrBdag_B);
    write(xml_out, "TrBdag_Eta", TrBdag_Eta);
    write(xml_out, "TrEtadag_Eta", TrEtadag_Eta);
    write(xml_out, "n_congrd", n_congrd);
    pop(xml_out);
  }

    
  /* Normalize */
  ddummy1 = TO_DOUBLE(1) / TO_DOUBLE(nhit * vol_cb);
  TrBdag_B_avg = TrBdag_B_avg * ddummy1;
  for(i=0; i < numMass; ++i)
  {
    TrBdag_Eta_avg[i] = TrBdag_Eta_avg[i] * ddummy1;
    TrEtadag_Eta_avg[i] = TrEtadag_Eta_avg[i] * ddummy1;
  }

  /* For niceties and redundancies, turn into the usual traces */
  /* The traces are linear combinations of the inner products, so */
  /* the average goes through */
  /* Use these to make */
  /*   Tr D^{-1} = m * < b_dag * eta > */
  /*   Tr (D^dag D)^{-1} = Tr H^{-2} = < b_dag * eta > */
  /*   Omega = Tr (D^dag D)^{-1}  + Tr D^{-2} */
  /*         = 2 * m^2 * < eta^dag * eta > */
  /*   Tr D^{-2} = Omega - Tr (D^dag D)^{-1} */
  /*   Tr D^{-2} = Omega - Tr H^{-2} */
      
  two = 2;

  for(i=0; i < numMass; ++i)
  {
    mass = FLOAT(Mass[i]);
    ddummy1 = two * mass * mass;
    TrHinv_sq[i] = TrBdag_Eta_avg[i];
    TrDinv[i] = TrHinv_sq[i] * mass;
    TrDinv_sq[i] = -TrHinv_sq[i];
    TrDinv_sq[i] += TrEtadag_Eta_avg[i] * ddummy1;
  }
    
  push(xml_out,"Staggered_trace_avg");
  write(xml_out, "nhit", nhit);
  write(xml_out, "TrBdag_B_avg", TrBdag_B_avg);
  write(xml_out, "TrBdag_Eta_avg", TrBdag_Eta_avg);
  write(xml_out, "TrEtadag_Eta_avg", TrEtadag_Eta_avg);
  write(xml_out, "TrDinv", TrDinv);
  write(xml_out, "TrHinv_sq", TrHinv_sq);
  write(xml_out, "TrDinv_sq", TrDinv_sq);
  write(xml_out, "n_congrd_tot", n_congrd_tot);
  pop(xml_out);

      
        
  END_CODE();
}

}  // end namespace Chroma
