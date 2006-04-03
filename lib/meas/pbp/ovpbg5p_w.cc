// $Id: ovpbg5p_w.cc,v 3.0 2006-04-03 04:59:04 edwards Exp $
/*! \file
 *  \brief Calculates noise estimator for the overlap trace 
 */

#error "Converted but not tested"

#include "chromabase.h"
#include "meas/pbp/ovpbg5p_w.h"


namespace Chroma {

//! OVPBG5P - Calculates noise estimator for the overlap trace
/*!
 * OVPBG5P - Calculates noise estimator for the overlap trace used in the
 *           following constructions
 *   <psibar_psi(x)> = Tr D^{-1}
 *   <psibar psi(x) psibar psi(y)> = Tr D^{-2} - Tr D{-1}  
 *   <psibar gamma_5 psi psibar gamma_5 psi> = = Tr H^{-2}    

 * What is computed are the traces so iso(flavor) triplets and singlets can be
 * reconstructed  

 * This routine is specific to Wilson and really Overlap fermions! 

 * state       -- gauge field ( Read ) 
 * nhit        -- number of Gaussian hits ( Read ) 
 * n_zero      -- the topology of the field ( Read ) 
 * Mass        -- array of quark mass values ( Read ) 
 * numMass     -- number of quark masses ( Read ) 
 * RsdCG       -- the CG accuracy used for multishift CG ( Read ) 
 */


void OvPbg5p(XMLWriter& xml_out,
	     const OverlapFermActBase& S_f,
	     Handle<const ConnectState> state,
	     const multi1d<Real>& Mass, 
	     enum InvType invType,
	     const multi1d<Real>& RsdCG,
	     int MaxCG, 
	     int n_zero, 
	     int nhit)
{
  START_CODE();

  const int numMass = Mass.size();

  multi1d<Double> TrDinv(numMass);
  multi1d<Double> TrHinv_sq(numMass);
  multi1d<Double> TrDinv_sq(numMass);
  multi1d<Double> TrBdag_Eta(numMass);
  multi1d<Double> TrEtadag_Eta(numMass);
  multi1d<Double> TrBdag_Eta_avg(numMass);
  multi1d<Double> TrEtadag_Eta_avg(numMass);
  multi2d<LatticeFermion> eta(Nsubl, numMass);
  LatticeFermion b;

  Double ddummy1;
  Double ddummy2;

  int G5 = Ns*Ns - 1;
        
  Double TrBdag_B;
  Double TrBdag_B_avg = 0;
  TrBdag_Eta_avg = 0;
  TrEtadag_Eta_avg = 0;
  int n_congrd = 0;
  int n_congrd_tot = 0;

    
  for(int ihit = 1; ihit <= nhit; ++ihit)
  {
    /* Fill with random gaussian noise such that < b_dagger * b > = 1 per d.o.f. */
    gaussian(b);

#if 0
    /* If using Schroedinger functional, zero out the boundaries */
    S_f.getFermBC().modifyF(b);
#endif

    /* Make the source chiral */
    if (n_zero > 0)
    {
      /* Have pos. chirality zero modes, so use neg. chirality for source */
      LatticeFermion tmp = chiralProjectMinus(b);
      b = tmp;
    }
    else
    {
      /* Have either neg. chirality or no zero modes, so use pos. chirality for source */
      LatticeFermion tmp = chiralProjectPlus(b);
      b = tmp;
    }
    
    /* Keep track of actual norm of noise */
    TrBdag_B = norm2(b);
       

    /* Prop. solution */
    /* eta = (D*D_dag)^(-1) * b   */
    /* The inverse of the un-preconditioned matrix is needed here. This */
    /* is given by MOvQprop! */
    n_congrd = 0;
    S_f.multiQprop(eta, Mass, state, b, invType, RsdCG, 4, MaxCG, n_congrd);


    /* Construct   */
    /*   Tr Bdag_B     = < b_dag * eta >    */
    /*   Tr Bdag_Eta   = < b_dag * eta >     */
    /*   Tr Etadag_Eta = < eta_dag * eta >     */
    /* Use these to make */
    /*   Tr D^{-1} = 2 * (mu/(1-mu^2)) * <b_dag * (eta - b)> */
    /*   Tr H^{-2} = 2 * (1/(1-mu^2)) * <b_dag * (eta - b)> */
    /*   Omega = Tr H^{-2}  + Tr D^{-2} */
    /*         = 4 * (mu^2/(1-mu^2)^2) * <(eta - b)^dag * (eta - b)> */
    /*   Tr D^{-2} = Omega - Tr H^{-2} */
    
    TrBdag_Eta = 0;
    TrEtadag_Eta = 0;

    for(int i=0; i < numMass; ++i)
    {
      TrBdag_Eta[i] += real(innerProductReal(b,eta[i]));
      TrEtadag_Eta[i] += norm2(eta[i]);
    }

    /* Add on to running average */
    TrBdag_B_avg += TrBdag_B;
    TrBdag_Eta_avg += TrBdag_Eta;
    TrEtadag_Eta_avg += TrEtadag_Eta;
    n_congrd_tot += n_congrd;

    /* Normalize and print */
    ddummy1 = Double(1) / Double(Layout::vol());
    TrBdag_B *= ddummy1;
    for(int i=0; i < numMass; ++i)
    {
      TrBdag_Eta[i] *= ddummy1;
      TrEtadag_Eta[i] *= ddummy1;
    }

    if (nhit > 1)
    {
      push(xml_out,"Overlap_trace_hit");
      write(xml_out, "ihit", ihit);
      write(xml_out, "nhit", nhit);
      write(xml_out, "TrBdag_B", TrBdag_B);
      write(xml_out, "TrBdag_Eta", TrBdag_Eta);
      write(xml_out, "TrEtadag_Eta", TrEtadag_Eta);
      write(xml_out, "n_congrd", n_congrd);
      pop(xml_out);
    }
  }

    
  /* Normalize */
  ddummy1 = Double(1) / Double(nhit * Layout::vol());
  TrBdag_B_avg *= TrBdag_B_avg;
  for(int i=0; i < numMass; ++i)
  {
    TrBdag_Eta_avg[i] *= ddummy1;
    TrEtadag_Eta_avg[i] *= ddummy1;
  }

  /* For niceties and redundancies, turn into the usual traces */
  /* The traces are linear combinations of the inner products, so */
  /* the average goes through */
  /* Use these to make */
  /*   Tr D^{-1} = 2 * (mu/(1-mu^2)) * <b_dag * (eta - b)> */
  /*   Tr H^{-2} = 2 * (1/(1-mu^2)) * <b_dag * (eta - b)> */
  /*   Omega = Tr H^{-2}  + Tr D^{-2} */
  /*         = 4 * (mu^2/(1-mu^2)^2) * <(eta - b)^dag * (eta - b)> */
  /*         = 4 * (mu^2/(1-mu^2)^2) * <b_dag*b - 2*b_dag*eta + eta_dag*eta> */
  /*   Tr D^{-2} = Omega - Tr H^{-2} */
      
  for(int i=0; i < numMass; ++i)
  {
    Real mass = Mass[i];
    ddummy1 = 2.0 / Double(1.0 - mass*mass);
    ddummy2 = TrBdag_Eta_avg[i] - TrBdag_B_avg;

    TrHinv_sq[i] = ddummy1 * ddummy2;
    TrDinv[i] = TrHinv_sq[i] * mass;

    ddummy1 *= mass;
    ddummy2 = TrEtadag_Eta_avg[i] - 2.0*TrBdag_Eta_avg[i] + TrBdag_B_avg;
    ddummy2 *= ddummy1*ddummy1;
    TrDinv_sq[i] = ddummy2 - TrHinv_sq[i];
  }
    
  push(xml_out,"Overlap_trace_avg");
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
