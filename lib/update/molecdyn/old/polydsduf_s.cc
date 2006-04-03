/* $Id: polydsduf_s.cc,v 3.0 2006-04-03 04:59:10 edwards Exp $ ($Date: 2006-04-03 04:59:10 $)								   */
/* polydsduf								   */

/* This routine is specific to staggered fermions!                         */

/* subroutine for pseudofermionic force calculation                        */
/* in the polynomial hybrid Monte Carlo and related algorithms.            */
/* (S. Sint, 9/97)  						           */

/* The pseudofermionic action is of the form  S_f = chi^dag P_n(Q) chi     */
/* where Q is the rescaled staggered kernel: 				   */

/*     Q   = mdagm_resc * M^dagM                                           */     
/*     M^dagM      = 1 - Kappa^2 D'  D'                                    */
/*     mdagm_resc  = PolyArgResc/(1+4*Nd^2*Kappa^2)                        */

/* and P_n(x) is a polynomial of degree n which approximates the function  */
/* x^(-r/s) in the interval [eps,1] with r = PolyPowNum and s = PolyPowDen.*/

/* The routine has been checked against the standard HMC force routine     */
/* dsduf for r=s=1  (S. Sint, 12/97)                                       */ 

/* polydsduf -- computes the derivative of the pseudofermionic action with */ 
/* respect to the link field                                               */


/* u -- gauge field           ( Read )   */

/*         |  dS      dS_f               */
/* ds_u -- | ----   + -----   ( Modify ) */
/*         |  dU       dU                */

/* psi -- pseudofermion       ( Read )   */

void polydsduf(multi1d<LatticeColorMatrix>& ds_u,
	       const multi1d<LatticeColorMatrix>& u,
	       const LatticeFermion psi)
{
  multi1d<LatticeFermion> chi(PolyDegHalf+1);
  multi1d<LatticeFermion> v_mu(Nd);
  multi1d<LatticeFermion> w_mu(Nd);
  LatticeFermion v1;
  LatticeFermion v2;
  LatticeFermion phi;

  Real const1;
  Real zero_real;
  Real mdagm_fact;
  Double rho_norm;

  Complex root;
  Complex dummy;

  int mu;
  int cb;
  int Ncb;
  int j;
  int jj;
  int nh;

  START_CODE("subroutine");;

/* All fields live on 1 checkerboard only: */

  Ncb=1;

                

  mdagm_resc=PolyArgResc/(TO_REAL(1)+TO_REAL(4*Nd*Nd)*KappaMD*KappaMD);
  mdagm_fact=PolyNorm*mdagm_resc;

  const1=KappaMD*KappaMD*mdagm_fact;

  zero_real = 0;

/* compute the PolyDegHalf vectors chi_j,   j= 0,..,PolyDegHalf    */
/* Note: C indices start at zero and the index of PolyRoots is     */
/* thus shifted by one as compared to the notes.                   */


  chi[0] = psi;

  for (j=1; j <= PolyDegHalf; j++)
  {  
    jj=j-1;
    root = PolyRoots[jj];
    mdagm (u, chi[j-1], KappaMD, chi[j], Ncb);
    chi[j] = mdagm_fact * chi[j];
    v1 = chi[j-1];
    chi[j] -= root * v1;

/*--------------------------------------------------------*/
/* the following is not accepted by the macro compiler    */
/*							  */
/*      chi[j] -= root * chi[j-1];          */
/*							  */
/* because it cannot tell that j-1 and j are different    */
/* and in place multiplication of a complex number times  */
/* a lattice fermion (complex components) is              */
/* of course not possible.                                */
/*--------------------------------------------------------*/

  }


/* start force computation */

  nh=PolyDegHalf;
  phi = chi[nh];

  for (j=1; j <= PolyDegHalf; j++)
  {
    for(mu = 0;mu  < ( Nd); ++mu )
    {
      cb = 1;
      v1(rb[cb]) = shift(chi[nh-j], FORWARD, mu);
      v_mu[mu] = u[mu][cb] * v1;            
    }

/* add v_mu(mu) over mu and subtract backward part of Dslash to obtain */
/* v1=Dslash*chi(nh-j)        					       */

    v1 = v_mu[0];
    for (mu = 1; mu  < ( Nd); ++mu )
    {
      v1 += v_mu[mu];
    }
    for( mu = 0; mu  < ( Nd); ++mu )
    {
      cb = 0;
      v1 -= shift(adj[u[mu][cb]], (1-cb), BACKWARD, mu) * shift(chi[nh-j], (1-cb), BACKWARD, mu);
    }
 

/* compute w_mu */

    for ( mu = 0;mu  < ( Nd); ++mu )
    {
      cb = 0;
      v2(rb[cb]) = shift(v1, FORWARD, mu);
      w_mu[mu] = u[mu][cb] * v2;            
    }

       
/* contribution of the first term (out of 4) to the force            */
/* (current w_mu not needed thereafter, can be multiplied by const1) */

    cb = 0;
    for ( mu = 0;mu  < ( Nd); ++mu )
    {
      w_mu[mu] = w_mu[mu] * const1;
      ds_u[mu][cb] += w_mu[mu] * adj(phi);
    }

/* compute new w_mu */

    for ( mu = 0;mu  < ( Nd); ++mu )
    {
      cb = 1;
      v2(rb[cb]) = shift(phi, FORWARD, mu);
      w_mu[mu] = u[mu][cb] * v2;            
    }      

/* contribution of the second term (out of 4) to the force        */
/*(current v1 not needed thereafter, can be multiplied by const1) */

    cb = 1;
    v1 = v1 * const1;
    for ( mu = 0;mu  < ( Nd); ++mu )
    {
      ds_u[mu][cb] -= w_mu[mu] * adj(v1);
    }


/* add w_mu(mu) over mu and subtract backward part of Dslash to obtain  */
/* v1=Dslash*phi        						*/

    v1 = w_mu[0];
    for (mu = 1; mu  < ( Nd); ++mu )
    {
      v1 += w_mu[mu];
    }
    for( mu = 0; mu  < ( Nd); ++mu )
    {
      cb = 0;
      v1 -= shift(adj[u[mu][cb]], (1-cb), BACKWARD, mu) * shift(phi, (1-cb), BACKWARD, mu);
    }

/* contribution of the third term (out of 4) to the force            */
/* (current v_mu not needed thereafter, can be multiplied by const1) */

    cb = 1;
    for(mu = 0; mu < Nd; ++mu)
    {
      v_mu[mu] = v_mu[mu] * const1;
      ds_u[mu][cb] -= v_mu[mu] * adj(v1);
    }

/* compute new v_mu */

    for(mu = 0; mu < Nd; ++mu)
    {
      cb = 0;
      v2(rb[cb]) = shift(v1, FORWARD, mu);
      v_mu[mu] = u[mu][cb] * v2;            
    }      

/* contribution of the fourth term (out of 4) to the force                 */
/* (current chi(nh-j) not needed thereafter, can be multiplied by const1) */

    cb = 0;
    chi[nh-j] = chi[nh-j] * const1;
    for(mu = 0; mu < Nd; ++mu)
    {
      ds_u[mu][cb] += v_mu[mu] * adj(chi[nh-j]);
    }

    if (j< PolyDegHalf)
    {

/* construct new phi = chi_{n/2+j}, using v_mu */


      v2 = v_mu[0];
      for(mu = 1; mu < Nd; ++mu)
      {
	v2 += v_mu[mu];
      }
      for(mu = 0; mu < Nd; ++mu)
      {
	cb = 1;
	v2 -= shift(adj[u[mu][cb]], (1-cb), BACKWARD, mu) * shift(v1, (1-cb), BACKWARD, mu);
      }
 


/* dummy (complex) = mdagm_fact (real) minus root (complex)  */

      root = adj(PolyRoots[nh-j]);
      dummy = 0;
      dummy = -root
	dummy += cmplx(mdagm_fact,zero_real);

      v1 = phi * dummy;
      phi = v1;
      phi -= v2 * const1;

    }

/* end of loop over j */

  }  

                
  END_CODE("subroutine");;
}
