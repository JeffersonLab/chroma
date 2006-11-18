/*  $Id: asq_dsl_s.cc,v 3.1 2006-11-18 02:33:03 kostas Exp $  */

#include "chromabase.h"
#include "actions/ferm/linop/asq_dsl_s.h"


namespace Chroma 
{ 
  //! Creation routine
  /*!
   *
   * This routine is specific to staggered fermions!
   *
   * Description:
   *
   * This routine applies the "asq" or "asqtad" operator D' to Psi, 
   * putting the result in Chi.
   *
   *	       Nd-1
   *	       ---
   *	       \                     F                     
   *   chi(x)  :=  >  isign eta  (x) [U  (x) psi(x+mu)
   *	       /             mu	     mu
   *	       ---
   *	       mu=0
   *
   *			+ c_3 U  (x) U  (x+mu) U  (x+2mu) psi(x+3mu) ]
   *                             mu     mu        mu
   *
   *	             Nd-1
   *	             ---
   *	             \                      +F
   *                -    >  isign eta  (x)  [U  (x-mu) psi(x-mu)
   *	             /             mu	    mu
   *	             ---
   *	             mu=0
   *
   *                             +      +          +
   *			+ c_3 U  (x) U  (x-2mu) U  (x-3mu) psi(x-3mu) ]
   *                             mu     mu         mu
   * Note the KS phase factors are already included in the U's!
   *
   * Arguments:
   *
   *
   *  U_fat     'fat link' Gauge field				(Read)
   *  U_triple    'triple-links' UUU				(Read)
   *  Psi	      Pseudofermion field				(Read)
   *  Chi	      Pseudofermion field				(Write)
   *		      +
   *  ISign      D' or D'  ( +1 | -1 ) respectively		(Read)
   *  CB	      Checkerboard of OUTPUT vector			(Read)
   *
   * NOTE: the coefficient c_3 is included in u_triple! 
   *
   *
   *  recoded by mcneile
   */
  void QDPStaggeredDslash::create(Handle<AsqtadConnectStateBase> state_)
  { 
    START_CODE();

    state = state_;

    //  XMLFileWriter xml_out("output2.xml");
    //  push(xml_out, "more_tests");
    //  write(xml_out, "u_fat", u_fat);
    //  write(xml_out,"u_triple", u_triple);
    //  pop(xml_out);

    END_CODE();
  }

  void QDPStaggeredDslash::apply (LatticeStaggeredFermion& chi, const LatticeStaggeredFermion& psi, enum PlusMinus isign, int cb) const
  {
    START_CODE();

    const multi1d<LatticeColorMatrix>& u_fat = state->getFatLinks();
    const multi1d<LatticeColorMatrix>& u_triple = state->getTripleLinks();

    // need convention on isign
    //
    // isign == PLUS is normal isign == MINUS is daggered
    //

    LatticeStaggeredFermion tmp_0 = zero;
    LatticeStaggeredFermion tmp_1 = zero;
    LatticeStaggeredFermion tmp_2 = zero;

    int mu;

    /* Forward one-hop and three-hop neigbhors */
    /* Note the KS phase factors are already included in the U's! */
    mu = 0;
    tmp_0 = shift(psi, FORWARD, mu);
    chi[rb[cb]] = u_fat[mu] * tmp_0; // Color matrix X color vector : 66 flops
    //  NEIGHBOUR(tmp_0, tmp_1, FORWARD, mu, 1);
    tmp_1 = shift(tmp_0, FORWARD, mu);
    tmp_2 = shift(tmp_1, FORWARD, mu);
    chi[rb[cb]] += u_triple[mu] * tmp_2; // +66 + 6  flops

    for(mu = 1;mu  <= ( Nd - 1); ++mu )
    {
      tmp_0 = shift(psi, FORWARD, mu);
      chi[rb[cb]] += u_fat[mu] * tmp_0; //+66 + 6 flops
      //    NEIGHBOUR(tmp_0, tmp_1, FORWARD, mu, 1);
      tmp_1 = shift(tmp_0, FORWARD, mu);
      tmp_2 = shift(tmp_1, FORWARD, mu);
      chi[rb[cb]] += u_triple[mu] * tmp_2; //+66 + 6 flops
    }

    // Up to now 570 flops

    /* Backward one-hop and three-hop neigbhors */
    /* Note the KS phase factors are already included in the U's! */
    for(mu = 0;mu  <= ( Nd - 1); ++mu )
    {
      chi[rb[cb]] -= shift(adj(u_fat[mu]), BACKWARD, mu) * shift(psi,BACKWARD, mu);
      // +66 + 6 flops
      tmp_0 = shift(adj(u_triple[mu]),  BACKWARD, mu) * shift(psi, BACKWARD,  mu);
      // +66 flops
      // NEIGHBOUR(tmp_0, tmp_1, BACKWARD, mu, 1);
      tmp_1 = shift(tmp_0, BACKWARD, mu);
      tmp_2 = shift(tmp_1, BACKWARD, mu);

      chi[rb[cb]] -= tmp_2; //+6 flops
    }


    // this requires more thought
    // No it doesnt
    // const int MINUS = -1 ;
     
    if(isign == MINUS)
      chi = -chi;

    END_CODE();
  }

} // End Namespace Chroma
