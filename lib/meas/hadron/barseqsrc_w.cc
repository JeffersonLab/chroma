// $Id: barseqsrc_w.cc,v 1.8 2005-01-14 18:42:35 edwards Exp $
/*! \file
 *  \brief Construct baryon sequential sources.
 */

#include "chromabase.h"
#include "util/ft/sftmom.h"
#include "meas/hadron/barseqsrc_w.h"

namespace Chroma {

//! Construct baryon sequential sources
/*!
 * \ingroup hadron
 *
 * This routine is specific to Wilson fermions!
 *
 * Construct baryon sequential sources.
 *
 * Note: only equal baryons at source and sink, for now.
 *
 * \param quark_propagator_1   first (u) quark propagator ( Read )
 * \param quark_propagator_2   second (d) quark propagator ( Read )
 * \param seq_src_prop         sequential source as propagator ( Write )
 * \param t_sink               time coordinate of the sink ( Read )
 * \param sink_mom             sink baryon momentum ( Read )
 * \param j_decay              direction of the exponential decay ( Read )
 * \param seq_src              flag indicating the type of the sequential source ( Read ) 
 */

void barSeqSource(const LatticePropagator& quark_propagator_1, 
		  const LatticePropagator& quark_propagator_2,
		  LatticePropagator& seq_src_prop, 
		  int t_sink, const multi1d<int>& sink_mom, 
		  int j_decay, int seq_src)
{
  START_CODE();

  LatticePropagator src_prop_tmp;
  LatticePropagator q1_tmp;
  LatticePropagator q2_tmp;
  LatticePropagator di_quark;
  LatticeColorMatrix col_mat;
  
  if ( Ns != 4 || Nc != 3 )		/* Code is specific to Ns=4 and Nc=3. */
  {
    END_CODE();
    return;
  }

  SpinMatrix g_one = 1.0;

  /* C = Gamma(10) */
  SpinMatrix Cgm = 0.5 * (Gamma(10) * (Gamma(2) * g_one  +  timesI(Gamma(1) * g_one)));

  /* C g_5 NR = (1/2)*C gamma_5 * ( 1 + g_4 ) */ 
  SpinMatrix Cg5NR = 0.5 * (Gamma(5) * (g_one + (g_one * Gamma(8))));


  switch (seq_src)
  {
  case 0:
    /* "\bar u O u" insertion in proton, ie. "(u C gamma_5 d) u" */
    /* T = (1 + gamma_4) / 2 = (1 + Gamma(8)) / 2 */
    {
      /* C gamma_5 = Gamma(5) = - (C gamma_5)^T */
      di_quark = quarkContract24(quark_propagator_1 * Gamma(5), 
				 Gamma(5) * quark_propagator_2);
      q1_tmp = 1;
      src_prop_tmp = (di_quark + Gamma(8)*di_quark) 
	 + traceSpin(di_quark)*(q1_tmp + Gamma(8)*q1_tmp);

      q1_tmp = Gamma(5) * quark_propagator_2 * Gamma(5);
      q2_tmp = quark_propagator_1 + quark_propagator_1*Gamma(8);
      src_prop_tmp -= quarkContract13(q1_tmp, q2_tmp) + quarkContract12(q2_tmp, q1_tmp);
      src_prop_tmp *= 0.5;
    }
    break;

  case 1:
    /* "\bar d O d" insertion in proton, ie. "(u C gamma_5 d) u" */
    /* T = (1 + gamma_4) / 2 = (1 + Gamma(8)) / 2 */
    {
      /* C gamma_5 = Gamma(5) = - (C gamma_5)^T */
      q2_tmp = quark_propagator_1 * Gamma(5);
      q1_tmp = q2_tmp + Gamma(8) * q2_tmp;
      q2_tmp = Gamma(5) * quark_propagator_1;
      src_prop_tmp = -quarkContract14(q1_tmp, q2_tmp);

      q1_tmp = q2_tmp * Gamma(5);
      q2_tmp = quark_propagator_1 + quark_propagator_1 * Gamma(8);
      src_prop_tmp -= quarkContract12(q2_tmp, q1_tmp);

      src_prop_tmp *= 0.5;
    }
    break;

  case 2:
    /* "\bar u O u" insertion in proton, ie. "(u C gamma_5 d) u" */
    /* T = \Sigma_3 (1 + gamma_4) / 2 = -i (Gamma(3) + Gamma(11)) / 2 */
    {
      /* C gamma_5 = Gamma(5) = - (C gamma_5)^T */
      q1_tmp = quark_propagator_1 * Gamma(5);
      q2_tmp = Gamma(5) * quark_propagator_2;
      di_quark = quarkContract24(q1_tmp, q2_tmp);
      src_prop_tmp = Gamma(3) * di_quark;
      src_prop_tmp += Gamma(11) * di_quark;

      col_mat = traceSpin(di_quark);
      q1_tmp = 1;
      di_quark = Gamma(3) * q1_tmp;
      di_quark += Gamma(11) * q1_tmp;
      src_prop_tmp += col_mat * di_quark;

      q1_tmp = q2_tmp * Gamma(5);
      q2_tmp = quark_propagator_1 * Gamma(3);
      q2_tmp += quark_propagator_1 * Gamma(11);
      src_prop_tmp -= quarkContract13(q1_tmp, q2_tmp) + quarkContract12(q2_tmp, q1_tmp);

      q1_tmp = 0.5 * timesMinusI(src_prop_tmp);
      src_prop_tmp = q1_tmp;
    }
    break;

  case 3:
    /* "\bar d O d" insertion in proton, ie. "(u C gamma_5 d) u" */
    /* T = \Sigma_3 (1 + gamma_4) / 2 = -i (Gamma(3) + Gamma(11)) / 2 */
    {
      /* C gamma_5 = Gamma(5) = - (C gamma_5)^T */
      q2_tmp = quark_propagator_1 * Gamma(5);
      q1_tmp = Gamma(3) * q2_tmp + Gamma(11) * q2_tmp;
      q2_tmp = Gamma(5) * quark_propagator_1;
      src_prop_tmp = -quarkContract14(q1_tmp, q2_tmp);

      q1_tmp = q2_tmp * Gamma(5);
      q2_tmp = quark_propagator_1 * Gamma(3) + quark_propagator_1 * Gamma(11);
      src_prop_tmp -= quarkContract12(q2_tmp, q1_tmp);

      q1_tmp = 0.5 * timesMinusI(src_prop_tmp);
      src_prop_tmp = q1_tmp;
    }
    break;

  case 4:
    /* "\bar u O u" insertion in Delta^+,
       ie. "2*(u C gamma_- d) u + (u C gamma_- u) d" */
    /* T = (1 + gamma_4) / 2 = (1 + Gamma(8)) / 2 */
    {
      /* C gamma_- = Cgm = (C gamma_-)^T */
      q1_tmp = quark_propagator_1 * Cgm;
      q2_tmp = Cgm * quark_propagator_2;
      di_quark = quarkContract24(q1_tmp, q2_tmp);
      src_prop_tmp = di_quark + Gamma(8) * di_quark;

      col_mat = traceSpin(di_quark);
      q1_tmp = 1;
      di_quark = q1_tmp + Gamma(8) * q1_tmp;
      src_prop_tmp += col_mat * di_quark;

      q1_tmp = q2_tmp * Cgm;
      q2_tmp = quark_propagator_1 + quark_propagator_1 * Gamma(8);
      src_prop_tmp += quarkContract13(q1_tmp, q2_tmp) + quarkContract12(q2_tmp, q1_tmp);

      q1_tmp = Cgm * quark_propagator_1;
      q2_tmp = q1_tmp * Cgm;
      q1_tmp = quark_propagator_2 + Gamma(8) * quark_propagator_2;
      src_prop_tmp += quarkContract12(q1_tmp, q2_tmp) + quarkContract24(q1_tmp, q2_tmp);

      q2_tmp = q1_tmp * Cgm;
      q1_tmp = Cgm * quark_propagator_1;
      src_prop_tmp += quarkContract14(q2_tmp, q1_tmp);

      q1_tmp = Cgm * quark_propagator_2;
      q2_tmp = q1_tmp + q1_tmp * Gamma(8);
      q1_tmp = quark_propagator_1 * Cgm;
      src_prop_tmp += quarkContract14(q1_tmp, q2_tmp) + quarkContract13(q1_tmp, q2_tmp);
      src_prop_tmp *= 2;
    }
    break;

  case 5:
    /* "\bar d O d" insertion in Delta^+,
       ie. "2*(u C gamma_- d) u + (u C gamma_- u) d" */
    /* T = (1 + gamma_4) / 2 = (1 + Gamma(8)) / 2 */
    {
      /* C gamma_- = Cgm = (C gamma_-)^T */
      q2_tmp = quark_propagator_1 * Cgm;
      q1_tmp = q2_tmp + Gamma(8) * q2_tmp;
      q2_tmp = Cgm * quark_propagator_1;
      src_prop_tmp = quarkContract14(q1_tmp, q2_tmp);

      q1_tmp = q2_tmp * Cgm;
      q2_tmp = quark_propagator_1 + quark_propagator_1 * Gamma(8);
      src_prop_tmp += quarkContract12(q2_tmp, q1_tmp);

      q1_tmp = quark_propagator_1 * Cgm;
      q2_tmp = Cgm * quark_propagator_1;
      di_quark = quarkContract24(q1_tmp, q2_tmp);
      src_prop_tmp += di_quark + Gamma(8) * di_quark;

      di_quark = quarkContract13(q1_tmp, q2_tmp);
      src_prop_tmp += di_quark + di_quark * Gamma(8);
      src_prop_tmp *= 2;

      col_mat = traceSpin(di_quark);
      q1_tmp = 1;
      q2_tmp = q1_tmp + Gamma(8) * q1_tmp;
      src_prop_tmp += col_mat * q2_tmp;
    }
    break;

  case 6:
    /* "\bar u O u" insertion in NR proton, ie. 
     * "(u C gamma_5 (1/2)(1 + gamma_4)  d) u" */
    /* T = (1 + gamma_4) / 2 = (1 + Gamma(8)) / 2 */
    {
      /* Use precomputed Cg5NR = C gamma_5 (1/2) ( 1 + g_4 ) */
      q1_tmp = quark_propagator_1 * Cg5NR;
      q2_tmp = Cg5NR * quark_propagator_2;
      di_quark = quarkContract24(q1_tmp, q2_tmp);
      src_prop_tmp = di_quark + Gamma(8) * di_quark;

      col_mat = traceSpin(di_quark);
      q1_tmp = 1;
      di_quark = q1_tmp + Gamma(8) * q1_tmp;
      src_prop_tmp += col_mat * di_quark;

      q1_tmp = q2_tmp * Cg5NR;
      q2_tmp = quark_propagator_1 + quark_propagator_1 * Gamma(8);
      src_prop_tmp -= quarkContract13(q1_tmp, q2_tmp) + quarkContract12(q2_tmp, q1_tmp);
      src_prop_tmp *= 0.5;
    }
    break;

  case 7:
    /* "\bar d O d" insertion in NR proton, ie. 
     * "(u C gamma_5 (1/2)(1 + gamma_4) d) u" */
    /* T = (1 + gamma_4) / 2 = (1 + Gamma(8)) / 2 */
    {
      /* C gamma_5 (1/2)( 1 + g_4)  = Cg5NR  */
      q2_tmp = quark_propagator_1 * Cg5NR;
      q1_tmp = q2_tmp + Gamma(8) * q2_tmp;
      q2_tmp = Cg5NR * quark_propagator_1;
      src_prop_tmp = -quarkContract14(q1_tmp, q2_tmp);

      q1_tmp = q2_tmp * Cg5NR;
      q2_tmp = quark_propagator_1 + quark_propagator_1 * Gamma(8);
      src_prop_tmp -= quarkContract12(q2_tmp, q1_tmp);
      src_prop_tmp *= 0.5;
    }
    break;

  case 8:
    /* "\bar u O u" insertion in NR proton, ie. 
     * "(u C gamma_5 (1/2)(1+gamma_4) d) u" */
    /* T = \Sigma_3 (1 + gamma_4) / 2 = -i (Gamma(3) + Gamma(11)) / 2 */
    {
      /* C gamma_5 (1/2) (1 + gamma_4) = Cg5NR */
      q1_tmp = quark_propagator_1 * Cg5NR;
      q2_tmp = Cg5NR * quark_propagator_2;
      di_quark = quarkContract24(q1_tmp, q2_tmp);
      src_prop_tmp = Gamma(3) * di_quark + Gamma(11) * di_quark;

      col_mat = traceSpin(di_quark);
      q1_tmp = 1;
      di_quark = Gamma(3) * q1_tmp  +  Gamma(11) * q1_tmp;
      src_prop_tmp += col_mat * di_quark;

      q1_tmp = q2_tmp * Cg5NR;
      q2_tmp = quark_propagator_1 * Gamma(3)  +  quark_propagator_1 * Gamma(11);
      src_prop_tmp -= quarkContract13(q1_tmp, q2_tmp) + quarkContract12(q2_tmp, q1_tmp);

      q1_tmp = 0.5 * timesMinusI(src_prop_tmp);
      src_prop_tmp = q1_tmp;
    }
    break;

  case 9:
    /* "\bar d O d" insertion in proton, ie. "(u C gamma_5 (1/2)(1+gamma_4)d) u" */
    /* T = \Sigma_3 (1 + gamma_4) / 2 = -i (Gamma(3) + Gamma(11)) / 2 */
    {
      /* C gamma_5 (1/2)(1+gamma_4)= Cg5NR */
      q2_tmp = quark_propagator_1 * Cg5NR;
      q1_tmp = Gamma(3) * q2_tmp  +  Gamma(11) * q2_tmp;
      q2_tmp = Cg5NR * quark_propagator_1;
      src_prop_tmp = -quarkContract14(q1_tmp, q2_tmp);

      q1_tmp = q2_tmp * Cg5NR;
      q2_tmp = quark_propagator_1 * Gamma(3)  +  quark_propagator_1 * Gamma(11);
      src_prop_tmp -= quarkContract12(q2_tmp, q1_tmp);

      q1_tmp = 0.5 * timesMinusI(src_prop_tmp);
      src_prop_tmp = q1_tmp;
    }
    break;

  case 21:
    /* "\bar u O u" insertion in NR proton, ie. 
     * "(u C gamma_5 (1/2)(1 + gamma_4)  d) u" */
    /* T = (1 + \Sigma_3)*(1 + gamma_4) / 2 
       = (1 + Gamma(8) - i G(3) - i G(11)) / 2 */

    /*
     *  Note that this is precisely src 6 + src 8, and corresponds
     *  to the combination of propagators used by the MIT group
     */

    {
      /* Use precomputed Cg5NR = C gamma_5 (1/2) ( 1 + g_4 ) */
      q1_tmp = quark_propagator_1 * Cg5NR;
      q2_tmp = Cg5NR * quark_propagator_2;
      di_quark = quarkContract24(q1_tmp, q2_tmp);

      /*
       *  We begin with the polarized piece
       */

      q1_tmp = Gamma(3) * di_quark  +  Gamma(11) * di_quark;
      /*
       *  Now multiply by -i
       */

      src_prop_tmp = timesMinusI(q1_tmp);

      /*
       * Now add the unpolarized piece
       */

      src_prop_tmp += di_quark  +  Gamma(8) * di_quark;

      /*
       *  Now the second term
       */

      col_mat = traceSpin(di_quark);
      q1_tmp = 1;

      /*
       *  First the polarized part
       */

      di_quark = Gamma(3) * q1_tmp  +  Gamma(11) * q1_tmp;
      q1_tmp = timesMinusI(di_quark); /*  Multiply by -i */
      src_prop_tmp += col_mat * q1_tmp;

      /*
       *  Now the unpolarized part
       */

      q1_tmp = 1;
      di_quark = q1_tmp  +  Gamma(8) * q1_tmp;
      src_prop_tmp += col_mat * di_quark;

      /*
       *  The third term...
       */

      q1_tmp = q2_tmp * Cg5NR;

      /*  First the polarized piece */

      di_quark = quark_propagator_1 * Gamma(3)  +  quark_propagator_1 * Gamma(11);
      q2_tmp = timesMinusI(di_quark);	/* Muliply by -i */

      /* Add the unpolarized piece */

      q2_tmp += quark_propagator_1  +  quark_propagator_1 * Gamma(8);

      src_prop_tmp -= quarkContract13(q1_tmp, q2_tmp) + quarkContract12(q2_tmp, q1_tmp);

      /*  Finally, multiply everything by 1/2 */

      src_prop_tmp *= 0.5;

    }
    break;

  case 22:
    /* "\bar d O d" insertion in proton, ie. "(u C gamma_5 (1/2)(1+gamma_4)d) u" */
    /* T = \Sigma_3 (1 + gamma_4) / 2 = -i (Gamma(3) + Gamma(11)) / 2 */
    {
      /* C gamma_5 (1/2)(1+gamma_4)= Cg5NR */
      q2_tmp = quark_propagator_1 * Cg5NR;

      /*
       *  First the polarized piece
       */

      di_quark = Gamma(3) * q2_tmp   +  Gamma(11) * q2_tmp;
      q1_tmp = timesMinusI(di_quark);

      /*
       *  Now add the unpolarized piece
       */

      q1_tmp += q2_tmp  +  Gamma(8) * q2_tmp;


      q2_tmp = Cg5NR * quark_propagator_1;
      src_prop_tmp = -quarkContract14(q1_tmp, q2_tmp);

      /* Now the second term */

      q1_tmp = q2_tmp * Cg5NR;

      /*
       *  First the polarized piece
       */

      di_quark = quark_propagator_1 * Gamma(3) + quark_propagator_1 * Gamma(11);
      q2_tmp = timesMinusI(di_quark);

      /*
       *  Now add the unpolarized piece
       */

      q2_tmp += quark_propagator_1 + quark_propagator_1 * Gamma(8);

      src_prop_tmp -= quarkContract12(q2_tmp, q1_tmp);
      src_prop_tmp *= 0.5;

    }
    break;

  default:
    QDP_error_exit("Unknown sequential source type", seq_src);
  }

  /* Now take hermitian conjugate and multiply on both sides
     with gamma_5 = Gamma(15) */
  {
    /* This bit does the hermitian conjugation. */
    
    q2_tmp = adj(src_prop_tmp);
    src_prop_tmp = q2_tmp;
  }
   
  /* src_prop_tmp now holds src^{dagger} */
  /* If we are working with half inversions I need to do an 
     extra projection here */

  /* Now slap the gamma_5 on either end */
  q1_tmp = src_prop_tmp * Gamma(15);
  src_prop_tmp = Gamma(15) * q1_tmp;
        
  /*
   *  We now inject momentum at sink if required
   */
  multi1d<LatticeInteger> my_coord(Nd);

  bool nonzero = false;
  for(int mu=0, j=0; mu < Nd; mu++)
  {
    my_coord[mu] = Layout::latticeCoordinate(mu);	/* Obtains the muth coordinate */

    if (mu != j_decay)
    {
      if(sink_mom[j] != 0)
	nonzero = true;

      j++;
    }
  }

  // multiply in the phase if required
  if (nonzero)
  {
    LatticeReal p_dot_x = 0;
    for(int mu=0, j=0; mu < Nd; ++mu)
    {
      if (mu == j_decay)
        continue;

      p_dot_x += my_coord[mu] * sink_mom[j] * twopi / Real(Layout::lattSize()[mu]);
      j++;
    }
            
    src_prop_tmp *= cmplx(cos(p_dot_x),sin(p_dot_x));
  }

  /*
   * Now mask out all but sink time slice
   */
  seq_src_prop = where(my_coord[j_decay] == t_sink,
		       src_prop_tmp,
		       LatticePropagator(zero));
        
  END_CODE();
}

}  // end namespace Chroma
