// $Id: mesonseqsrc_w.cc,v 1.1 2003-12-17 03:56:46 edwards Exp $
/*! \file
 *  \brief Construct meson sequential sources.
 */

#include "chromabase.h"
#include "util/ft/sftmom.h"
#include "meas/hadron/mesonseqsrc_w.h"

using namespace QDP;

//! Construct a meson sequential source.
/*!
 * \ingroup hadron
 *
 * This routine is specific to Wilson fermions!
 *
 *  delta(tz-tx) exp(i p.z) \gamma_5 G \gamma_5
 *
 * \param quark_propagator   input quark propagator ( Read ) 
 * \param seq_src_prop       sequential source as propagator ( Write ) 
 * \param t_sink             time coordinate of the sink ( Read ) 
 * \param sink_mom           sink pion momentum ( Read ) 
 * \param j_decay            direction of the exponential decay ( Read ) 
 * \param seq_src            the particular type of source ( Read )
 */

void mesonSeqSource(const LatticePropagator& quark_propagator,
		    LatticePropagator& seq_src_prop, 
		    int t_sink, multi1d<int>& sink_mom, 
		    int j_decay, int seq_src)
{
  LatticePropagator src_prop_tmp;

  START_CODE("mesonSeqSource");

  if ( Ns != 4 || Nc != 3 )		// Code is specific to Ns=4 and Nc=3
  {
    END_CODE("mesonSeqSource");
    return;
  }

  int G5 = Ns*Ns-1;
  
  switch (seq_src)
  {
  case 10:
  {
    /*
     *  First we consider the case of the pion
     */

    /*
     *  Begin by evaluating q2_tmp = gamma_5 * D * gamma_5
     *  where in this encoding gamma_5 = gamma(15)
     */

#if 0
    if( HalfInvertP == 1 ) { 
      QDP_error_exit("HaltInvertP is not yet implemented for mesonic sequential sources\n");
    }
#endif

    src_prop_tmp = Gamma(G5) * quark_propagator * Gamma(G5);
  }
  break;

  default:
    QDP_error_exit("Unknown sequential source type", seq_src);
  }


  /*
   *  We now inject momentum at sink if required
   */
  bool nonzero = false;
  for(int mu=0, j=0; mu < Nd; mu++)
  {
    if (mu != j_decay)
    {
      if(sink_mom[j] != 0)
      {
	nonzero = true;
	break;
      }
      j++;
    }
  }

  // multiply in the phase if required
  if (nonzero)
  {
    multi1d<LatticeInteger> my_coord(Nd);
    for(int mu=0; mu < Nd; mu++)
      my_coord[mu] = Layout::latticeCoordinate(mu);	/* Obtains the muth coordinate */

    LatticeReal p_dot_x = 0;
    for(int mu=0, j=0; mu < Nd; ++mu)
    {
      if (mu == j_decay)
        continue;

      p_dot_x += my_coord[mu] * sink_mom[j] * twopi / real(Layout::lattSize()[mu]);
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
        
  END_CODE("mesonSeqSource");
}
