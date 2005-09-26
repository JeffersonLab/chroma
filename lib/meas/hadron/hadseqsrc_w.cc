// $Id: hadseqsrc_w.cc,v 2.1 2005-09-26 04:48:35 edwards Exp $
/*! \file
 *  \brief Construct hadron sequential sources
 */

#include "chromabase.h"
#include "util/ft/sftmom.h"
#include "meas/hadron/seqsrc_funcmap_w.h"
#include "meas/hadron/hadseqsrc_w.h"

namespace Chroma 
{

  //! Construct hadron sequential sources
  /*!
   * \ingroup hadron
   *
   * This routine is specific to Wilson fermions!
   *
   * Construct hadronic sequential sources.
   *
   * \param quark_propagators    array of quark propagators ( Read )
   * \param seq_src_prop         sequential source as propagator ( Write )
   * \param t_sink               time coordinate of the sink ( Read )
   * \param sink_mom             sink baryon momentum ( Read )
   * \param j_decay              direction of the exponential decay ( Read )
   * \param seq_src_name         string name of the sequential source ( Read )
   *
   * \return Sequential source propagator
   */

  LatticePropagator hadSeqSource(const multi1d<LatticePropagator>& quark_propagators, 
				 int t_sink, const multi1d<int>& sink_mom, 
				 int j_decay, 
				 const string& seq_src_name)
  {
    START_CODE();

    LatticePropagator seq_src_prop;
    LatticePropagator src_prop_tmp;
  
    if ( Ns != 4 || Nc != 3 )		/* Code is specific to Ns=4 and Nc=3. */
    {
      END_CODE();
      QDPIO::cerr << __func__ << ": only Ns=4 and Nc=3 supported" << endl;
      QDP_abort(1);
    }

    try 
    {
      src_prop_tmp = 
	TheSeqSourceFuncMap::Instance().callFunction(seq_src_name,
						     quark_propagators);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << "Sequential source map call failed: " << e << endl;
      QDP_abort(1);
    }

    /* Now take hermitian conjugate and multiply on both sides
       with gamma_5 = Gamma(15) */
    {
      /* This bit does the hermitian conjugation. */
      LatticePropagator q1_tmp = adj(src_prop_tmp);
      src_prop_tmp = q1_tmp;
    }
   
    /* src_prop_tmp now holds src^{dagger} */
    /* If we are working with half inversions I need to do an 
       extra projection here */
    {
      /* Now slap the gamma_5 on either end */
      LatticePropagator q1_tmp = src_prop_tmp * Gamma(15);
      src_prop_tmp = Gamma(15) * q1_tmp;
    }
        
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

    return seq_src_prop;
  }

}  // end namespace Chroma
