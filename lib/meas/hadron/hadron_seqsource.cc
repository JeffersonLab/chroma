// $Id: hadron_seqsource.cc,v 3.0 2006-04-03 04:58:59 edwards Exp $
/*! \file
 *  \brief Construct hadron sequential sources
 */

#include "chromabase.h"
#include "meas/hadron/hadron_seqsource.h"

namespace Chroma 
{

  //! Anonymous namespace
  /*! @ingroup hadron */
  namespace
  {
    //! Construct hadron sequential sources
    /*! @ingroup hadron */
    template<typename T> 
    T   hadSeqSourceProject(const T& source_prop,
			    const multi1d<int>& sink_mom, 
			    int t_sink, int j_decay)
    {
      START_CODE();

      if (j_decay < 0 || j_decay >= Nd)
      {
	QDPIO::cerr << __func__ << ": j_decay out of bounds" << endl;
	QDP_abort(1);
      }

      /* Now take hermitian conjugate and multiply on both sides
	 with gamma_5 = Gamma(15) */
      T seq_src_tmp = Gamma(15) * adj(source_prop) * Gamma(15);
        
      /*
       *  We now inject momentum at sink if required
       */
      LatticeReal p_dot_x = zero;
      bool nonzeroP = false;
      for(int mu=0, j=0; mu < Nd; mu++)
      {
	if (mu != j_decay)
	{
	  if (sink_mom[j] != 0)
	  {
	    nonzeroP = true;
	    p_dot_x += Layout::latticeCoordinate(mu) * sink_mom[j] 
	      * twopi / Real(Layout::lattSize()[mu]);
	  }
	  j++;
	}
      }

      if (nonzeroP)
	seq_src_tmp *= cmplx(cos(p_dot_x),sin(p_dot_x));

      /*
       * Now mask out all but sink time slice
       */
      T seq_src_prop = where(Layout::latticeCoordinate(j_decay) == t_sink,
			     seq_src_tmp,
			     LatticePropagator(zero));
        
      END_CODE();

      return seq_src_prop;
    }
  }


  // Default versions
  template<>
  LatticePropagator
  HadronSeqSource<LatticePropagator>::project(const LatticePropagator& src_prop_tmp,
					      const multi1d<int>& sink_mom, 
					      int t_sink, int j_decay) const
  {
    return hadSeqSourceProject<LatticePropagator>(src_prop_tmp,
						  sink_mom,
						  t_sink, j_decay);
  }

  // Any other fermion types, like staggered can go here


}  // end namespace Chroma
