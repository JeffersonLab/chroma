// $Id: single_phase.cc,v 3.1 2006-11-27 20:07:58 edwards Exp $
/*! \file
 *  \brief Compute a single phase factor
 */

#include "chromabase.h"
#include "util/ft/single_phase.h"

namespace Chroma 
{

  //! A single exp(ip.x) phase used in hadron construction
  /*! @ingroup ft */
  LatticeComplex singlePhase(const multi1d<int>& t_srce, 
			     const multi1d<int>& sink_mom, 
			     int j_decay)
  {
    START_CODE();

    // Sanity checks
    if (j_decay < 0 || j_decay >= Nd)
    {
      QDPIO::cerr << __func__ << ": j_decay out of bounds" << endl;
      QDP_abort(1);
    }

    if (t_srce.size() != Nd)
    {
      QDPIO::cerr << __func__ << ": t_srce not of size = " << Nd << endl;
      QDP_abort(1);
    }

    if (sink_mom.size() != Nd-1)
    {
      QDPIO::cerr << __func__ << ": mom not of size = " << Nd-1 << endl;
      QDP_abort(1);
    }

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
	  p_dot_x += (Layout::latticeCoordinate(mu) - t_srce[mu]) * sink_mom[j] 
	    * twopi / Real(Layout::lattSize()[mu]);
	}
	j++;
      }
    }

    LatticeComplex phase;

    if (nonzeroP)
      phase = cmplx(cos(p_dot_x),sin(p_dot_x));
    else
      phase = 1.0;

    END_CODE();

    return phase;
  }

}  // end namespace Chroma
