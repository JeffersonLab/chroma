// $Id: hadron_seqsource.cc,v 3.5 2006-10-11 13:27:07 edwards Exp $
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
			    const multi1d<int>& t_srce, 
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
	    p_dot_x += (Layout::latticeCoordinate(mu) - t_srce[mu]) * sink_mom[j] 
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


  //! Get source location
  multi1d<int>
  hadSeqSourceGetTSrce(const multi1d<ForwardProp_t>& forward_headers)
  {
    multi1d<int> t_srce = forward_headers[0].source_header.getTSrce();

    for(int loop=1; loop < forward_headers.size(); ++loop)
    {
      multi1d<int> t_srce_b = forward_headers[loop].source_header.getTSrce();

      // Bummer, I wish qdp++ had a multi1d.operator!=()
      bool same = true;
      for(int i=0; i < t_srce.size(); ++i)
      {
	if (t_srce_b[i] != t_srce[i]) 
	  same = false;
      }
      
      if (! same)
      {
	QDPIO::cerr << __func__ << ": the t_srce in the forward props are not all equal"
		    << endl;
	QDP_abort(1);
      }
    }

    return t_srce;
  }


  //! Get fermion bc
  multi1d<int>
  hadSeqSourceGetBC(const multi1d<ForwardProp_t>& forward_headers)
  {
    multi1d<int> bc = getFermActBoundary(forward_headers[0].prop_header.fermact);

    for(int loop=1; loop < forward_headers.size(); ++loop)
    {
      multi1d<int> bc_b = getFermActBoundary(forward_headers[loop].prop_header.fermact);
    
      // Bummer, I wish qdp++ had a multi1d.operator!=()
      bool same = true;
      for(int i=0; i < bc.size(); ++i)
      {
	if (bc_b[i] != bc[i]) 
	  same = false;
      }
      
      if (! same)
      {
	QDPIO::cerr << __func__ << ": the bc in the forward props are not all equal"
		    << endl;
	QDP_abort(1);
      }
    }

    return bc;
  }




  // Default versions
  template<>
  LatticePropagator
  HadronSeqSource<LatticePropagator>::project(const LatticePropagator& src_prop_tmp,
					      const multi1d<int>& t_srce, 
					      const multi1d<int>& sink_mom, 
					      int t_sink, int j_decay) const
  {
    return hadSeqSourceProject<LatticePropagator>(src_prop_tmp,
						  t_srce,
						  sink_mom,
						  t_sink, j_decay);
  }


  // Default versions
  template<>
  multi1d<int>
  HadronSeqSource<LatticePropagator>::getTSrce(const multi1d<ForwardProp_t>& forward_headers) const
  {
    return hadSeqSourceGetTSrce(forward_headers);
  }


  // Default versions
  template<>
  multi1d<int>
  HadronSeqSource<LatticePropagator>::getBC(const multi1d<ForwardProp_t>& forward_headers) const
  {
    return hadSeqSourceGetBC(forward_headers);
  }


  // Any other fermion types, like staggered can go here


}  // end namespace Chroma
