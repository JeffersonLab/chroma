// $Id: baryon_seqsrc_w.cc,v 3.1 2006-11-27 04:33:35 edwards Exp $
/*! \file
 *  \brief Construct baryon sequential sources.
 */

#include "meas/hadron/baryon_seqsrc_w.h"
#include "meas/hadron/seqsource_factory_w.h"

namespace Chroma 
{


  //! Anonymous namespace
  /*! @ingroup hadron */
  namespace
  {
    //! Construct hadron sequential sources
    /*! @ingroup hadron */
    template<typename T> 
    T   barSeqSourceProject(const T& source_prop,
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
       * Now mask out all but sink time slice
       */
      T seq_src_prop = where(Layout::latticeCoordinate(j_decay) == t_sink,
			     seq_src_tmp,
			     LatticePropagator(zero));
        
      END_CODE();

      return seq_src_prop;
    }

    //! Get fermion bc
    multi1d<int>
    barSeqSourceGetBC(const multi1d<ForwardProp_t>& forward_headers)
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

  }



  // Default versions
  void
  BaryonSeqSourceBase::setBC(const multi1d<ForwardProp_t>& forward_headers)
  {
    getBC() = barSeqSourceGetBC(forward_headers);
  }



  // Combine projection with time-ordering
  LatticePropagator 
  BaryonSeqSourceBase::projectBaryon(const LatticePropagator& src_prop_tmp,
				     const multi1d<ForwardProp_t>& forward_headers)
  {
    setTSrce(forward_headers);
    setBC(forward_headers);

    // Multiply in time-ordering phase and project
    return timeOrder() * phases() * project(src_prop_tmp);
  }


  // Time-ordering phase of source and sink hadron states
  Complex BaryonSeqSourceBase::timeOrder() const
  {
    Complex phase;

    int j_decay  = getDecayDir();
    int t_sink   = getTSink();
    int t_source = getTSrce()[j_decay];
    int bc_spec  = getBC()[j_decay];

    if ( (bc_spec < 0) && (t_source > t_sink) )
      phase = -1;
    else
      phase = 1;

    return phase;
  }

}  // end namespace Chroma
