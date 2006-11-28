// $Id: hadron_seqsource.cc,v 3.8 2006-11-28 19:28:57 edwards Exp $
/*! \file
 *  \brief Construct hadron sequential sources
 */

#include "util/ft/single_phase.h"
#include "meas/hadron/hadron_seqsource.h"
#include "util/ferm/gamma5_herm_w.h"

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
			    int t_sink, int j_decay)
    {
      START_CODE();

      if (j_decay < 0 || j_decay >= Nd)
      {
	QDPIO::cerr << __func__ << ": j_decay out of bounds" << endl;
	QDP_abort(1);
      }

      /*
       * Now mask out all but sink time slice
       */
      T seq_src_prop = where(Layout::latticeCoordinate(j_decay) == t_sink,
			     source_prop,
			     T(zero));
        
      END_CODE();

      return seq_src_prop;
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


  }


  // Default versions
  template<>
  Complex
  HadronSeqSource<LatticePropagator>::tieBack(const multi1d<LatticeColorMatrix>& u,
					      const SequentialProp_t& seqprop_header,
					      const LatticePropagator& seqprop, 
					      int gamma_insertion)
  {
    getTSrce() = hadSeqSourceGetTSrce(seqprop_header.forward_props);
    LatticeComplex tr = trace(gamma5Herm(seqprop) * Gamma(gamma_insertion));
    Complex seq_src_value = peekSite(tr, getTSrce());
    return seq_src_value;
  }


  // Default versions
  template<>
  LatticePropagator
  HadronSeqSource<LatticePropagator>::project(const LatticePropagator& src_prop_tmp) const
  {
    return hadSeqSourceProject<LatticePropagator>(src_prop_tmp,
						  getTSink(), getDecayDir());
  }


  // Default versions
  template<>
  LatticeComplex
  HadronSeqSource<LatticePropagator>::phases() const
  {
    return singlePhase(getTSrce(), getSinkMom(), getDecayDir());
  }


  // Default versions
  template<>
  void
  HadronSeqSource<LatticePropagator>::setTSrce(const multi1d<ForwardProp_t>& forward_headers)
  {
    getTSrce() = hadSeqSourceGetTSrce(forward_headers);
  }

  // Any other fermion types, like staggered can go here


}  // end namespace Chroma
