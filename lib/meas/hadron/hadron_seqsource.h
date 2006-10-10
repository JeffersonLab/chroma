// -*- C++ -*-
// $Id: hadron_seqsource.h,v 3.3 2006-10-10 18:32:06 edwards Exp $ 
/*! \file
 *  \brief Construct hadron sequential sources
 */

#ifndef __hadron_seqsource_h__
#define __hadron_seqsource_h__

#include "chromabase.h"
#include "io/qprop_io.h"

namespace Chroma
{
  //! Construct hadron sequential sources
  /*! @ingroup hadron
   *
   * Supports creation of hadron sequential sources
   */
  template<typename T>
  class HadronSeqSource
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~HadronSeqSource() {}

    //! Construct the source
    virtual T operator()(const multi1d<LatticeColorMatrix>& u,
			 const multi1d<ForwardProp_t>& forward_headers,
			 const multi1d<T>& forward_props) const = 0;

  protected:
    //! Default method to put source on a fixed time slice and momenta
    /*!
     * \param src_prop_tmp     Sequential source before projection to a specific momenta ( Read )
     */
    virtual T project(const T& src_prop_tmp,
		      const multi1d<int>& t_srce, 
		      const multi1d<int>& sink_mom, 
		      int t_sink, int j_decay) const;

    //! Time-ordering phase of source and sink hadron states
    /*! Default is 1 */
    virtual Complex timeOrder(const multi1d<int>& t_srce, int t_sink, int j_decay) const
    {
      return Real(1);
    }

    //! Convenience function to yank the source location from the forward prop headers
    virtual multi1d<int> getTSrce(const multi1d<ForwardProp_t>& forward_headers) const;
  };

}


#endif
