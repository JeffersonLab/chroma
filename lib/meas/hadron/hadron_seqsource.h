// -*- C++ -*-
// $Id: hadron_seqsource.h,v 2.1 2006-02-09 02:25:25 edwards Exp $ 
/*! \file
 *  \brief Construct hadron sequential sources
 */

#ifndef __hadron_seqsource_h__
#define __hadron_seqsource_h__

#include "chromabase.h"

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
			 const multi1d<T>& forward_props) const = 0;

  protected:
    //! Default method to put source on a fixed time slice and momenta
    /*!
     * \param src_prop_tmp     Sequential source before projection to a specific momenta ( Read )
     */
    virtual T project(const T& src_prop_tmp,
		      const multi1d<int>& sink_mom, 
		      int t_sink, int j_decay) const;

  };

}


#endif
