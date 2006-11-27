// -*- C++ -*-
// $Id: hadron_seqsource.h,v 3.5 2006-11-27 04:33:35 edwards Exp $ 
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
			 const multi1d<T>& forward_props) = 0;

  protected:
    //! Project onto a definite time-slice
    virtual T project(const LatticePropagator& src_prop_tmp) const;

    //! Construct phases
    virtual LatticeComplex phases() const;

    //! Convenience function to yank the source location from the forward prop headers
    virtual void setTSrce(const multi1d<ForwardProp_t>& forward_headers);

    //! Set t_srce
    virtual multi1d<int>& getTSrce() = 0;

    //! Get t_srce
    virtual const multi1d<int>& getTSrce() const = 0;

    //! Get t_sink
    virtual int getTSink() const = 0;

    //! Get sink_mom
    virtual const multi1d<int>& getSinkMom() const = 0;

    //! Get decay_dir
    virtual const int getDecayDir() const = 0;
  };

}


#endif
