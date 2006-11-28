// -*- C++ -*-
// $Id: hadron_seqsource.h,v 3.6 2006-11-28 19:28:57 edwards Exp $ 
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

    //! Evaluate the seqprop back at the source - this is the 2-pt at the sink
    /*!
     *  For the case of a meson, we have evaluated as the sequential source
     *
     *  H(y, 0; tx, p) = \sum exp{ip.x} U(y,x) \gamma_5\Gamma_f^\dag\gamma_5 D(x,0) 
     *
     *  H^\dag(y, 0; tx, p) = \sum_x exp{-ip.x} \gamma_5 D(0,x) \Gamma_f U(x,y) \gamma_5
     *
     *  Thus we can see that 
     *
     *  Tr[ \gamma_5 H^\dag(0,0; tx, p)\gamma_5 \Gamma_i] = 
     *                         \sum_x exp{-ip.x} Tr[ D(0,x)\Gamma_f U(x,0) \Gamma_i ]
     *
     *  which is the desired meson correlator at momentum p and timslice tx
     */
    virtual Complex tieBack(const multi1d<LatticeColorMatrix>& u,
			    const SequentialProp_t& seqprop_header,
			    const T& seqprop, 
			    int gamma_insertion);

    //! Compute the 2-pt at the sink
    virtual Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			      const multi1d<ForwardProp_t>& forward_headers,
			      const multi1d<T>& forward_props,
			      int gamma_insertion) = 0;

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
