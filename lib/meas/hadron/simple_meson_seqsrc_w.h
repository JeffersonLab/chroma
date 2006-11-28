// -*- C++ -*-
// $Id: simple_meson_seqsrc_w.h,v 3.2 2006-11-28 19:28:57 edwards Exp $
/*! \file
 *  \brief Construct meson sequential sources.
 */

#ifndef __simple_meson_seqsrc_w_h__
#define __simple_meson_seqsrc_w_h__

#include "meas/hadron/meson_seqsrc_w.h"

namespace Chroma 
{

  //! Name and registration
  /*! @ingroup hadron */
  namespace SimpleMesonSeqSourceEnv
  {
    bool registerAll();

  
    //! Simple meson sequential source parameters
    /*! @ingroup hadron */
    struct Params
    {
      Params();
      Params(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;

      multi1d<int>     sink_mom;        /*!< sink momentum */
      int              t_sink;          /*!< time slice of sink */
      int              j_decay;         /*!< Decay direction */
    };


    //! Simple meson sequential source construction
    /*! @ingroup hadron
     *
     * Create a simple meson sequential propagator source
     */
    class SimpleMesonSeqSource : public MesonSeqSourceBase
    {
    public:
      //! Full constructor
      SimpleMesonSeqSource(const Params& p, int gamma) : params(p), gamma_sink(gamma) {}

      //! Default destructor
      ~SimpleMesonSeqSource() {}
      
      //! Construct the source
      /*!
       * \param u                    gauge field ( Read )
       * \param forward_props        array of quark propagators ( Read )
       * \param forward_headers      corresponding array of quark propagators ( Read )
       *
       * \return \f$\gamma_5 * \Gamma(gamma_sink)^dag * \gamma_5 * F\$
       */
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);

    protected:
      //! Set t_srce
      multi1d<int>& getTSrce() {return t_srce;}

      //! Get t_srce
      const multi1d<int>& getTSrce() const {return t_srce;}

      //! Get t_sink
      int getTSink() const {return params.t_sink;}

      //! Get sink_mom
      const multi1d<int>& getSinkMom() const {return params.sink_mom;}

      //! Get decay_dir
      const int getDecayDir() const {return params.j_decay;}

    private:
      //! Hide partial constructor
      SimpleMesonSeqSource() {}

    private:
      multi1d<int>  t_srce;   /*<! Must come from propagator headers */
      Params  params;   /*!< Seqsource params */
      int gamma_sink;   /*!< The gamma matrix at the sink */
    };


    //! Simple meson sequential source construction
    /*! @ingroup hadron
     *
     * Create a simple meson sequential propagator source
     */
    class PionPionSeqSource : public MesonSeqSourceBase
    {
    public:
      //! Full constructor
      PionPionSeqSource(const Params& p) : params(p) {}

      //! Default destructor
      ~PionPionSeqSource() {}
      
      //! Construct the source
      /*!
       * \param u                    gauge field ( Read )
       * \param forward_props        array of quark propagators ( Read )
       * \param forward_headers      corresponding array of quark propagators ( Read )
       *
       * \return \f$\gamma_5 * \Gamma(15)^dag * \gamma_5 * F * \gamma_5\$
       */
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);

    protected:
      //! Set t_srce
      multi1d<int>& getTSrce() {return t_srce;}

      //! Get t_srce
      const multi1d<int>& getTSrce() const {return t_srce;}

      //! Get t_sink
      int getTSink() const {return params.t_sink;}

      //! Get sink_mom
      const multi1d<int>& getSinkMom() const {return params.sink_mom;}

      //! Get decay_dir
      const int getDecayDir() const {return params.j_decay;}

    private:
      //! Hide partial constructor
      PionPionSeqSource() {}

    private:
      multi1d<int>  t_srce;   /*<! Must come from propagator headers */
      Params  params;   /*!< Seqsource params */
    };

  }  // end namespace


  //! Reader
  /*! @ingroup hadron */
  void read(XMLReader& xml, const string& path, SimpleMesonSeqSourceEnv::Params& param);

  //! Writer
  /*! @ingroup hadron */
  void write(XMLWriter& xml, const string& path, const SimpleMesonSeqSourceEnv::Params& param);


}  // end namespace Chroma

#endif
