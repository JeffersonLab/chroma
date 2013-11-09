// -*- C++ -*-
// $Id: photon_seqsrc_w.h,v 3.5 2006-11-28 19:28:57 edwards Exp $
/*! \file
 *  \brief Construct a photon sequential sources via LSZ reduction
 */

#ifndef __photon_seqsrc_w_h__
#define __photon_seqsrc_w_h__

#include "meas/hadron/meson_seqsrc_w.h"

namespace Chroma 
{

  //! Name and registration
  /*! @ingroup hadron */
  namespace PhotonRhoSeqSourceEnv
  {
    bool registerAll();

  
    //! Construct a photon sequential sources via LSZ reduction
    /*! @ingroup hadron */
    struct Params
    {
      Params();
      Params(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;

      int              pol_dir;         /*!< Polarization direction */
      int              t_sink_start;    /*!< time slice of sink (inclusive) to start integration */
      int              t_sink_end;      /*!< time slice of sink (inclusive) to end integration */
      Real             Q_sq;            /*!< Photon virtuality, \f$Q_f^2 = c^2|\vec{p_f}|^2 - \omega_f^2 */
      Real             c_sq;            /*!< Photon speed of light */
      Real             xi;              /*!< Renormalized anisotropy */

      multi1d<int>     sink_mom;        /*!< sink momentum */
      int              j_decay;         /*!< Decay direction */

      int              t_sink;          /*!< will always be -1 here */
    };


    //! Construct a photon sequential sources via LSZ reduction
    /*! @ingroup hadron
     *
     *  Photon source via a rho
     */
    class PhotonRhoSeqSource : public MesonSeqSourceBase
    {
    public:
      //! Full constructor
      PhotonRhoSeqSource(const Params& p) : params(p) {}

      //! Construct the source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion)
	{
	  QDPIO::cerr << __func__ << ": not implemented" << endl;
	  QDP_abort(1);
	  return 0.0;
	}

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
      multi1d<int>  t_srce;   /*<! Must come from propagator headers */
      Params  params;   /*!< Seqsource params */
    };



    //! Construct a photon sequential sources via LSZ reduction
    /*! @ingroup hadron
     *
     *  Photon source via a point-split rho
     */
    class PointSplitPhotonRhoSeqSource : public MesonSeqSourceBase
    {
    public:
      //! Full constructor
      PointSplitPhotonRhoSeqSource(const Params& p) : params(p) {}

      //! Construct the source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion)
	{
	  QDPIO::cerr << __func__ << ": not implemented" << endl;
	  QDP_abort(1);
	  return 0.0;
	}

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
      multi1d<int>  t_srce;   /*<! Must come from propagator headers */
      Params  params;   /*!< Seqsource params */
    };

  }  // end namespace


  //! Reader
  /*! @ingroup hadron */
  void read(XMLReader& xml, const string& path, PhotonRhoSeqSourceEnv::Params& param);

  //! Writer
  /*! @ingroup hadron */
  void write(XMLWriter& xml, const string& path, const PhotonRhoSeqSourceEnv::Params& param);


}  // end namespace Chroma

#endif
