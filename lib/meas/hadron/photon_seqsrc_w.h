// -*- C++ -*-
// $Id: photon_seqsrc_w.h,v 3.0 2006-04-03 04:59:00 edwards Exp $
/*! \file
 *  \brief Construct a photon sequential sources via LSZ reduction
 */

#ifndef __photon_seqsrc_w_h__
#define __photon_seqsrc_w_h__

#include "meas/hadron/hadron_seqsource.h"

namespace Chroma 
{

  //! Name and registration
  /*! @ingroup hadron */
  namespace PhotonRhoSeqSourceEnv
  {
    extern const bool registered;

  
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
    class PhotonRhoSeqSource : public HadronSeqSource<LatticePropagator>
    {
    public:
      //! Full constructor
      PhotonRhoSeqSource(const Params& p) : params(p) {}

      //! Construct the source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<LatticePropagator>& forward_props) const;

    private:
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
