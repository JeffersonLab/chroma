// -*- C++ -*-
// $Id: derivmesonseqsrc_w.h,v 2.2 2006-02-09 03:50:06 edwards Exp $
/*! \file
 *  \brief Construct derivative meson sequential sources.
 */

#ifndef __derivmesonseqsrc_w_h__
#define __derivmesonseqsrc_w_h__

#include "meas/hadron/hadron_seqsource.h"

namespace Chroma 
{

  //! Name and registration
  /*! @ingroup hadron */
  namespace DerivMesonSeqSourceEnv
  {
    extern const bool registered;

  
    //! Deriv meson sequential source parameters
    /*! @ingroup hadron */
    struct Params
    {
      Params();
      Params(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;

      int              deriv_length;    /*!< Displacement length in derivative */

      multi1d<int>     sink_mom;        /*!< sink momentum */
      int              t_sink;          /*!< time slice of sink */
      int              j_decay;         /*!< Decay direction */
    };

    //! Deriv meson sequential source parameters
    /*! @ingroup hadron */
    struct ParamsDir
    {
      ParamsDir();
      ParamsDir(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;

      int              deriv_dir;       /*!< Polarization direction */
      int              deriv_length;    /*!< Displacement length in derivative */

      multi1d<int>     sink_mom;        /*!< sink momentum */
      int              t_sink;          /*!< time slice of sink */
      int              j_decay;         /*!< Decay direction */
    };


    //! Construct pion_1-(A2=A1xD_T2) sequential source
    /*!
     * \ingroup hadron
     *
     * Comes from a1xD_T2 in Manke's hep-lat/0210030 .
     * The sink interpolator is   \f$\Gamma_f \equiv \gamma_5\epsilon_{ijk}\gamma_j D_k\f$  
     * where   \f$D_i = s_{ijk}\nabla_j\nabla_k\f$
     * and     \f$\nabla_\mu f(x) = U_\mu(x)f(x+\mu) - U_{-\mu}(x)f(x-\mu)\f$
     */
    class MesPionA1xDT2SeqSrc : public HadronSeqSource<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionA1xDT2SeqSrc(const ParamsDir& p) : params(p) {}

      //! Default destructor
      ~MesPionA1xDT2SeqSrc() {}
      
      //! Construct pion_1-(A2=A1xD_T2) sequential source
      /*!
       * \ingroup hadron
       *
       * \param forward_props    array of quark propagators ( Read )
       *
       * \return \f$F^\dagger \gamma_5 \Gamma_f\f$
       */
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<LatticePropagator>& forward_props) const;

    private:
      //! Hide partial constructor
      MesPionA1xDT2SeqSrc() {}

    private:
      ParamsDir  params;   /*!< Seqsource params */
    };

  }  // end namespace


  //! Reader
  /*! @ingroup hadron */
  void read(XMLReader& xml, const string& path, DerivMesonSeqSourceEnv::Params& param);

  //! Writer
  /*! @ingroup hadron */
  void write(XMLWriter& xml, const string& path, const DerivMesonSeqSourceEnv::Params& param);


  //! Reader
  /*! @ingroup hadron */
  void read(XMLReader& xml, const string& path, DerivMesonSeqSourceEnv::ParamsDir& param);

  //! Writer
  /*! @ingroup hadron */
  void write(XMLWriter& xml, const string& path, const DerivMesonSeqSourceEnv::ParamsDir& param);


}  // end namespace Chroma

#endif
