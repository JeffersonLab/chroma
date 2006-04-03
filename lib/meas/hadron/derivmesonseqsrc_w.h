// -*- C++ -*-
// $Id: derivmesonseqsrc_w.h,v 3.0 2006-04-03 04:58:59 edwards Exp $
/*! \file
 *  \brief Construct derivative meson sequential sources.
 *
 *  These operators come from Liao and Manke, hep-lat/0210030 .
 *  There are a few variants of derivatives. They are
 *
 *  \f$D_i = s_{ijk}\nabla_j\nabla_k\f$
 *  \f$B_i = \epsilon_{ijk}\nabla_j\nabla_k\f$
 *  \f$\nabla_\mu f(x) = U_\mu(x)f(x+\mu) - U_{-\mu}(x)f(x-\mu)\f$
 *
 *  where
 *   
 *  \f$s_{ijk} = |\epsilon_{ijk}|\f$
 *  \f$S_{\alpha jk} = 0\quad j\ne k, S_{111}=S_{222}=1, S_{122}=S_{233}=-1\f$
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


    //! Construct pion_1-(PionxNabla_T1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  Pion x Nabla_T1
     * The sink interpolator structure is
     * \f$\Gamma_f \equiv \gamma_5\nabla_i\f$
     */
    class MesPionPionxNablaT1SeqSrc : public HadronSeqSource<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionPionxNablaT1SeqSrc(const ParamsDir& p) : params(p) {}

      //! Default destructor
      ~MesPionPionxNablaT1SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<LatticePropagator>& forward_props) const;

    private:
      //! Hide partial constructor
      MesPionPionxNablaT1SeqSrc() {}

    private:
      ParamsDir  params;   /*!< Seqsource params */
    };


    //! Construct pion_1-(A0xNabla_T1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  a0 x nabla_T1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \nabla_i\f$
     */
    class MesPionA0xNablaT1SeqSrc : public HadronSeqSource<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionA0xNablaT1SeqSrc(const ParamsDir& p) : params(p) {}

      //! Default destructor
      ~MesPionA0xNablaT1SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<LatticePropagator>& forward_props) const;

    private:
      //! Hide partial constructor
      MesPionA0xNablaT1SeqSrc() {}

    private:
      ParamsDir  params;   /*!< Seqsource params */
    };


    //! Construct pion_1-(A0_2xNabla_T1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  a0_2 x nabla_T1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_4 \nabla_i\f$
     */
    class MesPionA02xNablaT1SeqSrc : public HadronSeqSource<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionA02xNablaT1SeqSrc(const ParamsDir& p) : params(p) {}

      //! Default destructor
      ~MesPionA02xNablaT1SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<LatticePropagator>& forward_props) const;

    private:
      //! Hide partial constructor
      MesPionA02xNablaT1SeqSrc() {}

    private:
      ParamsDir  params;   /*!< Seqsource params */
    };


    //! Construct pion_1-(RhoxNabla_A1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  rho x nabla_A1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_i\nabla_i\f$
     */
    class MesPionRhoxNablaA1SeqSrc : public HadronSeqSource<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionRhoxNablaA1SeqSrc(const Params& p) : params(p) {}

      //! Default destructor
      ~MesPionRhoxNablaA1SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<LatticePropagator>& forward_props) const;

    private:
      //! Hide partial constructor
      MesPionRhoxNablaA1SeqSrc() {}

    private:
      Params  params;   /*!< Seqsource params */
    };


    //! Construct pion_1-(RhoxNabla_T1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  rho x nabla_T1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \epsilon_{ijk}\gamma_j \nabla_k\f$  
     */
    class MesPionRhoxNablaT1SeqSrc : public HadronSeqSource<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionRhoxNablaT1SeqSrc(const ParamsDir& p) : params(p) {}

      //! Default destructor
      ~MesPionRhoxNablaT1SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<LatticePropagator>& forward_props) const;

    private:
      //! Hide partial constructor
      MesPionRhoxNablaT1SeqSrc() {}

    private:
      ParamsDir  params;   /*!< Seqsource params */
    };


    //! Construct pion_1-(RhoxNabla_T2) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  rho x nabla_T2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv s_{ijk}\gamma_j D_k\f$  
     */
    class MesPionRhoxNablaT2SeqSrc : public HadronSeqSource<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionRhoxNablaT2SeqSrc(const ParamsDir& p) : params(p) {}

      //! Default destructor
      ~MesPionRhoxNablaT2SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<LatticePropagator>& forward_props) const;

    private:
      //! Hide partial constructor
      MesPionRhoxNablaT2SeqSrc() {}

    private:
      ParamsDir  params;   /*!< Seqsource params */
    };


    //! Construct pion_1-(A1xNabla_A1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  a1 x D_A1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_5\gamma_i \nabla_i\f$  
     */
    class MesPionA1xNablaA1SeqSrc : public HadronSeqSource<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionA1xNablaA1SeqSrc(const Params& p) : params(p) {}

      //! Default destructor
      ~MesPionA1xNablaA1SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<LatticePropagator>& forward_props) const;

    private:
      //! Hide partial constructor
      MesPionA1xNablaA1SeqSrc() {}

    private:
      Params  params;   /*!< Seqsource params */
    };


    //! Construct pion_1-(A1xNabla_T2) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  a1 x nabla_T2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_5 s_{ijk}\gamma_j \nabla_k\f$  
     */
    class MesPionA1xNablaT2SeqSrc : public HadronSeqSource<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionA1xNablaT2SeqSrc(const ParamsDir& p) : params(p) {}

      //! Default destructor
      ~MesPionA1xNablaT2SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<LatticePropagator>& forward_props) const;

    private:
      //! Hide partial constructor
      MesPionA1xNablaT2SeqSrc() {}

    private:
      ParamsDir  params;   /*!< Seqsource params */
    };


    //! Construct pion_1-(A1xNabla_E) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  a1 x nabla_E
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_5 S_{\alpha jk}\gamma_j \nabla_k\f$  
     */
    class MesPionA1xNablaESeqSrc : public HadronSeqSource<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionA1xNablaESeqSrc(const ParamsDir& p) : params(p) {}

      //! Default destructor
      ~MesPionA1xNablaESeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<LatticePropagator>& forward_props) const;

    private:
      //! Hide partial constructor
      MesPionA1xNablaESeqSrc() {}

    private:
      ParamsDir  params;   /*!< Seqsource params */
    };


    //! Construct pion_1-(B1xNabla_T1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  b1 x nabla_T1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_4\gamma_5\epsilon_{ijk}\gamma_j \nabla_k\f$  
     */
    class MesPionB1xNablaT1SeqSrc : public HadronSeqSource<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionB1xNablaT1SeqSrc(const ParamsDir& p) : params(p) {}

      //! Default destructor
      ~MesPionB1xNablaT1SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<LatticePropagator>& forward_props) const;

    private:
      //! Hide partial constructor
      MesPionB1xNablaT1SeqSrc() {}

    private:
      ParamsDir  params;   /*!< Seqsource params */
    };


    //! Construct pion_1-(A0_2xD_T2) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  a0_2 x D_T2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_4 D_i\f$  
     */
    class MesPionA02xDT2SeqSrc : public HadronSeqSource<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionA02xDT2SeqSrc(const ParamsDir& p) : params(p) {}

      //! Default destructor
      ~MesPionA02xDT2SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<LatticePropagator>& forward_props) const;

    private:
      //! Hide partial constructor
      MesPionA02xDT2SeqSrc() {}

    private:
      ParamsDir  params;   /*!< Seqsource params */
    };


    //! Construct pion_1-(A1xD_A2) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  a1 x D_A2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_5\gamma_i D_i\f$  
     */
    class MesPionA1xDA2SeqSrc : public HadronSeqSource<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionA1xDA2SeqSrc(const Params& p) : params(p) {}

      //! Default destructor
      ~MesPionA1xDA2SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<LatticePropagator>& forward_props) const;

    private:
      //! Hide partial constructor
      MesPionA1xDA2SeqSrc() {}

    private:
      Params  params;   /*!< Seqsource params */
    };


    //! Construct pion_1-(A1xD_E) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  a1 x D_E
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_5 S_{\alpha jk}\gamma_j D_k\f$  
     */
    class MesPionA1xDESeqSrc : public HadronSeqSource<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionA1xDESeqSrc(const ParamsDir& p) : params(p) {}

      //! Default destructor
      ~MesPionA1xDESeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<LatticePropagator>& forward_props) const;

    private:
      //! Hide partial constructor
      MesPionA1xDESeqSrc() {}

    private:
      ParamsDir  params;   /*!< Seqsource params */
    };


    //! Construct pion_1-(A1xD_T1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  a1 x D_T1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_5 s_{ijk}\gamma_j D_k\f$  
     */
    class MesPionA1xDT1SeqSrc : public HadronSeqSource<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionA1xDT1SeqSrc(const ParamsDir& p) : params(p) {}

      //! Default destructor
      ~MesPionA1xDT1SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<LatticePropagator>& forward_props) const;

    private:
      //! Hide partial constructor
      MesPionA1xDT1SeqSrc() {}

    private:
      ParamsDir  params;   /*!< Seqsource params */
    };


    //! Construct pion_1-(A1xD_T2) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  a1 x D_T2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_5\epsilon_{ijk}\gamma_j D_k\f$  
     */
    class MesPionA1xDT2SeqSrc : public HadronSeqSource<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionA1xDT2SeqSrc(const ParamsDir& p) : params(p) {}

      //! Default destructor
      ~MesPionA1xDT2SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<LatticePropagator>& forward_props) const;

    private:
      //! Hide partial constructor
      MesPionA1xDT2SeqSrc() {}

    private:
      ParamsDir  params;   /*!< Seqsource params */
    };


    //! Construct pion_1-(B1xD_A2) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  b1 x D_A2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_4\gamma_5 \gamma_i D_i\f$  
     */
    class MesPionB1xDA2SeqSrc : public HadronSeqSource<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionB1xDA2SeqSrc(const Params& p) : params(p) {}

      //! Default destructor
      ~MesPionB1xDA2SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<LatticePropagator>& forward_props) const;

    private:
      //! Hide partial constructor
      MesPionB1xDA2SeqSrc() {}

    private:
      Params  params;   /*!< Seqsource params */
    };

    //! Construct pion_1-(B1xD_E) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  b1 x D_E
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_4\gamma_5 S_{\alpha jk}\gamma_j D_k\f$  
     */
    class MesPionB1xDESeqSrc : public HadronSeqSource<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionB1xDESeqSrc(const ParamsDir& p) : params(p) {}

      //! Default destructor
      ~MesPionB1xDESeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<LatticePropagator>& forward_props) const;

    private:
      //! Hide partial constructor
      MesPionB1xDESeqSrc() {}

    private:
      ParamsDir  params;   /*!< Seqsource params */
    };


    //! Construct pion_1-(B1xD_T1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  b1 x D_T1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_4\gamma_5 s_{ijk}\gamma_j D_k\f$  
     */
    class MesPionB1xDT1SeqSrc : public HadronSeqSource<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionB1xDT1SeqSrc(const ParamsDir& p) : params(p) {}

      //! Default destructor
      ~MesPionB1xDT1SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<LatticePropagator>& forward_props) const;

    private:
      //! Hide partial constructor
      MesPionB1xDT1SeqSrc() {}

    private:
      ParamsDir  params;   /*!< Seqsource params */
    };


    //! Construct pion_1-(B1xD_T2) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  b1 x D_T2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_4\gamma_5 \epsilon_{ijk}\gamma_j D_k\f$  
     */
    class MesPionB1xDT2SeqSrc : public HadronSeqSource<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionB1xDT2SeqSrc(const ParamsDir& p) : params(p) {}

      //! Default destructor
      ~MesPionB1xDT2SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<LatticePropagator>& forward_props) const;

    private:
      //! Hide partial constructor
      MesPionB1xDT2SeqSrc() {}

    private:
      ParamsDir  params;   /*!< Seqsource params */
    };


    //! Construct pion_1-(RhoxD_A2) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  rho x D_A2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_i D_i\f$  
     */
    class MesPionRhoxDA2SeqSrc : public HadronSeqSource<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionRhoxDA2SeqSrc(const Params& p) : params(p) {}

      //! Default destructor
      ~MesPionRhoxDA2SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<LatticePropagator>& forward_props) const;

    private:
      //! Hide partial constructor
      MesPionRhoxDA2SeqSrc() {}

    private:
      Params  params;   /*!< Seqsource params */
    };


    //! Construct pion_1-(RhoxD_T1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  rho x D_T1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv s_{ijk}\gamma_j D_k\f$  
     */
    class MesPionRhoxDT1SeqSrc : public HadronSeqSource<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionRhoxDT1SeqSrc(const ParamsDir& p) : params(p) {}

      //! Default destructor
      ~MesPionRhoxDT1SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<LatticePropagator>& forward_props) const;

    private:
      //! Hide partial constructor
      MesPionRhoxDT1SeqSrc() {}

    private:
      ParamsDir  params;   /*!< Seqsource params */
    };


    //! Construct pion_1-(RhoxD_T2) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  rho x D_T2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \epsilon_{ijk}\gamma_j D_k\f$  
     */
    class MesPionRhoxDT2SeqSrc : public HadronSeqSource<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionRhoxDT2SeqSrc(const ParamsDir& p) : params(p) {}

      //! Default destructor
      ~MesPionRhoxDT2SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<LatticePropagator>& forward_props) const;

    private:
      //! Hide partial constructor
      MesPionRhoxDT2SeqSrc() {}

    private:
      ParamsDir  params;   /*!< Seqsource params */
    };


    //! Construct pion_1-(PionxD_T2) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  pion x D_T2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_4\gamma_5 D_i\f$  
     */
    class MesPionPionxDT2SeqSrc : public HadronSeqSource<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionPionxDT2SeqSrc(const ParamsDir& p) : params(p) {}

      //! Default destructor
      ~MesPionPionxDT2SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<LatticePropagator>& forward_props) const;

    private:
      //! Hide partial constructor
      MesPionPionxDT2SeqSrc() {}

    private:
      ParamsDir  params;   /*!< Seqsource params */
    };


    //! Construct pion_1-(PionxB_T1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  pion x B_T1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_5 B_i\f$  
     */
    class MesPionPionxBT1SeqSrc : public HadronSeqSource<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionPionxBT1SeqSrc(const ParamsDir& p) : params(p) {}

      //! Default destructor
      ~MesPionPionxBT1SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<LatticePropagator>& forward_props) const;

    private:
      //! Hide partial constructor
      MesPionPionxBT1SeqSrc() {}

    private:
      ParamsDir  params;   /*!< Seqsource params */
    };


    //! Construct pion_1-(RhoxB_T1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is rho x B_T1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \epsilon_{ijk}\gamma_j B_k\f$  
     */
    class MesPionRhoxBT1SeqSrc : public HadronSeqSource<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionRhoxBT1SeqSrc(const ParamsDir& p) : params(p) {}

      //! Default destructor
      ~MesPionRhoxBT1SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<LatticePropagator>& forward_props) const;

    private:
      //! Hide partial constructor
      MesPionRhoxBT1SeqSrc() {}

    private:
      ParamsDir  params;   /*!< Seqsource params */
    };


    //! Construct pion_1-(RhoxB_T2) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  rho x B_T2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv s_{ijk}\gamma_j D_k\f$  
     */
    class MesPionRhoxBT2SeqSrc : public HadronSeqSource<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionRhoxBT2SeqSrc(const ParamsDir& p) : params(p) {}

      //! Default destructor
      ~MesPionRhoxBT2SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<LatticePropagator>& forward_props) const;

    private:
      //! Hide partial constructor
      MesPionRhoxBT2SeqSrc() {}

    private:
      ParamsDir  params;   /*!< Seqsource params */
    };


    //! Construct pion_1-(A1xB_A1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  a1 x B_A1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_5 \gamma_i B_i\f$  
     */
    class MesPionA1xBA1SeqSrc : public HadronSeqSource<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionA1xBA1SeqSrc(const Params& p) : params(p) {}

      //! Default destructor
      ~MesPionA1xBA1SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<LatticePropagator>& forward_props) const;

    private:
      //! Hide partial constructor
      MesPionA1xBA1SeqSrc() {}

    private:
      Params  params;   /*!< Seqsource params */
    };


    //! Construct pion_1-(RhoxB_T1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  a11 x B_T1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_5 \epsilon_{ijk}\gamma_j B_k\f$  
     */
    class MesPionA1xBT1SeqSrc : public HadronSeqSource<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionA1xBT1SeqSrc(const ParamsDir& p) : params(p) {}

      //! Default destructor
      ~MesPionA1xBT1SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<LatticePropagator>& forward_props) const;

    private:
      //! Hide partial constructor
      MesPionA1xBT1SeqSrc() {}

    private:
      ParamsDir  params;   /*!< Seqsource params */
    };


    //! Construct pion_1-(A1xB_T2) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  a1 x B_T2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_5 s_{ijk}\gamma_j B_k\f$  
     */
    class MesPionA1xBT2SeqSrc : public HadronSeqSource<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionA1xBT2SeqSrc(const ParamsDir& p) : params(p) {}

      //! Default destructor
      ~MesPionA1xBT2SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<LatticePropagator>& forward_props) const;

    private:
      //! Hide partial constructor
      MesPionA1xBT2SeqSrc() {}

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
