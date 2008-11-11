// -*- C++ -*-
// $Id: deriv_meson_seqsrc_w.h,v 3.7 2008-11-11 21:27:42 edwards Exp $
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
 *
 *  The sequential sources return
 *    \return \f$\gamma_5 * \Gamma_f^\dag * \gamma_5 * op(F)\$
 *
 */

#ifndef __deriv_meson_seqsrc_w_h__
#define __deriv_meson_seqsrc_w_h__

#include "meas/hadron/meson_seqsrc_w.h"

namespace Chroma 
{

  //! Name and registration
  /*! @ingroup hadron */
  namespace DerivMesonSeqSourceEnv
  {
    bool registerAll();

  
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



    //! Base class for meson deriv-sequential source construction
    /*! @ingroup hadron
     */
    class DerivMesonSeqSourceBase : public MesonSeqSourceBase
    {
    public:
      //! Default destructor
      virtual ~DerivMesonSeqSourceBase() {}

      //! Construct the source
      virtual LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
					   const multi1d<ForwardProp_t>& forward_headers,
					   const multi1d<LatticePropagator>& forward_props) = 0;

      //! Compute the 2-pt at the sink
      virtual Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
				const multi1d<ForwardProp_t>& forward_headers,
				const multi1d<LatticePropagator>& forward_props,
				int gamma_insertion) = 0;

    protected:
      //! Get deriv_length
      virtual int getDerivLength() const = 0;

      //! Apply first deriv (nabla) to the right onto source
      virtual LatticePropagator nabla(const LatticePropagator& forward_prop,
				      const multi1d<LatticeColorMatrix>& u,
				      int mu) const;
      
      //! Apply right "D_i" operator onto source
      virtual LatticePropagator D(const LatticePropagator& forward_prop,
				  const multi1d<LatticeColorMatrix>& u,
				  int mu) const;
      
      //! Apply right "B_i" operator onto source
      virtual LatticePropagator B(const LatticePropagator& forward_prop,
				  const multi1d<LatticeColorMatrix>& u,
				  int mu) const;
      
      //! Apply left and right "nabla_i" onto the source
      /*!
       * \f$\nabla_\mu f(x) = U_\mu(x)f(x+\mu) - U_{-\mu}(x)f(x-\mu)\f$
       *
       * \return $\f \nabla_\mu F(x,0) = U_\mu(x) F(x+\mu)  - U_{x-\mu}^\dag(x-\mu) F(x-\mu)\f$
       */
      virtual LatticePropagator threePtNabla(const LatticePropagator& forward_prop,
					     const multi1d<LatticeColorMatrix>& u,
					     int mu) const;
      
      //! Apply left and right "nabla_i" onto the source
      /*!
       * \f$\nabla_\mu f(x) = U_\mu(x)f(x+\mu) - U_{-\mu}(x)f(x-\mu)\f$
       *
       * \return $\f \nabla_\mu F(x,0) = U_\mu(x) F(x+\mu)  - U_{x-\mu}^\dag(x-\mu) F(x-\mu)\f$
       */
      virtual multi1d<LatticePropagator> threePtNablaVector(const LatticePropagator& forward_prop,
							    const multi1d<LatticeColorMatrix>& u) const;
      
      //! Apply left and right "D_i" operator onto source
      /*!
       * \f$D_i = s_{ijk}\nabla_j\nabla_k\f$
       *
       * where  \f$s_{ijk} = +1 \quad\forall i\ne j, j\ne k, i \ne k\f$
       * 
       * \return $\f D_\mu F(x,0\f$
       */
      virtual LatticePropagator threePtD(const LatticePropagator& forward_prop,
					 const multi1d<LatticeColorMatrix>& u,
					 int mu) const;
      
      //! Apply left and right "D_i" operator onto source
      /*!
       * \f$D_i = s_{ijk}\nabla_j\nabla_k\f$
       *
       * where  \f$s_{ijk} = +1 \quad\forall i\ne j, j\ne k, i \ne k\f$
       * 
       * \return $\f D_\mu F(x,0\f$
       */
      virtual multi1d<LatticePropagator> threePtDVector(const LatticePropagator& forward_prop,
							const multi1d<LatticeColorMatrix>& u) const;
      
      //! Apply left and right "B_i" operator onto source
      /*!
       * \f$B_i = \epsilon_{ijk}\nabla_j\nabla_k\f$
       *
       * \return $\fB_\mu F(z,0)\f$
       */
      virtual LatticePropagator threePtB(const LatticePropagator& forward_prop,
					 const multi1d<LatticeColorMatrix>& u,
					 int mu) const;

      //! Apply left and right "B_i" operator onto source
      /*!
       * \f$B_i = \epsilon_{ijk}\nabla_j\nabla_k\f$
       *
       * \return $\fB_\mu F(z,0)\f$
       */
      virtual multi1d<LatticePropagator> threePtBVector(const LatticePropagator& forward_prop,
							const multi1d<LatticeColorMatrix>& u) const;

      //! Apply left and right "nabla_i" onto the source
      /*!
       * \f$\nabla_\mu f(x) = U_\mu(x)f(x+\mu) - U_{-\mu}(x)f(x-\mu)\f$
       *
       * \return $\f \nabla_\mu F(x,0) = U_\mu(x) F(x+\mu)  - U_{x-\mu}^\dag(x-\mu) F(x-\mu)\f$
       */
      virtual LatticeComplex twoPtNabla(const LatticePropagator& forward_prop,
					const multi1d<LatticeColorMatrix>& u,
					int mu, int g, int gamma_insertion) const;
      
      //! Apply left and right "D_i" operator onto source
      /*!
       * \f$D_i = s_{ijk}\nabla_j\nabla_k\f$
       *
       * where  \f$s_{ijk} = +1 \quad\forall i\ne j, j\ne k, i \ne k\f$
       * 
       * \return $\f D_\mu F(x,0\f$
       */
      virtual LatticeComplex twoPtD(const LatticePropagator& forward_prop,
				    const multi1d<LatticeColorMatrix>& u,
				    int mu, int g, int gamma_insertion) const;
      
      //! Apply left and right "B_i" operator onto source
      /*!
       * \f$B_i = \epsilon_{ijk}\nabla_j\nabla_k\f$
       *
       * \return $\fB_\mu F(z,0)\f$
       */
      virtual LatticeComplex twoPtB(const LatticePropagator& forward_prop,
				    const multi1d<LatticeColorMatrix>& u,
				    int mu, int g, int insertion) const;

      //! Project onto the fixed sink-momentum and return the 2-pt at the sink
      virtual Complex momentumProject(const LatticeComplex& corr_fn) const;
    };


    //--------------------------------------------------------------------------------------------
    //! Base class for meson deriv-sequential source construction
    /*! @ingroup hadron
     */
    class DerivMesonSeqSourceBaseNoDir : public DerivMesonSeqSourceBase
    {
    public:
      //! Default destructor
      DerivMesonSeqSourceBaseNoDir(const Params& p) : params(p) {}

      //! Default destructor
      virtual ~DerivMesonSeqSourceBaseNoDir() {}

      //! Construct the source
      virtual LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
					   const multi1d<ForwardProp_t>& forward_headers,
					   const multi1d<LatticePropagator>& forward_props) = 0;

      //! Compute the 2-pt at the sink
      virtual Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
				const multi1d<ForwardProp_t>& forward_headers,
				const multi1d<LatticePropagator>& forward_props,
				int gamma_insertion) = 0;

    protected:
      //! Set t_srce
      virtual multi1d<int>& getTSrce() {return t_srce;}

      //! Get t_srce
      virtual const multi1d<int>& getTSrce() const {return t_srce;}

      //! Get t_sink
      virtual int getTSink() const {return params.t_sink;}

      //! Get sink_mom
      virtual const multi1d<int>& getSinkMom() const {return params.sink_mom;}
      
      //! Get decay_dir
      virtual const int getDecayDir() const {return params.j_decay;}

      //! Get deriv_length
      virtual int getDerivLength() const {return params.deriv_length;}

    private:
      //! Hide partial constructor
      DerivMesonSeqSourceBaseNoDir() {}

    private:
      multi1d<int>  t_srce;   /*!< Must come from propagator headers */
      Params        params;   /*!< Seqsource params */
    };



    //! Base class for meson deriv-sequential source construction
    /*! @ingroup hadron
     */
    class DerivMesonSeqSourceBaseDir : public DerivMesonSeqSourceBase
    {
    public:
      //! Default destructor
      DerivMesonSeqSourceBaseDir(const ParamsDir& p) : params(p) {}

      //! Default destructor
      virtual ~DerivMesonSeqSourceBaseDir() {}

      //! Construct the source
      virtual LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
					   const multi1d<ForwardProp_t>& forward_headers,
					   const multi1d<LatticePropagator>& forward_props) = 0;

      //! Compute the 2-pt at the sink
      virtual Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
				const multi1d<ForwardProp_t>& forward_headers,
				const multi1d<LatticePropagator>& forward_props,
				int gamma_insertion) = 0;

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

      //! Get deriv_length
      int getDerivLength() const {return params.deriv_length;}

      //! Get deriv_dir
      virtual const int getDerivDir() const {return params.deriv_dir;}

    private:
      //! Hide partial constructor
      DerivMesonSeqSourceBaseDir() {}

    private:
      multi1d<int>  t_srce;   /*!< Must come from propagator headers */
      ParamsDir     params;   /*!< Seqsource params */
    };


    //--------------------------------------------------------------------------------------------
    //! Construct a0-(pionxNabla_T1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  pion x Nabla_T1
     * The sink interpolator structure is
     * \f$\Gamma_f \equiv \gamma_5\nabla_i\f$
     */
    class MesA0PionxNablaT1SeqSrc : public DerivMesonSeqSourceBaseDir
    {
    public:
      //! Full constructor
      MesA0PionxNablaT1SeqSrc(const ParamsDir& p) : DerivMesonSeqSourceBaseDir(p) {}

      //! Default destructor
      ~MesA0PionxNablaT1SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(a0xNabla_T1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  a0 x nabla_T1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \nabla_i\f$
     */
    class MesA0A0xNablaT1SeqSrc : public DerivMesonSeqSourceBaseDir
    {
    public:
      //! Full constructor
      MesA0A0xNablaT1SeqSrc(const ParamsDir& p) : DerivMesonSeqSourceBaseDir(p) {}

      //! Default destructor
      ~MesA0A0xNablaT1SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(a0_2xNabla_T1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  a0_2 x nabla_T1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_4 \nabla_i\f$
     */
    class MesA0A02xNablaT1SeqSrc : public DerivMesonSeqSourceBaseDir
    {
    public:
      //! Full constructor
      MesA0A02xNablaT1SeqSrc(const ParamsDir& p) : DerivMesonSeqSourceBaseDir(p) {}

      //! Default destructor
      ~MesA0A02xNablaT1SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(pion_2xNabla_T1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  pion_2 x nabla_T1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_4 \gamma_5 \nabla_i\f$
     */
    class MesA0Pion2xNablaT1SeqSrc : public DerivMesonSeqSourceBaseDir
    {
    public:
      //! Full constructor
      MesA0Pion2xNablaT1SeqSrc(const ParamsDir& p) : DerivMesonSeqSourceBaseDir(p) {}

      //! Default destructor
      ~MesA0Pion2xNablaT1SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(rhoxNabla_A1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  rho x nabla_A1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_i\nabla_i\f$
     */
    class MesA0RhoxNablaA1SeqSrc : public DerivMesonSeqSourceBaseNoDir
    {
    public:
      //! Full constructor
      MesA0RhoxNablaA1SeqSrc(const Params& p) : DerivMesonSeqSourceBaseNoDir(p) {}

      //! Default destructor
      ~MesA0RhoxNablaA1SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(rhoxNabla_T1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  rho x nabla_T1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \epsilon_{ijk}\gamma_j \nabla_k\f$  
     */
    class MesA0RhoxNablaT1SeqSrc : public DerivMesonSeqSourceBaseDir
    {
    public:
      //! Full constructor
      MesA0RhoxNablaT1SeqSrc(const ParamsDir& p) : DerivMesonSeqSourceBaseDir(p) {}

      //! Default destructor
      ~MesA0RhoxNablaT1SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(rhoxNabla_T2) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  rho x nabla_T2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv s_{ijk}\gamma_j D_k\f$  
     */
    class MesA0RhoxNablaT2SeqSrc : public DerivMesonSeqSourceBaseDir
    {
    public:
      //! Full constructor
      MesA0RhoxNablaT2SeqSrc(const ParamsDir& p) : DerivMesonSeqSourceBaseDir(p) {}

      //! Default destructor
      ~MesA0RhoxNablaT2SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(rhoxNabla_E) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  rho x nabla_E
     * The sink interpolator is   
     * \f$\Gamma_f \equiv S_{\alpha jk}\gamma_j D_k\f$  
     */
    class MesA0RhoxNablaESeqSrc : public DerivMesonSeqSourceBaseDir
    {
    public:
      //! Full constructor
      MesA0RhoxNablaESeqSrc(const ParamsDir& p) : DerivMesonSeqSourceBaseDir(p) {}

      //! Default destructor
      ~MesA0RhoxNablaESeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(rho_2xNabla_A1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  rho_2 x nabla_A1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_4 \gamma_i\nabla_i\f$
     */
    class MesA0Rho2xNablaA1SeqSrc : public DerivMesonSeqSourceBaseNoDir
    {
    public:
      //! Full constructor
      MesA0Rho2xNablaA1SeqSrc(const Params& p) : DerivMesonSeqSourceBaseNoDir(p) {}

      //! Default destructor
      ~MesA0Rho2xNablaA1SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(rho_2xNabla_T1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  rho_2 x nabla_T1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \epsilon_{ijk}\gamma_4\gamma_j \nabla_k\f$  
     */
    class MesA0Rho2xNablaT1SeqSrc : public DerivMesonSeqSourceBaseDir
    {
    public:
      //! Full constructor
      MesA0Rho2xNablaT1SeqSrc(const ParamsDir& p) : DerivMesonSeqSourceBaseDir(p) {}

      //! Default destructor
      ~MesA0Rho2xNablaT1SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(rho_2xNabla_T2) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  rho_2 x nabla_T2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv s_{ijk}\gamma_4\gamma_j D_k\f$  
     */
    class MesA0Rho2xNablaT2SeqSrc : public DerivMesonSeqSourceBaseDir
    {
    public:
      //! Full constructor
      MesA0Rho2xNablaT2SeqSrc(const ParamsDir& p) : DerivMesonSeqSourceBaseDir(p) {}

      //! Default destructor
      ~MesA0Rho2xNablaT2SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(rho_2xNabla_E) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  rho_2 x nabla_E
     * The sink interpolator is   
     * \f$\Gamma_f \equiv S_{\alpha jk}\gamma_4\gamma_j D_k\f$  
     */
    class MesA0Rho2xNablaESeqSrc : public DerivMesonSeqSourceBaseDir
    {
    public:
      //! Full constructor
      MesA0Rho2xNablaESeqSrc(const ParamsDir& p) : DerivMesonSeqSourceBaseDir(p) {}

      //! Default destructor
      ~MesA0Rho2xNablaESeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(a1xNabla_A1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  a1 x D_A1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_5\gamma_i \nabla_i\f$  
     */
    class MesA0A1xNablaA1SeqSrc : public DerivMesonSeqSourceBaseNoDir
    {
    public:
      //! Full constructor
      MesA0A1xNablaA1SeqSrc(const Params& p) : DerivMesonSeqSourceBaseNoDir(p) {}

      //! Default destructor
      ~MesA0A1xNablaA1SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(a1xNabla_T1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  a1 x nabla_T1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_5 \epsilon_{ijk}\gamma_j \nabla_k\f$  
     */
    class MesA0A1xNablaT1SeqSrc : public DerivMesonSeqSourceBaseDir
    {
    public:
      //! Full constructor
      MesA0A1xNablaT1SeqSrc(const ParamsDir& p) : DerivMesonSeqSourceBaseDir(p) {}

      //! Default destructor
      ~MesA0A1xNablaT1SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(a1xNabla_T2) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  a1 x nabla_T2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_5 s_{ijk}\gamma_j \nabla_k\f$  
     */
    class MesA0A1xNablaT2SeqSrc : public DerivMesonSeqSourceBaseDir
    {
    public:
      //! Full constructor
      MesA0A1xNablaT2SeqSrc(const ParamsDir& p) : DerivMesonSeqSourceBaseDir(p) {}

      //! Default destructor
      ~MesA0A1xNablaT2SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(a1xNabla_E) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  a1 x nabla_E
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_5 S_{\alpha jk}\gamma_j \nabla_k\f$  
     */
    class MesA0A1xNablaESeqSrc : public DerivMesonSeqSourceBaseDir
    {
    public:
      //! Full constructor
      MesA0A1xNablaESeqSrc(const ParamsDir& p) : DerivMesonSeqSourceBaseDir(p) {}

      //! Default destructor
      ~MesA0A1xNablaESeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(b1xNabla_A1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  b1 x nabla_A1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_4\gamma_5 \gamma_i \nabla_i\f$  
     */
    class MesA0B1xNablaA1SeqSrc : public DerivMesonSeqSourceBaseNoDir
    {
    public:
      //! Full constructor
      MesA0B1xNablaA1SeqSrc(const Params& p) : DerivMesonSeqSourceBaseNoDir(p) {}

      //! Default destructor
      ~MesA0B1xNablaA1SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(b1xNabla_T1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  b1 x nabla_T1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_4\gamma_5\epsilon_{ijk}\gamma_j \nabla_k\f$  
     */
    class MesA0B1xNablaT1SeqSrc : public DerivMesonSeqSourceBaseDir
    {
    public:
      //! Full constructor
      MesA0B1xNablaT1SeqSrc(const ParamsDir& p) : DerivMesonSeqSourceBaseDir(p) {}

      //! Default destructor
      ~MesA0B1xNablaT1SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(b1xNabla_T2) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  b1 x nabla_T2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_4\gamma_5 s_{ijk}\gamma_j \nabla_k\f$  
     */
    class MesA0B1xNablaT2SeqSrc : public DerivMesonSeqSourceBaseDir
    {
    public:
      //! Full constructor
      MesA0B1xNablaT2SeqSrc(const ParamsDir& p) : DerivMesonSeqSourceBaseDir(p) {}

      //! Default destructor
      ~MesA0B1xNablaT2SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(b1xNabla_E) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  b1 x nabla_E
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_4\gamma_5 S_{\alpha jk}\gamma_j \nabla_k\f$  
     */
    class MesA0B1xNablaESeqSrc : public DerivMesonSeqSourceBaseDir
    {
    public:
      //! Full constructor
      MesA0B1xNablaESeqSrc(const ParamsDir& p) : DerivMesonSeqSourceBaseDir(p) {}

      //! Default destructor
      ~MesA0B1xNablaESeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(pionxD_T2) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  pion x D_T2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_5 D_i\f$  
     */
    class MesA0PionxDT2SeqSrc : public DerivMesonSeqSourceBaseDir
    {
    public:
      //! Full constructor
      MesA0PionxDT2SeqSrc(const ParamsDir& p) : DerivMesonSeqSourceBaseDir(p) {}

      //! Default destructor
      ~MesA0PionxDT2SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(a0xD_T2) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  a0 x D_T2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv D_i\f$  
     */
    class MesA0A0xDT2SeqSrc : public DerivMesonSeqSourceBaseDir
    {
    public:
      //! Full constructor
      MesA0A0xDT2SeqSrc(const ParamsDir& p) : DerivMesonSeqSourceBaseDir(p) {}

      //! Default destructor
      ~MesA0A0xDT2SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(a0_2xD_T2) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  a0_2 x D_T2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_4 D_i\f$  
     */
    class MesA0A02xDT2SeqSrc : public DerivMesonSeqSourceBaseDir
    {
    public:
      //! Full constructor
      MesA0A02xDT2SeqSrc(const ParamsDir& p) : DerivMesonSeqSourceBaseDir(p) {}

      //! Default destructor
      ~MesA0A02xDT2SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(pion_2xD_T2) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  pion_2 x D_T2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_4\gamma_5 D_i\f$  
     */
    class MesA0Pion2xDT2SeqSrc : public DerivMesonSeqSourceBaseDir
    {
    public:
      //! Full constructor
      MesA0Pion2xDT2SeqSrc(const ParamsDir& p) : DerivMesonSeqSourceBaseDir(p) {}

      //! Default destructor
      ~MesA0Pion2xDT2SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(a1xD_A2) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  a1 x D_A2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_5\gamma_i D_i\f$  
     */
    class MesA0A1xDA2SeqSrc : public DerivMesonSeqSourceBaseNoDir
    {
    public:
      //! Full constructor
      MesA0A1xDA2SeqSrc(const Params& p) : DerivMesonSeqSourceBaseNoDir(p) {}

      //! Default destructor
      ~MesA0A1xDA2SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(a1xD_T1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  a1 x D_T1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_5 s_{ijk}\gamma_j D_k\f$  
     */
    class MesA0A1xDT1SeqSrc : public DerivMesonSeqSourceBaseDir
    {
    public:
      //! Full constructor
      MesA0A1xDT1SeqSrc(const ParamsDir& p) : DerivMesonSeqSourceBaseDir(p) {}

      //! Default destructor
      ~MesA0A1xDT1SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(a1xD_T2) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  a1 x D_T2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_5\epsilon_{ijk}\gamma_j D_k\f$  
     */
    class MesA0A1xDT2SeqSrc : public DerivMesonSeqSourceBaseDir
    {
    public:
      //! Full constructor
      MesA0A1xDT2SeqSrc(const ParamsDir& p) : DerivMesonSeqSourceBaseDir(p) {}

      //! Default destructor
      ~MesA0A1xDT2SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(a1xD_E) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  a1 x D_E
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_5 S_{\alpha jk}\gamma_j D_k\f$  
     */
    class MesA0A1xDESeqSrc : public DerivMesonSeqSourceBaseDir
    {
    public:
      //! Full constructor
      MesA0A1xDESeqSrc(const ParamsDir& p) : DerivMesonSeqSourceBaseDir(p) {}

      //! Default destructor
      ~MesA0A1xDESeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(b1xD_A2) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  b1 x D_A2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_4\gamma_5 \gamma_i D_i\f$  
     */
    class MesA0B1xDA2SeqSrc : public DerivMesonSeqSourceBaseNoDir
    {
    public:
      //! Full constructor
      MesA0B1xDA2SeqSrc(const Params& p) : DerivMesonSeqSourceBaseNoDir(p) {}

      //! Default destructor
      ~MesA0B1xDA2SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(b1xD_T1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  b1 x D_T1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_4\gamma_5 s_{ijk}\gamma_j D_k\f$  
     */
    class MesA0B1xDT1SeqSrc : public DerivMesonSeqSourceBaseDir
    {
    public:
      //! Full constructor
      MesA0B1xDT1SeqSrc(const ParamsDir& p) : DerivMesonSeqSourceBaseDir(p) {}

      //! Default destructor
      ~MesA0B1xDT1SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(b1xD_T2) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  b1 x D_T2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_4\gamma_5 \epsilon_{ijk}\gamma_j D_k\f$  
     */
    class MesA0B1xDT2SeqSrc : public DerivMesonSeqSourceBaseDir
    {
    public:
      //! Full constructor
      MesA0B1xDT2SeqSrc(const ParamsDir& p) : DerivMesonSeqSourceBaseDir(p) {}

      //! Default destructor
      ~MesA0B1xDT2SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(b1xD_E) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  b1 x D_E
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_4\gamma_5 S_{\alpha jk}\gamma_j D_k\f$  
     */
    class MesA0B1xDESeqSrc : public DerivMesonSeqSourceBaseDir
    {
    public:
      //! Full constructor
      MesA0B1xDESeqSrc(const ParamsDir& p) : DerivMesonSeqSourceBaseDir(p) {}

      //! Default destructor
      ~MesA0B1xDESeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(rhoxD_A2) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  rho x D_A2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_i D_i\f$  
     */
    class MesA0RhoxDA2SeqSrc : public DerivMesonSeqSourceBaseNoDir
    {
    public:
      //! Full constructor
      MesA0RhoxDA2SeqSrc(const Params& p) : DerivMesonSeqSourceBaseNoDir(p) {}

      //! Default destructor
      ~MesA0RhoxDA2SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(rhoxD_T1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  rho x D_T1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv s_{ijk}\gamma_j D_k\f$  
     */
    class MesA0RhoxDT1SeqSrc : public DerivMesonSeqSourceBaseDir
    {
    public:
      //! Full constructor
      MesA0RhoxDT1SeqSrc(const ParamsDir& p) : DerivMesonSeqSourceBaseDir(p) {}

      //! Default destructor
      ~MesA0RhoxDT1SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(rhoxD_T2) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  rho x D_T2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \epsilon_{ijk}\gamma_j D_k\f$  
     */
    class MesA0RhoxDT2SeqSrc : public DerivMesonSeqSourceBaseDir
    {
    public:
      //! Full constructor
      MesA0RhoxDT2SeqSrc(const ParamsDir& p) : DerivMesonSeqSourceBaseDir(p) {}

      //! Default destructor
      ~MesA0RhoxDT2SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(rhoxD_E) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  rho x D_E
     * The sink interpolator is   
     * \f$\Gamma_f \equiv S_{\alpha jk}\gamma_j D_k\f$  
     */
    class MesA0RhoxDESeqSrc : public DerivMesonSeqSourceBaseDir
    {
    public:
      //! Full constructor
      MesA0RhoxDESeqSrc(const ParamsDir& p) : DerivMesonSeqSourceBaseDir(p) {}

      //! Default destructor
      ~MesA0RhoxDESeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(rho_2xD_A2) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  rho_2 x D_A2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_4 \gamma_i D_i\f$  
     */
    class MesA0Rho2xDA2SeqSrc : public DerivMesonSeqSourceBaseNoDir
    {
    public:
      //! Full constructor
      MesA0Rho2xDA2SeqSrc(const Params& p) : DerivMesonSeqSourceBaseNoDir(p) {}

      //! Default destructor
      ~MesA0Rho2xDA2SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(rho_2xD_T1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  rho_2 x D_T1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv s_{ijk} \gamma_4 \gamma_j D_k\f$  
     */
    class MesA0Rho2xDT1SeqSrc : public DerivMesonSeqSourceBaseDir
    {
    public:
      //! Full constructor
      MesA0Rho2xDT1SeqSrc(const ParamsDir& p) : DerivMesonSeqSourceBaseDir(p) {}

      //! Default destructor
      ~MesA0Rho2xDT1SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(rho_2xD_T2) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  rho_2 x D_T2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \epsilon_{ijk} \gamma_4 \gamma_j D_k\f$  
     */
    class MesA0Rho2xDT2SeqSrc : public DerivMesonSeqSourceBaseDir
    {
    public:
      //! Full constructor
      MesA0Rho2xDT2SeqSrc(const ParamsDir& p) : DerivMesonSeqSourceBaseDir(p) {}

      //! Default destructor
      ~MesA0Rho2xDT2SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(rho_2xD_E) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  rho_2 x D_E
     * The sink interpolator is   
     * \f$\Gamma_f \equiv S_{\alpha jk} \gamma_4 \gamma_j D_k\f$  
     */
    class MesA0Rho2xDESeqSrc : public DerivMesonSeqSourceBaseDir
    {
    public:
      //! Full constructor
      MesA0Rho2xDESeqSrc(const ParamsDir& p) : DerivMesonSeqSourceBaseDir(p) {}

      //! Default destructor
      ~MesA0Rho2xDESeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(pionxB_T1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  pion x B_T1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_5 B_i\f$  
     */
    class MesA0PionxBT1SeqSrc : public DerivMesonSeqSourceBaseDir
    {
    public:
      //! Full constructor
      MesA0PionxBT1SeqSrc(const ParamsDir& p) : DerivMesonSeqSourceBaseDir(p) {}

      //! Default destructor
      ~MesA0PionxBT1SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(a0xB_T1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  a0 x B_T1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv B_i\f$  
     */
    class MesA0A0xBT1SeqSrc : public DerivMesonSeqSourceBaseDir
    {
    public:
      //! Full constructor
      MesA0A0xBT1SeqSrc(const ParamsDir& p) : DerivMesonSeqSourceBaseDir(p) {}

      //! Default destructor
      ~MesA0A0xBT1SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(a0_2xB_T1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  a0_2 x B_T1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_4 B_i\f$  
     */
    class MesA0A02xBT1SeqSrc : public DerivMesonSeqSourceBaseDir
    {
    public:
      //! Full constructor
      MesA0A02xBT1SeqSrc(const ParamsDir& p) : DerivMesonSeqSourceBaseDir(p) {}

      //! Default destructor
      ~MesA0A02xBT1SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(pion_2xB_T1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  pion_2 x B_T1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_4\gamma_5 B_i\f$  
     */
    class MesA0Pion2xBT1SeqSrc : public DerivMesonSeqSourceBaseDir
    {
    public:
      //! Full constructor
      MesA0Pion2xBT1SeqSrc(const ParamsDir& p) : DerivMesonSeqSourceBaseDir(p) {}

      //! Default destructor
      ~MesA0Pion2xBT1SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(rhoxB_A1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is rho x B_A1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_i B_i\f$  
     */
    class MesA0RhoxBA1SeqSrc : public DerivMesonSeqSourceBaseNoDir
    {
    public:
      //! Full constructor
      MesA0RhoxBA1SeqSrc(const Params& p) : DerivMesonSeqSourceBaseNoDir(p) {}

      //! Default destructor
      ~MesA0RhoxBA1SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(rhoxB_T1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is rho x B_T1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \epsilon_{ijk}\gamma_j B_k\f$  
     */
    class MesA0RhoxBT1SeqSrc : public DerivMesonSeqSourceBaseDir
    {
    public:
      //! Full constructor
      MesA0RhoxBT1SeqSrc(const ParamsDir& p) : DerivMesonSeqSourceBaseDir(p) {}

      //! Default destructor
      ~MesA0RhoxBT1SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(rhoxB_T2) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  rho x B_T2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv s_{ijk}\gamma_j D_k\f$  
     */
    class MesA0RhoxBT2SeqSrc : public DerivMesonSeqSourceBaseDir
    {
    public:
      //! Full constructor
      MesA0RhoxBT2SeqSrc(const ParamsDir& p) : DerivMesonSeqSourceBaseDir(p) {}

      //! Default destructor
      ~MesA0RhoxBT2SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(rhoxB_E) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  rho x B_E
     * The sink interpolator is   
     * \f$\Gamma_f \equiv S_{\alpha jk}\gamma_j D_k\f$  
     */
    class MesA0RhoxBESeqSrc : public DerivMesonSeqSourceBaseDir
    {
    public:
      //! Full constructor
      MesA0RhoxBESeqSrc(const ParamsDir& p) : DerivMesonSeqSourceBaseDir(p) {}

      //! Default destructor
      ~MesA0RhoxBESeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(rho_2xB_A1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is rho_2 x B_A1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_4 \gamma_i B_i\f$  
     */
    class MesA0Rho2xBA1SeqSrc : public DerivMesonSeqSourceBaseNoDir
    {
    public:
      //! Full constructor
      MesA0Rho2xBA1SeqSrc(const Params& p) : DerivMesonSeqSourceBaseNoDir(p) {}

      //! Default destructor
      ~MesA0Rho2xBA1SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(rho_2xB_T1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is rho_2 x B_T1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \epsilon_{ijk} \gamma_4 \gamma_j B_k\f$  
     */
    class MesA0Rho2xBT1SeqSrc : public DerivMesonSeqSourceBaseDir
    {
    public:
      //! Full constructor
      MesA0Rho2xBT1SeqSrc(const ParamsDir& p) : DerivMesonSeqSourceBaseDir(p) {}

      //! Default destructor
      ~MesA0Rho2xBT1SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(rho_2xB_T2) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  rho_2 x B_T2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv s_{ijk} \gamma_4 \gamma_j D_k\f$  
     */
    class MesA0Rho2xBT2SeqSrc : public DerivMesonSeqSourceBaseDir
    {
    public:
      //! Full constructor
      MesA0Rho2xBT2SeqSrc(const ParamsDir& p) : DerivMesonSeqSourceBaseDir(p) {}

      //! Default destructor
      ~MesA0Rho2xBT2SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(rho_2xB_E) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  rho_2 x B_E
     * The sink interpolator is   
     * \f$\Gamma_f \equiv S_{\alpha jk} \gamma_4 \gamma_j D_k\f$  
     */
    class MesA0Rho2xBESeqSrc : public DerivMesonSeqSourceBaseDir
    {
    public:
      //! Full constructor
      MesA0Rho2xBESeqSrc(const ParamsDir& p) : DerivMesonSeqSourceBaseDir(p) {}

      //! Default destructor
      ~MesA0Rho2xBESeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(a1xB_A1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  a1 x B_A1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_5 \gamma_i B_i\f$  
     */
    class MesA0A1xBA1SeqSrc : public DerivMesonSeqSourceBaseNoDir
    {
    public:
      //! Full constructor
      MesA0A1xBA1SeqSrc(const Params& p) : DerivMesonSeqSourceBaseNoDir(p) {}

      //! Default destructor
      ~MesA0A1xBA1SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(rhoxB_T1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  a11 x B_T1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_5 \epsilon_{ijk}\gamma_j B_k\f$  
     */
    class MesA0A1xBT1SeqSrc : public DerivMesonSeqSourceBaseDir
    {
    public:
      //! Full constructor
      MesA0A1xBT1SeqSrc(const ParamsDir& p) : DerivMesonSeqSourceBaseDir(p) {}

      //! Default destructor
      ~MesA0A1xBT1SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(a1xB_T2) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  a1 x B_T2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_5 s_{ijk}\gamma_j B_k\f$  
     */
    class MesA0A1xBT2SeqSrc : public DerivMesonSeqSourceBaseDir
    {
    public:
      //! Full constructor
      MesA0A1xBT2SeqSrc(const ParamsDir& p) : DerivMesonSeqSourceBaseDir(p) {}

      //! Default destructor
      ~MesA0A1xBT2SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(a1xB_E) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  a1 x B_E
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_5 S_{\alpha jk}\gamma_j B_k\f$  
     */
    class MesA0A1xBESeqSrc : public DerivMesonSeqSourceBaseDir
    {
    public:
      //! Full constructor
      MesA0A1xBESeqSrc(const ParamsDir& p) : DerivMesonSeqSourceBaseDir(p) {}

      //! Default destructor
      ~MesA0A1xBESeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(b1xB_A1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  b1 x B_A1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_4 \gamma_5 \gamma_i B_i\f$  
     */
    class MesA0B1xBA1SeqSrc : public DerivMesonSeqSourceBaseNoDir
    {
    public:
      //! Full constructor
      MesA0B1xBA1SeqSrc(const Params& p) : DerivMesonSeqSourceBaseNoDir(p) {}

      //! Default destructor
      ~MesA0B1xBA1SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(b1xB_T1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  b1 x B_T1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_4 \gamma_5 \epsilon_{ijk}\gamma_j B_k\f$  
     */
    class MesA0B1xBT1SeqSrc : public DerivMesonSeqSourceBaseDir
    {
    public:
      //! Full constructor
      MesA0B1xBT1SeqSrc(const ParamsDir& p) : DerivMesonSeqSourceBaseDir(p) {}

      //! Default destructor
      ~MesA0B1xBT1SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(b1xB_T2) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  b1 x B_T2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_4 \gamma_5 s_{ijk}\gamma_j B_k\f$  
     */
    class MesA0B1xBT2SeqSrc : public DerivMesonSeqSourceBaseDir
    {
    public:
      //! Full constructor
      MesA0B1xBT2SeqSrc(const ParamsDir& p) : DerivMesonSeqSourceBaseDir(p) {}

      //! Default destructor
      ~MesA0B1xBT2SeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
    };


    //! Construct a0-(b1xB_E) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  b1 x B_E
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_4 \gamma_5 S_{\alpha jk}\gamma_j B_k\f$  
     */
    class MesA0B1xBESeqSrc : public DerivMesonSeqSourceBaseDir
    {
    public:
      //! Full constructor
      MesA0B1xBESeqSrc(const ParamsDir& p) : DerivMesonSeqSourceBaseDir(p) {}

      //! Default destructor
      ~MesA0B1xBESeqSrc() {}
      
      //! Construct sequential source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);
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
