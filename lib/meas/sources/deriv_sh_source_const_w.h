// -*- C++ -*-
// $Id: deriv_sh_source_const_w.h,v 1.2 2006-02-20 17:50:13 edwards Exp $
/*! \file
 *  \brief Construct derivative source construction
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

#ifndef __deriv_sh_source_const_w_h__
#define __deriv_sh_source_const_w_h__

#include "meas/sources/source_construction.h"

namespace Chroma 
{

  //! Name and registration
  /*! @ingroup hadron */
  namespace DerivShellSourceConstEnv
  {
    extern const bool registered;

  
    //! Deriv meson sequential source parameters
    /*! @ingroup hadron */
    struct Params
    {
      Params();
      Params(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;

      std::string      source_type;          /*!< source smearing type */

      std::string      quark_smearing;       /*!< xml string holding smearing params */
      std::string      quark_smearing_type;  /*!< quark smearing type name */

      int              j_decay;              /*!< Decay direction */
      int              deriv_length;         /*!< Displacement length in derivative */
      multi1d<int>     t_srce;               /*!< source location */

      std::string      link_smearing;        /*!< link smearing xml */
      std::string      link_smearing_type;   /*!< link smearing type name */
    };

    //! Deriv meson sequential source parameters
    /*! @ingroup hadron */
    struct ParamsDir
    {
      ParamsDir();
      ParamsDir(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;

      //! A type of conversion constructor
      ParamsDir(const Params& p) : source_type(p.source_type), quark_smearing(p.quark_smearing),
				   quark_smearing_type(p.quark_smearing_type), j_decay(p.j_decay),
				   deriv_dir(-1), deriv_length(p.deriv_length), t_srce(p.t_srce),
				   link_smearing(p.link_smearing), 
				   link_smearing_type(p.link_smearing_type)
	{}


      std::string      source_type;          /*!< source smearing type */

      std::string      quark_smearing;       /*!< xml string holding smearing params */
      std::string      quark_smearing_type;  /*!< quark smearing type name */

      int              j_decay;              /*!< Decay direction */
      int              deriv_dir;            /*!< Polarization direction */
      int              deriv_length;         /*!< Displacement length in derivative */
      multi1d<int>     t_srce;               /*!< source location */

      std::string      link_smearing;        /*!< link smearing xml */
      std::string      link_smearing_type;   /*!< link smearing type name */
    };


    //! Derivative base class source construction
    /*! @ingroup sources
     *
     * Abstract base class for derivative source construction
     */
    template<typename T>
    class DerivSourceConst : public QuarkSourceConstruction<T>
    {
    public:
      //! Construct the source
      T operator()(const multi1d<LatticeColorMatrix>& u) const;

    protected:
      //! Get the parameters
      virtual const ParamsDir& getParamsDir() const = 0;

      //! Do the spin actions on the initial point source
      virtual T deriv(const multi1d<LatticeColorMatrix>& u, const T& pt) const = 0;
    };



    //! Construct pion_1-(PionxNabla_T1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  Pion x Nabla_T1
     * The sink interpolator structure is
     * \f$\Gamma_f \equiv \gamma_5\nabla_i\f$
     */
    class MesPionPionxNablaT1SrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionPionxNablaT1SrcConst(const ParamsDir& p) : params(p) {}

    protected:
      //! Get the parameters
      const ParamsDir& getParamsDir() const {return params;}

      //! Do the spin actions on the initial point source
      LatticePropagator deriv(const multi1d<LatticeColorMatrix>& u, const LatticePropagator& pt) const;
      
    private:
      ParamsDir  params;   /*!< source params */
    };


    //! Construct pion_1-(A0xNabla_T1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  a0 x nabla_T1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \nabla_i\f$
     */
    class MesPionA0xNablaT1SrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionA0xNablaT1SrcConst(const ParamsDir& p) : params(p) {}

    protected:
      //! Get the parameters
      const ParamsDir& getParamsDir() const {return params;}

      //! Do the spin actions on the initial point source
      LatticePropagator deriv(const multi1d<LatticeColorMatrix>& u, const LatticePropagator& pt) const;

    private:
      ParamsDir  params;   /*!< source params */
    };


    //! Construct pion_1-(A0_2xNabla_T1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  a0_2 x nabla_T1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_4 \nabla_i\f$
     */
    class MesPionA02xNablaT1SrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionA02xNablaT1SrcConst(const ParamsDir& p) : params(p) {}

    protected:
      //! Get the parameters
      const ParamsDir& getParamsDir() const {return params;}

      //! Do the spin actions on the initial point source
      LatticePropagator deriv(const multi1d<LatticeColorMatrix>& u, const LatticePropagator& pt) const;
      
    private:
      ParamsDir  params;   /*!< source params */
    };


    //! Construct pion_1-(RhoxNabla_A1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  rho x nabla_A1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_i\nabla_i\f$
     */
    class MesPionRhoxNablaA1SrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionRhoxNablaA1SrcConst(const Params& p) : params(p) {}

    protected:
      //! Get the parameters
      const ParamsDir& getParamsDir() const {return params;}

      //! Do the spin actions on the initial point source
      LatticePropagator deriv(const multi1d<LatticeColorMatrix>& u, const LatticePropagator& pt) const;
      
    private:
      ParamsDir  params;   /*!< source params */
    };


    //! Construct pion_1-(RhoxNabla_T1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  rho x nabla_T1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \epsilon_{ijk}\gamma_j \nabla_k\f$  
     */
    class MesPionRhoxNablaT1SrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionRhoxNablaT1SrcConst(const ParamsDir& p) : params(p) {}

    protected:
      //! Get the parameters
      const ParamsDir& getParamsDir() const {return params;}

      //! Do the spin actions on the initial point source
      LatticePropagator deriv(const multi1d<LatticeColorMatrix>& u, const LatticePropagator& pt) const;

    private:
      ParamsDir  params;   /*!< source params */
    };


    //! Construct pion_1-(RhoxNabla_T2) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  rho x nabla_T2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv s_{ijk}\gamma_j D_k\f$  
     */
    class MesPionRhoxNablaT2SrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionRhoxNablaT2SrcConst(const ParamsDir& p) : params(p) {}

    protected:
      //! Get the parameters
      const ParamsDir& getParamsDir() const {return params;}

      //! Do the spin actions on the initial point source
      LatticePropagator deriv(const multi1d<LatticeColorMatrix>& u, const LatticePropagator& pt) const;

    private:
      ParamsDir  params;   /*!< source params */
    };


    //! Construct pion_1-(A1xNabla_A1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  a1 x D_A1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_5\gamma_i \nabla_i\f$  
     */
    class MesPionA1xNablaA1SrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionA1xNablaA1SrcConst(const Params& p) : params(p) {}

    protected:
      //! Get the parameters
      const ParamsDir& getParamsDir() const {return params;}

      //! Do the spin actions on the initial point source
      LatticePropagator deriv(const multi1d<LatticeColorMatrix>& u, const LatticePropagator& pt) const;

    private:
      ParamsDir  params;   /*!< source params */
    };


    //! Construct pion_1-(A1xNabla_T2) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  a1 x nabla_T2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_5 s_{ijk}\gamma_j \nabla_k\f$  
     */
    class MesPionA1xNablaT2SrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionA1xNablaT2SrcConst(const ParamsDir& p) : params(p) {}

    protected:
      //! Get the parameters
      const ParamsDir& getParamsDir() const {return params;}

      //! Do the spin actions on the initial point source
      LatticePropagator deriv(const multi1d<LatticeColorMatrix>& u, const LatticePropagator& pt) const;

    private:
      ParamsDir  params;   /*!< source params */
    };


    //! Construct pion_1-(A1xNabla_E) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  a1 x nabla_E
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_5 S_{\alpha jk}\gamma_j \nabla_k\f$  
     */
    class MesPionA1xNablaESrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionA1xNablaESrcConst(const ParamsDir& p) : params(p) {}

    protected:
      //! Get the parameters
      const ParamsDir& getParamsDir() const {return params;}

      //! Do the spin actions on the initial point source
      LatticePropagator deriv(const multi1d<LatticeColorMatrix>& u, const LatticePropagator& pt) const;

    private:
      ParamsDir  params;   /*!< source params */
    };


    //! Construct pion_1-(B1xNabla_T1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  b1 x nabla_T1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_4\gamma_5\epsilon_{ijk}\gamma_j \nabla_k\f$  
     */
    class MesPionB1xNablaT1SrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionB1xNablaT1SrcConst(const ParamsDir& p) : params(p) {}

    protected:
      //! Get the parameters
      const ParamsDir& getParamsDir() const {return params;}

      //! Do the spin actions on the initial point source
      LatticePropagator deriv(const multi1d<LatticeColorMatrix>& u, const LatticePropagator& pt) const;

    private:
      ParamsDir  params;   /*!< source params */
    };


    //! Construct pion_1-(A0_2xD_T2) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  a0_2 x D_T2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_4 D_i\f$  
     */
    class MesPionA02xDT2SrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionA02xDT2SrcConst(const ParamsDir& p) : params(p) {}

    protected:
      //! Get the parameters
      const ParamsDir& getParamsDir() const {return params;}

      //! Do the spin actions on the initial point source
      LatticePropagator deriv(const multi1d<LatticeColorMatrix>& u, const LatticePropagator& pt) const;

    private:
      ParamsDir  params;   /*!< source params */
    };


    //! Construct pion_1-(A1xD_A2) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  a1 x D_A2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_5\gamma_i D_i\f$  
     */
    class MesPionA1xDA2SrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionA1xDA2SrcConst(const Params& p) : params(p) {}

    protected:
      //! Get the parameters
      const ParamsDir& getParamsDir() const {return params;}

      //! Do the spin actions on the initial point source
      LatticePropagator deriv(const multi1d<LatticeColorMatrix>& u, const LatticePropagator& pt) const;

    private:
      ParamsDir  params;   /*!< source params */
    };


    //! Construct pion_1-(A1xD_E) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  a1 x D_E
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_5 S_{\alpha jk}\gamma_j D_k\f$  
     */
    class MesPionA1xDESrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionA1xDESrcConst(const ParamsDir& p) : params(p) {}

    protected:
      //! Get the parameters
      const ParamsDir& getParamsDir() const {return params;}

      //! Do the spin actions on the initial point source
      LatticePropagator deriv(const multi1d<LatticeColorMatrix>& u, const LatticePropagator& pt) const;

    private:
      ParamsDir  params;   /*!< source params */
    };


    //! Construct pion_1-(A1xD_T1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  a1 x D_T1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_5 s_{ijk}\gamma_j D_k\f$  
     */
    class MesPionA1xDT1SrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionA1xDT1SrcConst(const ParamsDir& p) : params(p) {}

    protected:
      //! Get the parameters
      const ParamsDir& getParamsDir() const {return params;}

      //! Do the spin actions on the initial point source
      LatticePropagator deriv(const multi1d<LatticeColorMatrix>& u, const LatticePropagator& pt) const;

    private:
      ParamsDir  params;   /*!< source params */
    };


    //! Construct pion_1-(A1xD_T2) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  a1 x D_T2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_5\epsilon_{ijk}\gamma_j D_k\f$  
     */
    class MesPionA1xDT2SrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionA1xDT2SrcConst(const ParamsDir& p) : params(p) {}

    protected:
      //! Get the parameters
      const ParamsDir& getParamsDir() const {return params;}

      //! Do the spin actions on the initial point source
      LatticePropagator deriv(const multi1d<LatticeColorMatrix>& u, const LatticePropagator& pt) const;

    private:
      ParamsDir  params;   /*!< source params */
    };


    //! Construct pion_1-(B1xD_A2) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  b1 x D_A2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_4\gamma_5 \gamma_i D_i\f$  
     */
    class MesPionB1xDA2SrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionB1xDA2SrcConst(const Params& p) : params(p) {}

    protected:
      //! Get the parameters
      const ParamsDir& getParamsDir() const {return params;}

      //! Do the spin actions on the initial point source
      LatticePropagator deriv(const multi1d<LatticeColorMatrix>& u, const LatticePropagator& pt) const;

    private:
      ParamsDir  params;   /*!< source params */
    };

    //! Construct pion_1-(B1xD_E) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  b1 x D_E
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_4\gamma_5 S_{\alpha jk}\gamma_j D_k\f$  
     */
    class MesPionB1xDESrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionB1xDESrcConst(const ParamsDir& p) : params(p) {}

    protected:
      //! Get the parameters
      const ParamsDir& getParamsDir() const {return params;}

      //! Do the spin actions on the initial point source
      LatticePropagator deriv(const multi1d<LatticeColorMatrix>& u, const LatticePropagator& pt) const;

    private:
      ParamsDir  params;   /*!< source params */
    };


    //! Construct pion_1-(B1xD_T1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  b1 x D_T1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_4\gamma_5 s_{ijk}\gamma_j D_k\f$  
     */
    class MesPionB1xDT1SrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionB1xDT1SrcConst(const ParamsDir& p) : params(p) {}

    protected:
      //! Get the parameters
      const ParamsDir& getParamsDir() const {return params;}

      //! Do the spin actions on the initial point source
      LatticePropagator deriv(const multi1d<LatticeColorMatrix>& u, const LatticePropagator& pt) const;

    private:
      ParamsDir  params;   /*!< source params */
    };


    //! Construct pion_1-(B1xD_T2) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  b1 x D_T2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_4\gamma_5 \epsilon_{ijk}\gamma_j D_k\f$  
     */
    class MesPionB1xDT2SrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionB1xDT2SrcConst(const ParamsDir& p) : params(p) {}

    protected:
      //! Get the parameters
      const ParamsDir& getParamsDir() const {return params;}

      //! Do the spin actions on the initial point source
      LatticePropagator deriv(const multi1d<LatticeColorMatrix>& u, const LatticePropagator& pt) const;

    private:
      ParamsDir  params;   /*!< source params */
    };


    //! Construct pion_1-(RhoxD_A2) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  rho x D_A2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_i D_i\f$  
     */
    class MesPionRhoxDA2SrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionRhoxDA2SrcConst(const Params& p) : params(p) {}

    protected:
      //! Get the parameters
      const ParamsDir& getParamsDir() const {return params;}

      //! Do the spin actions on the initial point source
      LatticePropagator deriv(const multi1d<LatticeColorMatrix>& u, const LatticePropagator& pt) const;

    private:
      ParamsDir  params;   /*!< source params */
    };


    //! Construct pion_1-(RhoxD_T1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  rho x D_T1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv s_{ijk}\gamma_j D_k\f$  
     */
    class MesPionRhoxDT1SrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionRhoxDT1SrcConst(const ParamsDir& p) : params(p) {}

    protected:
      //! Get the parameters
      const ParamsDir& getParamsDir() const {return params;}

      //! Do the spin actions on the initial point source
      LatticePropagator deriv(const multi1d<LatticeColorMatrix>& u, const LatticePropagator& pt) const;

    private:
      ParamsDir  params;   /*!< source params */
    };


    //! Construct pion_1-(RhoxD_T2) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  rho x D_T2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \epsilon_{ijk}\gamma_j D_k\f$  
     */
    class MesPionRhoxDT2SrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionRhoxDT2SrcConst(const ParamsDir& p) : params(p) {}

    protected:
      //! Get the parameters
      const ParamsDir& getParamsDir() const {return params;}

      //! Do the spin actions on the initial point source
      LatticePropagator deriv(const multi1d<LatticeColorMatrix>& u, const LatticePropagator& pt) const;

    private:
      ParamsDir  params;   /*!< source params */
    };


    //! Construct pion_1-(PionxD_T2) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  pion x D_T2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_4\gamma_5 D_i\f$  
     */
    class MesPionPionxDT2SrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionPionxDT2SrcConst(const ParamsDir& p) : params(p) {}

    protected:
      //! Get the parameters
      const ParamsDir& getParamsDir() const {return params;}

      //! Do the spin actions on the initial point source
      LatticePropagator deriv(const multi1d<LatticeColorMatrix>& u, const LatticePropagator& pt) const;

    private:
      ParamsDir  params;   /*!< source params */
    };


    //! Construct pion_1-(PionxB_T1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  pion x B_T1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_5 B_i\f$  
     */
    class MesPionPionxBT1SrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionPionxBT1SrcConst(const ParamsDir& p) : params(p) {}

    protected:
      //! Get the parameters
      const ParamsDir& getParamsDir() const {return params;}

      //! Do the spin actions on the initial point source
      LatticePropagator deriv(const multi1d<LatticeColorMatrix>& u, const LatticePropagator& pt) const;

    private:
      ParamsDir  params;   /*!< source params */
    };


    //! Construct pion_1-(RhoxB_T1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is rho x B_T1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \epsilon_{ijk}\gamma_j B_k\f$  
     */
    class MesPionRhoxBT1SrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionRhoxBT1SrcConst(const ParamsDir& p) : params(p) {}

    protected:
      //! Get the parameters
      const ParamsDir& getParamsDir() const {return params;}

      //! Do the spin actions on the initial point source
      LatticePropagator deriv(const multi1d<LatticeColorMatrix>& u, const LatticePropagator& pt) const;

    private:
      ParamsDir  params;   /*!< source params */
    };


    //! Construct pion_1-(RhoxB_T2) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  rho x B_T2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv s_{ijk}\gamma_j D_k\f$  
     */
    class MesPionRhoxBT2SrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionRhoxBT2SrcConst(const ParamsDir& p) : params(p) {}

    protected:
      //! Get the parameters
      const ParamsDir& getParamsDir() const {return params;}

      //! Do the spin actions on the initial point source
      LatticePropagator deriv(const multi1d<LatticeColorMatrix>& u, const LatticePropagator& pt) const;

    private:
      ParamsDir  params;   /*!< source params */
    };


    //! Construct pion_1-(A1xB_A1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  a1 x B_A1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_5 \gamma_i B_i\f$  
     */
    class MesPionA1xBA1SrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionA1xBA1SrcConst(const Params& p) : params(p) {}

    protected:
      //! Get the parameters
      const ParamsDir& getParamsDir() const {return params;}

      //! Do the spin actions on the initial point source
      LatticePropagator deriv(const multi1d<LatticeColorMatrix>& u, const LatticePropagator& pt) const;

    private:
      ParamsDir  params;   /*!< source params */
    };


    //! Construct pion_1-(RhoxB_T1) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  a11 x B_T1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_5 \epsilon_{ijk}\gamma_j B_k\f$  
     */
    class MesPionA1xBT1SrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionA1xBT1SrcConst(const ParamsDir& p) : params(p) {}

    protected:
      //! Get the parameters
      const ParamsDir& getParamsDir() const {return params;}

      //! Do the spin actions on the initial point source
      LatticePropagator deriv(const multi1d<LatticeColorMatrix>& u, const LatticePropagator& pt) const;

    private:
      ParamsDir  params;   /*!< source params */
    };


    //! Construct pion_1-(A1xB_T2) sequential source
    /*!
     * \ingroup hadron
     *
     * Operator is  a1 x B_T2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_5 s_{ijk}\gamma_j B_k\f$  
     */
    class MesPionA1xBT2SrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionA1xBT2SrcConst(const ParamsDir& p) : params(p) {}

    protected:
      //! Get the parameters
      const ParamsDir& getParamsDir() const {return params;}

      //! Do the spin actions on the initial point source
      LatticePropagator deriv(const multi1d<LatticeColorMatrix>& u, const LatticePropagator& pt) const;

    private:
      ParamsDir  params;   /*!< source params */
    };

  }  // end namespace


  //! Reader
  /*! @ingroup hadron */
  void read(XMLReader& xml, const string& path, DerivShellSourceConstEnv::Params& param);

  //! Writer
  /*! @ingroup hadron */
  void write(XMLWriter& xml, const string& path, const DerivShellSourceConstEnv::Params& param);


  //! Reader
  /*! @ingroup hadron */
  void read(XMLReader& xml, const string& path, DerivShellSourceConstEnv::ParamsDir& param);

  //! Writer
  /*! @ingroup hadron */
  void write(XMLWriter& xml, const string& path, const DerivShellSourceConstEnv::ParamsDir& param);


}  // end namespace Chroma

#endif
