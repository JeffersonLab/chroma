// -*- C++ -*-
// $Id: deriv_sh_source_const_w.h,v 1.3 2006-02-20 20:42:24 edwards Exp $
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
    class MesPionxNablaT1SrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionxNablaT1SrcConst(const ParamsDir& p) : params(p) {}

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
    class MesA0xNablaT1SrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesA0xNablaT1SrcConst(const ParamsDir& p) : params(p) {}

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
    class MesA02xNablaT1SrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesA02xNablaT1SrcConst(const ParamsDir& p) : params(p) {}

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
    class MesRhoxNablaA1SrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesRhoxNablaA1SrcConst(const Params& p) : params(p) {}

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
    class MesRhoxNablaT1SrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesRhoxNablaT1SrcConst(const ParamsDir& p) : params(p) {}

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
    class MesRhoxNablaT2SrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesRhoxNablaT2SrcConst(const ParamsDir& p) : params(p) {}

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
    class MesA1xNablaA1SrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesA1xNablaA1SrcConst(const Params& p) : params(p) {}

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
    class MesA1xNablaT2SrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesA1xNablaT2SrcConst(const ParamsDir& p) : params(p) {}

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
    class MesA1xNablaESrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesA1xNablaESrcConst(const ParamsDir& p) : params(p) {}

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
    class MesB1xNablaT1SrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesB1xNablaT1SrcConst(const ParamsDir& p) : params(p) {}

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
    class MesA02xDT2SrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesA02xDT2SrcConst(const ParamsDir& p) : params(p) {}

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
    class MesA1xDA2SrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesA1xDA2SrcConst(const Params& p) : params(p) {}

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
    class MesA1xDESrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesA1xDESrcConst(const ParamsDir& p) : params(p) {}

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
    class MesA1xDT1SrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesA1xDT1SrcConst(const ParamsDir& p) : params(p) {}

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
    class MesA1xDT2SrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesA1xDT2SrcConst(const ParamsDir& p) : params(p) {}

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
    class MesB1xDA2SrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesB1xDA2SrcConst(const Params& p) : params(p) {}

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
    class MesB1xDESrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesB1xDESrcConst(const ParamsDir& p) : params(p) {}

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
    class MesB1xDT1SrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesB1xDT1SrcConst(const ParamsDir& p) : params(p) {}

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
    class MesB1xDT2SrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesB1xDT2SrcConst(const ParamsDir& p) : params(p) {}

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
    class MesRhoxDA2SrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesRhoxDA2SrcConst(const Params& p) : params(p) {}

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
    class MesRhoxDT1SrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesRhoxDT1SrcConst(const ParamsDir& p) : params(p) {}

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
    class MesRhoxDT2SrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesRhoxDT2SrcConst(const ParamsDir& p) : params(p) {}

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
    class MesPionxDT2SrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionxDT2SrcConst(const ParamsDir& p) : params(p) {}

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
    class MesPionxBT1SrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesPionxBT1SrcConst(const ParamsDir& p) : params(p) {}

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
    class MesRhoxBT1SrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesRhoxBT1SrcConst(const ParamsDir& p) : params(p) {}

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
    class MesRhoxBT2SrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesRhoxBT2SrcConst(const ParamsDir& p) : params(p) {}

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
    class MesA1xBA1SrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesA1xBA1SrcConst(const Params& p) : params(p) {}

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
    class MesA1xBT1SrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesA1xBT1SrcConst(const ParamsDir& p) : params(p) {}

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
    class MesA1xBT2SrcConst : public DerivSourceConst<LatticePropagator>
    {
    public:
      //! Full constructor
      MesA1xBT2SrcConst(const ParamsDir& p) : params(p) {}

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
