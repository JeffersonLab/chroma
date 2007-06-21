// -*- C++ -*-
// $Id: deriv_quark_displacement_w.h,v 3.3 2007-06-21 19:18:34 edwards Exp $
/*! \file
 *  \brief Derivative displacements
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

#ifndef __deriv_quark_displacement_w_h__
#define __deriv_quark_displacement_w_h__

#include "meas/smear/quark_displacement.h"

namespace Chroma 
{

  //! Name and registration
  /*! @ingroup smear */
  namespace DerivQuarkDisplacementEnv
  {
    bool registerAll();

  
    //! Params for derivative quark displacement
    /*! @ingroup smear */
    struct Params
    {
      Params();
      Params(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;

      std::string      displacement_type;    /*!< displacement type */

      int              deriv_length;         /*!< Displacement length in derivative */
    };

    //! Deriv meson source parameters
    /*! @ingroup sources */
    struct ParamsDir
    {
      ParamsDir();
      ParamsDir(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;

      std::string      displacement_type;    /*!< displacement type */

      int              deriv_dir;            /*!< Polarization direction */
      int              deriv_length;         /*!< Displacement length in derivative */
    };


    //! Construct (right Nabla) source
    /*!
     * \ingroup sources
     *
     * Operator is  Nabla
     * The sink interpolator structure is
     * \f$\Gamma_f \equiv \nabla_i\f$
     */
    template<typename T>
    class RightNablaDisplace : public QuarkDisplacement<T>
    {
    public:
      //! Full constructor
      RightNablaDisplace(const ParamsDir& p) : params(p) {}

      //! Displace the quark
      void operator()(T& quark, 
		      const multi1d<LatticeColorMatrix>& u,
		      enum PlusMinus isign) const;

    private:
      ParamsDir  params;   /*!< source params */
    };

    //! Construct (right D) source
    /*!
     * \ingroup sources
     *
     * Operator is  D
     * The sink interpolator structure is
     * \f$\Gamma_f \equiv D_i\f$
     */
    template<typename T>
    class RightDDisplace : public QuarkDisplacement<T>
    {
    public:
      //! Full constructor
      RightDDisplace(const ParamsDir& p) : params(p) {}

      //! Displace the quark
      void operator()(T& quark, 
		      const multi1d<LatticeColorMatrix>& u,
		      enum PlusMinus isign) const;

    private:
      ParamsDir  params;   /*!< source params */
    };

    //! Construct (right B) source
    /*!
     * \ingroup sources
     *
     * Operator is  B
     * The sink interpolator structure is
     * \f$\Gamma_f \equiv B_i\f$
     */
    template<typename T>
    class RightBDisplace : public QuarkDisplacement<T>
    {
    public:
      //! Full constructor
      RightBDisplace(const ParamsDir& p) : params(p) {}

      //! Displace the quark
      void operator()(T& quark, 
		      const multi1d<LatticeColorMatrix>& u,
		      enum PlusMinus isign) const;

    private:
      ParamsDir  params;   /*!< source params */
    };

    //! Construct (right E) source
    /*!
     * \ingroup sources
     *
     * Operator is  E
     * The sink interpolator structure is
     * \f$\Gamma_f \equiv E_i\f$
     */
    template<typename T>
    class RightEDisplace : public QuarkDisplacement<T>
    {
    public:
      //! Full constructor
      RightEDisplace(const ParamsDir& p) : params(p) {}

      //! Displace the quark
      void operator()(T& quark, 
		      const multi1d<LatticeColorMatrix>& u,
		      enum PlusMinus isign) const;

    private:
      ParamsDir  params;   /*!< source params */
    };


    //! Construct (right Laplacian) source
    /*!
     * \ingroup sources
     *
     * Operator is  Laplacian
     * The sink interpolator structure is
     * \f$\Gamma_f \equiv Laplacian\f$
     */
    template<typename T>
    class RightLapDisplace : public QuarkDisplacement<T>
    {
    public:
      //! Full constructor
      RightLapDisplace(const Params& p) : params(p) {}

      //! Displace the quark
      void operator()(T& quark, 
		      const multi1d<LatticeColorMatrix>& u,
		      enum PlusMinus isign) const;

    private:
      Params  params;   /*!< source params */
    };



    //! Construct (PionxNabla_T1) source
    /*!
     * \ingroup sources
     *
     * Operator is  Pion x Nabla_T1
     * The sink interpolator structure is
     * \f$\Gamma_f \equiv \gamma_5\nabla_i\f$
     */
    template<typename T>
    class MesPionxNablaT1Displace : public QuarkDisplacement<T>
    {
    public:
      //! Full constructor
      MesPionxNablaT1Displace(const ParamsDir& p) : params(p) {}

      //! Displace the quark
      void operator()(T& quark, 
		      const multi1d<LatticeColorMatrix>& u,
		      enum PlusMinus isign) const;

    private:
      ParamsDir  params;   /*!< source params */
    };


    //! Construct (A0xNabla_T1) source
    /*!
     * \ingroup sources
     *
     * Operator is  a0 x nabla_T1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \nabla_i\f$
     */
    template<typename T>
    class MesA0xNablaT1Displace : public QuarkDisplacement<T>
    {
    public:
      //! Full constructor
      MesA0xNablaT1Displace(const ParamsDir& p) : params(p) {}

      //! Displace the quark
      void operator()(T& quark, 
		      const multi1d<LatticeColorMatrix>& u,
		      enum PlusMinus isign) const;

    private:
      ParamsDir  params;   /*!< source params */
    };


    //! Construct (A0_2xNabla_T1) source
    /*!
     * \ingroup sources
     *
     * Operator is  a0_2 x nabla_T1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_4 \nabla_i\f$
     */
    template<typename T>
    class MesA02xNablaT1Displace : public QuarkDisplacement<T>
    {
    public:
      //! Full constructor
      MesA02xNablaT1Displace(const ParamsDir& p) : params(p) {}

      //! Displace the quark
      void operator()(T& quark, 
		      const multi1d<LatticeColorMatrix>& u,
		      enum PlusMinus isign) const;

    private:
      ParamsDir  params;   /*!< source params */
    };


    //! Construct (RhoxNabla_A1) source
    /*!
     * \ingroup sources
     *
     * Operator is  rho x nabla_A1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_i\nabla_i\f$
     */
    template<typename T>
    class MesRhoxNablaA1Displace : public QuarkDisplacement<T>
    {
    public:
      //! Full constructor
      MesRhoxNablaA1Displace(const Params& p) : params(p) {}

      //! Displace the quark
      void operator()(T& quark, 
		      const multi1d<LatticeColorMatrix>& u,
		      enum PlusMinus isign) const;

    private:
      Params  params;   /*!< source params */
    };


    //! Construct (RhoxNabla_T1) source
    /*!
     * \ingroup sources
     *
     * Operator is  rho x nabla_T1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \epsilon_{ijk}\gamma_j \nabla_k\f$  
     */
    template<typename T>
    class MesRhoxNablaT1Displace : public QuarkDisplacement<T>
    {
    public:
      //! Full constructor
      MesRhoxNablaT1Displace(const ParamsDir& p) : params(p) {}

      //! Displace the quark
      void operator()(T& quark, 
		      const multi1d<LatticeColorMatrix>& u,
		      enum PlusMinus isign) const;

    private:
      ParamsDir  params;   /*!< source params */
    };


    //! Construct (RhoxNabla_T2) source
    /*!
     * \ingroup sources
     *
     * Operator is  rho x nabla_T2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv s_{ijk}\gamma_j D_k\f$  
     */
    template<typename T>
    class MesRhoxNablaT2Displace : public QuarkDisplacement<T>
    {
    public:
      //! Full constructor
      MesRhoxNablaT2Displace(const ParamsDir& p) : params(p) {}

      //! Displace the quark
      void operator()(T& quark, 
		      const multi1d<LatticeColorMatrix>& u,
		      enum PlusMinus isign) const;

    private:
      ParamsDir  params;   /*!< source params */
    };


    //! Construct (A1xNabla_A1) source
    /*!
     * \ingroup sources
     *
     * Operator is  a1 x D_A1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_5\gamma_i \nabla_i\f$  
     */
    template<typename T>
    class MesA1xNablaA1Displace : public QuarkDisplacement<T>
    {
    public:
      //! Full constructor
      MesA1xNablaA1Displace(const Params& p) : params(p) {}

      //! Displace the quark
      void operator()(T& quark, 
		      const multi1d<LatticeColorMatrix>& u,
		      enum PlusMinus isign) const;

    private:
      Params  params;   /*!< source params */
    };


    //! Construct (A1xNabla_T2) source
    /*!
     * \ingroup sources
     *
     * Operator is  a1 x nabla_T2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_5 s_{ijk}\gamma_j \nabla_k\f$  
     */
    template<typename T>
    class MesA1xNablaT2Displace : public QuarkDisplacement<T>
    {
    public:
      //! Full constructor
      MesA1xNablaT2Displace(const ParamsDir& p) : params(p) {}

      //! Displace the quark
      void operator()(T& quark, 
		      const multi1d<LatticeColorMatrix>& u,
		      enum PlusMinus isign) const;

    private:
      ParamsDir  params;   /*!< source params */
    };


    //! Construct (A1xNabla_E) source
    /*!
     * \ingroup sources
     *
     * Operator is  a1 x nabla_E
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_5 S_{\alpha jk}\gamma_j \nabla_k\f$  
     */
    template<typename T>
    class MesA1xNablaEDisplace : public QuarkDisplacement<T>
    {
    public:
      //! Full constructor
      MesA1xNablaEDisplace(const ParamsDir& p) : params(p) {}

      //! Displace the quark
      void operator()(T& quark, 
		      const multi1d<LatticeColorMatrix>& u,
		      enum PlusMinus isign) const;

    private:
      ParamsDir  params;   /*!< source params */
    };


    //! Construct (B1xNabla_T1) source
    /*!
     * \ingroup sources
     *
     * Operator is  b1 x nabla_T1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_4\gamma_5\epsilon_{ijk}\gamma_j \nabla_k\f$  
     */
    template<typename T>
    class MesB1xNablaT1Displace : public QuarkDisplacement<T>
    {
    public:
      //! Full constructor
      MesB1xNablaT1Displace(const ParamsDir& p) : params(p) {}

      //! Displace the quark
      void operator()(T& quark, 
		      const multi1d<LatticeColorMatrix>& u,
		      enum PlusMinus isign) const;

    private:
      ParamsDir  params;   /*!< source params */
    };


    //! Construct (A0_2xD_T2) source
    /*!
     * \ingroup sources
     *
     * Operator is  a0_2 x D_T2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_4 D_i\f$  
     */
    template<typename T>
    class MesA02xDT2Displace : public QuarkDisplacement<T>
    {
    public:
      //! Full constructor
      MesA02xDT2Displace(const ParamsDir& p) : params(p) {}

      //! Displace the quark
      void operator()(T& quark, 
		      const multi1d<LatticeColorMatrix>& u,
		      enum PlusMinus isign) const;

    private:
      ParamsDir  params;   /*!< source params */
    };


    //! Construct (A1xD_A2) source
    /*!
     * \ingroup sources
     *
     * Operator is  a1 x D_A2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_5\gamma_i D_i\f$  
     */
    template<typename T>
    class MesA1xDA2Displace : public QuarkDisplacement<T>
    {
    public:
      //! Full constructor
      MesA1xDA2Displace(const Params& p) : params(p) {}

      //! Displace the quark
      void operator()(T& quark, 
		      const multi1d<LatticeColorMatrix>& u,
		      enum PlusMinus isign) const;

    private:
      Params  params;   /*!< source params */
    };


    //! Construct (A1xD_E) source
    /*!
     * \ingroup sources
     *
     * Operator is  a1 x D_E
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_5 S_{\alpha jk}\gamma_j D_k\f$  
     */
    template<typename T>
    class MesA1xDEDisplace : public QuarkDisplacement<T>
    {
    public:
      //! Full constructor
      MesA1xDEDisplace(const ParamsDir& p) : params(p) {}

      //! Displace the quark
      void operator()(T& quark, 
		      const multi1d<LatticeColorMatrix>& u,
		      enum PlusMinus isign) const;

    private:
      ParamsDir  params;   /*!< source params */
    };


    //! Construct (A1xD_T1) source
    /*!
     * \ingroup sources
     *
     * Operator is  a1 x D_T1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_5 s_{ijk}\gamma_j D_k\f$  
     */
    template<typename T>
    class MesA1xDT1Displace : public QuarkDisplacement<T>
    {
    public:
      //! Full constructor
      MesA1xDT1Displace(const ParamsDir& p) : params(p) {}

      //! Displace the quark
      void operator()(T& quark, 
		      const multi1d<LatticeColorMatrix>& u,
		      enum PlusMinus isign) const;

    private:
      ParamsDir  params;   /*!< source params */
    };


    //! Construct (A1xD_T2) source
    /*!
     * \ingroup sources
     *
     * Operator is  a1 x D_T2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_5\epsilon_{ijk}\gamma_j D_k\f$  
     */
    template<typename T>
    class MesA1xDT2Displace : public QuarkDisplacement<T>
    {
    public:
      //! Full constructor
      MesA1xDT2Displace(const ParamsDir& p) : params(p) {}

      //! Displace the quark
      void operator()(T& quark, 
		      const multi1d<LatticeColorMatrix>& u,
		      enum PlusMinus isign) const;

    private:
      ParamsDir  params;   /*!< source params */
    };


    //! Construct (B1xD_A2) source
    /*!
     * \ingroup sources
     *
     * Operator is  b1 x D_A2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_4\gamma_5 \gamma_i D_i\f$  
     */
    template<typename T>
    class MesB1xDA2Displace : public QuarkDisplacement<T>
    {
    public:
      //! Full constructor
      MesB1xDA2Displace(const Params& p) : params(p) {}

      //! Displace the quark
      void operator()(T& quark, 
		      const multi1d<LatticeColorMatrix>& u,
		      enum PlusMinus isign) const;

    private:
      Params  params;   /*!< source params */
    };

    //! Construct (B1xD_E) source
    /*!
     * \ingroup sources
     *
     * Operator is  b1 x D_E
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_4\gamma_5 S_{\alpha jk}\gamma_j D_k\f$  
     */
    template<typename T>
    class MesB1xDEDisplace : public QuarkDisplacement<T>
    {
    public:
      //! Full constructor
      MesB1xDEDisplace(const ParamsDir& p) : params(p) {}

      //! Displace the quark
      void operator()(T& quark, 
		      const multi1d<LatticeColorMatrix>& u,
		      enum PlusMinus isign) const;

    private:
      ParamsDir  params;   /*!< source params */
    };


    //! Construct (B1xD_T1) source
    /*!
     * \ingroup sources
     *
     * Operator is  b1 x D_T1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_4\gamma_5 s_{ijk}\gamma_j D_k\f$  
     */
    template<typename T>
    class MesB1xDT1Displace : public QuarkDisplacement<T>
    {
    public:
      //! Full constructor
      MesB1xDT1Displace(const ParamsDir& p) : params(p) {}

      //! Displace the quark
      void operator()(T& quark, 
		      const multi1d<LatticeColorMatrix>& u,
		      enum PlusMinus isign) const;

    private:
      ParamsDir  params;   /*!< source params */
    };


    //! Construct (B1xD_T2) source
    /*!
     * \ingroup sources
     *
     * Operator is  b1 x D_T2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_4\gamma_5 \epsilon_{ijk}\gamma_j D_k\f$  
     */
    template<typename T>
    class MesB1xDT2Displace : public QuarkDisplacement<T>
    {
    public:
      //! Full constructor
      MesB1xDT2Displace(const ParamsDir& p) : params(p) {}

      //! Displace the quark
      void operator()(T& quark, 
		      const multi1d<LatticeColorMatrix>& u,
		      enum PlusMinus isign) const;

    private:
      ParamsDir  params;   /*!< source params */
    };


    //! Construct (RhoxD_A2) source
    /*!
     * \ingroup sources
     *
     * Operator is  rho x D_A2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_i D_i\f$  
     */
    template<typename T>
    class MesRhoxDA2Displace : public QuarkDisplacement<T>
    {
    public:
      //! Full constructor
      MesRhoxDA2Displace(const Params& p) : params(p) {}

      //! Displace the quark
      void operator()(T& quark, 
		      const multi1d<LatticeColorMatrix>& u,
		      enum PlusMinus isign) const;

    private:
      Params  params;   /*!< source params */
    };


    //! Construct (RhoxD_T1) source
    /*!
     * \ingroup sources
     *
     * Operator is  rho x D_T1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv s_{ijk}\gamma_j D_k\f$  
     */
    template<typename T>
    class MesRhoxDT1Displace : public QuarkDisplacement<T>
    {
    public:
      //! Full constructor
      MesRhoxDT1Displace(const ParamsDir& p) : params(p) {}

      //! Displace the quark
      void operator()(T& quark, 
		      const multi1d<LatticeColorMatrix>& u,
		      enum PlusMinus isign) const;

    private:
      ParamsDir  params;   /*!< source params */
    };


    //! Construct (RhoxD_T2) source
    /*!
     * \ingroup sources
     *
     * Operator is  rho x D_T2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \epsilon_{ijk}\gamma_j D_k\f$  
     */
    template<typename T>
    class MesRhoxDT2Displace : public QuarkDisplacement<T>
    {
    public:
      //! Full constructor
      MesRhoxDT2Displace(const ParamsDir& p) : params(p) {}

      //! Displace the quark
      void operator()(T& quark, 
		      const multi1d<LatticeColorMatrix>& u,
		      enum PlusMinus isign) const;

    private:
      ParamsDir  params;   /*!< source params */
    };


    //! Construct (PionxD_T2) source
    /*!
     * \ingroup sources
     *
     * Operator is  pion x D_T2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_4\gamma_5 D_i\f$  
     */
    template<typename T>
    class MesPionxDT2Displace : public QuarkDisplacement<T>
    {
    public:
      //! Full constructor
      MesPionxDT2Displace(const ParamsDir& p) : params(p) {}

      //! Displace the quark
      void operator()(T& quark, 
		      const multi1d<LatticeColorMatrix>& u,
		      enum PlusMinus isign) const;

    private:
      ParamsDir  params;   /*!< source params */
    };


    //! Construct (PionxB_T1) source
    /*!
     * \ingroup sources
     *
     * Operator is  pion x B_T1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_5 B_i\f$  
     */
    template<typename T>
    class MesPionxBT1Displace : public QuarkDisplacement<T>
    {
    public:
      //! Full constructor
      MesPionxBT1Displace(const ParamsDir& p) : params(p) {}

      //! Displace the quark
      void operator()(T& quark, 
		      const multi1d<LatticeColorMatrix>& u,
		      enum PlusMinus isign) const;

    private:
      ParamsDir  params;   /*!< source params */
    };


    //! Construct (RhoxB_T1) source
    /*!
     * \ingroup sources
     *
     * Operator is rho x B_T1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \epsilon_{ijk}\gamma_j B_k\f$  
     */
    template<typename T>
    class MesRhoxBT1Displace : public QuarkDisplacement<T>
    {
    public:
      //! Full constructor
      MesRhoxBT1Displace(const ParamsDir& p) : params(p) {}

      //! Displace the quark
      void operator()(T& quark, 
		      const multi1d<LatticeColorMatrix>& u,
		      enum PlusMinus isign) const;

    private:
      ParamsDir  params;   /*!< source params */
    };


    //! Construct (RhoxB_T2) source
    /*!
     * \ingroup sources
     *
     * Operator is  rho x B_T2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv s_{ijk}\gamma_j D_k\f$  
     */
    template<typename T>
    class MesRhoxBT2Displace : public QuarkDisplacement<T>
    {
    public:
      //! Full constructor
      MesRhoxBT2Displace(const ParamsDir& p) : params(p) {}

      //! Displace the quark
      void operator()(T& quark, 
		      const multi1d<LatticeColorMatrix>& u,
		      enum PlusMinus isign) const;

    private:
      ParamsDir  params;   /*!< source params */
    };


    //! Construct (A1xB_A1) source
    /*!
     * \ingroup sources
     *
     * Operator is  a1 x B_A1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_5 \gamma_i B_i\f$  
     */
    template<typename T>
    class MesA1xBA1Displace : public QuarkDisplacement<T>
    {
    public:
      //! Full constructor
      MesA1xBA1Displace(const Params& p) : params(p) {}

      //! Displace the quark
      void operator()(T& quark, 
		      const multi1d<LatticeColorMatrix>& u,
		      enum PlusMinus isign) const;

    private:
      Params  params;   /*!< source params */
    };


    //! Construct (RhoxB_T1) source
    /*!
     * \ingroup sources
     *
     * Operator is  a11 x B_T1
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_5 \epsilon_{ijk}\gamma_j B_k\f$  
     */
    template<typename T>
    class MesA1xBT1Displace : public QuarkDisplacement<T>
    {
    public:
      //! Full constructor
      MesA1xBT1Displace(const ParamsDir& p) : params(p) {}

      //! Displace the quark
      void operator()(T& quark, 
		      const multi1d<LatticeColorMatrix>& u,
		      enum PlusMinus isign) const;

    private:
      ParamsDir  params;   /*!< source params */
    };


    //! Construct (A1xB_T2) source
    /*!
     * \ingroup sources
     *
     * Operator is  a1 x B_T2
     * The sink interpolator is   
     * \f$\Gamma_f \equiv \gamma_5 s_{ijk}\gamma_j B_k\f$  
     */
    template<typename T>
    class MesA1xBT2Displace : public QuarkDisplacement<T>
    {
    public:
      //! Full constructor
      MesA1xBT2Displace(const ParamsDir& p) : params(p) {}

      //! Displace the quark
      void operator()(T& quark, 
		      const multi1d<LatticeColorMatrix>& u,
		      enum PlusMinus isign) const;

    private:
      ParamsDir  params;   /*!< source params */
    };

  }  // end namespace


  //! Reader
  /*! @ingroup sources */
  void read(XMLReader& xml, const string& path, DerivQuarkDisplacementEnv::Params& param);

  //! Writer
  /*! @ingroup sources */
  void write(XMLWriter& xml, const string& path, const DerivQuarkDisplacementEnv::Params& param);


  //! Reader
  /*! @ingroup sources */
  void read(XMLReader& xml, const string& path, DerivQuarkDisplacementEnv::ParamsDir& param);

  //! Writer
  /*! @ingroup sources */
  void write(XMLWriter& xml, const string& path, const DerivQuarkDisplacementEnv::ParamsDir& param);


}  // end namespace Chroma

#endif
