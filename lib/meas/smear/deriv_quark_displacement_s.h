// -*- C++ -*-
// $Id: deriv_quark_displacement_s.h,v 3.1 2006-11-17 02:55:11 edwards Exp $
/*! \file
 *  \brief Staggered Derivative displacements
 *
 *  There are a few variants of derivatives. One is
 *
 *  \f$\nabla_\mu f(x) = U_\mu(x)f(x+\mu) - U_{-\mu}(x)f(x-\mu)\f$
 */

#ifndef __deriv_quark_displacement_s_h__
#define __deriv_quark_displacement_s_h__

#include "meas/smear/quark_displacement.h"

namespace Chroma 
{

  //! Name and registration
  /*! @ingroup smear */
  namespace StaggeredDerivQuarkDisplacementEnv
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

  }  // end namespace


  //! Reader
  /*! @ingroup sources */
  void read(XMLReader& xml, const string& path, StaggeredDerivQuarkDisplacementEnv::Params& param);

  //! Writer
  /*! @ingroup sources */
  void write(XMLWriter& xml, const string& path, const StaggeredDerivQuarkDisplacementEnv::Params& param);


  //! Reader
  /*! @ingroup sources */
  void read(XMLReader& xml, const string& path, StaggeredDerivQuarkDisplacementEnv::ParamsDir& param);

  //! Writer
  /*! @ingroup sources */
  void write(XMLWriter& xml, const string& path, const StaggeredDerivQuarkDisplacementEnv::ParamsDir& param);


}  // end namespace Chroma

#endif
