// -*- C++ -*-
// $Id: aniso_sym_gaugeact_params.h,v 3.2 2008-05-21 17:07:50 bjoo Exp $
/*! \file
 *  \brief Params for Anisotropic Symanzik Gauge Action
 *
 *  Tree-level LW with tapole improvement, missing 1x2 in time, also including
 *  2-plaq term. Taken from Morningstar-Peardon, hep-lat/9911003
 */

#ifndef __aniso_sym_gaugeact_params_h__
#define __aniso_sym_gaugeact_params_h__

#include "gaugeact.h"
#include "gaugebc.h"
#include "io/aniso_io.h"
namespace Chroma
{


  //! Parameter structure
  /*! @ingroup gaugeacts */
  struct AnisoSymGaugeActParams 
  {
    // Base Constructor
    AnisoSymGaugeActParams() {};
    
    // Read params from some root path
    AnisoSymGaugeActParams(XMLReader& xml_in, const std::string& path);

    Real beta;  //!< The beta  coupling
    Real u_s;   //!< Spatial Tadpole coupling
    Real u_t;   //!< Temporal Tadpole coupling
    bool use_subtraction; //!< Whether to use subtraction trick
    Real sub_zero; //!< Arbitrary constant (zero point energy)
    AnisoParam_t aniso;
  };
  
  /*! @ingroup gaugeacts */
  void read(XMLReader& xml, const string& path, AnisoSymGaugeActParams& param);
  
  /*! @ingroup gaugeacts */
  void write(XMLWriter& xml, const string& path, const AnisoSymGaugeActParams& param);

  //! Parameter structure
  /*! @ingroup gaugeacts */
  struct AnisoSymSpatialGaugeActParams 
  {
    // Base Constructor
    AnisoSymSpatialGaugeActParams() {};
    
    // Read params from some root path
    AnisoSymSpatialGaugeActParams(XMLReader& xml_in, const std::string& path);

    Real beta;  //!< The beta  coupling
    Real u_s;   //!< Spatial Tadpole coupling
    AnisoParam_t aniso; //!< The anisotropy parameters
    bool use_subtraction; //!< Whether to use subtraction trick
    Real sub_zero; //!< Arbitrary constant (zero point energy)
  };
  
  /*! @ingroup gaugeacts */
  void read(XMLReader& xml, const string& path, AnisoSymSpatialGaugeActParams& param);
  
  /*! @ingroup gaugeacts */
  void write(XMLWriter& xml, const string& path, const AnisoSymSpatialGaugeActParams& param);
  


};


#endif
