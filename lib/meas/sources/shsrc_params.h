// -*- C++ -*-
// $Id: shsrc_params.h,v 1.1 2005-10-28 21:06:41 edwards Exp $
/*! \file
 *  \brief Shell source parameters
 *
 * This code is common to different fermions and propagators
 */

#ifndef __shsrc_params_h__
#define __shsrc_params_h__

#include "io/smearing_io.h"
#include "io/enum_io/enum_io.h"
#include "meas/sources/srcsnktype.h"
#include "meas/sources/wavetype.h"
#include "meas/smear/wvfkind.h"

namespace Chroma
{

  //! Point source parameters
  /*! @ingroup sources */
  struct ShellSourceConstParams
  {
    ShellSourceConstParams() {}
    ShellSourceConstParams(XMLReader& in, const std::string& path);
    
    WaveStateType    wave_state;    // S-wave or P-wave
    SmearingParam_t  sourceSmearParam;

    std::string      source_type;
    // wvf-function smearing type (Gaussian, Exponential, etc.)
    // smearing width
    // number of iteration for smearing
    int              j_decay;         // Decay direction
    int              direction;       // S-wave;   P-wave x(0) y(1) z(2)
    int              laplace_power;   // power=1 implies 1 laplacian
    int              disp_length;     // displacement length
    int              disp_dir;        // x(0), y(1), z(2)
    multi1d<int>     t_source;        // source location

    std::string      link_smearing;   /*!< link smearing xml */
  };


  //! Reader
  /*! @ingroup sources */
  void read(XMLReader& xml, const string& path, ShellSourceConstParams& param);

  //! Writer
  /*! @ingroup sources */
  void write(XMLWriter& xml, const string& path, const ShellSourceConstParams& param);

}  // end namespace Chroma


#endif
