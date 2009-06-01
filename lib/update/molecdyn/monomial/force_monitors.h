// -*- C++ -*-
/*! @file
 * @brief Helper function for calculating forces
 */

#ifndef __force_monitors_h__
#define __force_monitors_h__

#include "chromabase.h"

namespace Chroma
{
  //! Diagnostics about the forces, per direction and total
  /*! @ingroup monomial */
  struct ForceMonitors
  {
    Real   F_sq;     /*!< sum norm2(F) */
    Real   F_avg;    /*!< sum sqrt(norm2(F)) */
    Real   F_max;    /*!< max(localNorm2(F)) */
    multi1d<Real> F_sq_dir;
    multi1d<Real> F_avg_dir;
    multi1d<Real> F_max_dir;
  };


  //! Writes a ForceCalc_t
  /*! @ingroup monomial */
  void write(XMLWriter& xml_out, const string& path, const ForceMonitors& param);

  //! Helper function for monitoring forces
  /*! @ingroup monomial */
  void forceMonitorCalc(const multi1d<LatticeColorMatrix>& F, ForceMonitors& forces);

  //! Calculate and write out forces
  /*! @ingroup monomial */
  void monitorForces(XMLWriter& xml_out, const string& path, const multi1d<LatticeColorMatrix>& F);

  void setForceMonitoring(bool monitorP);


}
#endif
