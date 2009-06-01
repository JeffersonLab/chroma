// $Id: force_monitors.cc,v 3.4 2009-06-01 16:24:54 bjoo Exp $
/*! @file
 * @brief Helper function for calculating forces
 */

#include "update/molecdyn/monomial/force_monitors.h"
#include "util/gauge/taproj.h"

namespace Chroma 
{ 

 
  //! Writes a ForceCalc_t
  /*! @ingroup monomial */
  void write(XMLWriter& xml_out, const string& path, const ForceMonitors& param)
  {
    push(xml_out, path);
    
    write(xml_out, "F_sq_per_direction", param.F_sq_dir);
    write(xml_out, "F_avg_per_direction", param.F_avg_dir);
    write(xml_out, "F_max_per_direction", param.F_max_dir);

    write(xml_out, "F_sq", param.F_sq);
    write(xml_out, "F_avg", param.F_avg);
    write(xml_out, "F_max", param.F_max);

    pop(xml_out);
  }


  //! Helper function for calculating forces
  /*! @ingroup monomial */
  inline
  void forceMonitorCalc(const multi1d<LatticeColorMatrix>& F, ForceMonitors& forces)
  {
    START_CODE();

    //    ForceCalc_t forces;

    forces.F_sq  = zero;
    forces.F_avg = zero;
    forces.F_max = zero;

    forces.F_sq_dir.resize(Nd);
    forces.F_avg_dir.resize(Nd);
    forces.F_max_dir.resize(Nd);

    // Precompute these
    multi1d<LatticeReal> f2(F.size());
    multi1d<LatticeReal> f1(F.size());
    Real num_sites = QDP::Layout::vol();

    for(int mu=0; mu < F.size(); ++mu)
    {
      // This may be expensive
      LatticeColorMatrix F_tmp = F[mu];

      // Force will get taproj-ed at the end...
      // Make sure we look at the right hing
      taproj(F_tmp);

      // Jansens measure
      // -2 Tr F^2 = -2 ( -F^\dagger F ) = 2 F^\dagger F
      // the 2 counteracts taproj normalization maybe?
      f2[mu] = Real(-2)*real(trace(F_tmp*F_tmp));

      f1[mu] = sqrt(f2[mu]);
    }

    // Standard kind of sums
    for(int mu=0; mu < F.size(); ++mu)
    {
      // Get Square norms for direction mu - divide to get 'per site'
      forces.F_sq_dir[mu] = sum(f2[mu])/num_sites;

      // Get average norm for direction mu - divide to get 'per site'
      forces.F_avg_dir[mu] = sum(f1[mu])/num_sites;

      // Get max force for direction mu - this is already 'per site'
      forces.F_max_dir[mu] = globalMax(f1[mu]);
    
      //Sum up squares and averages
      forces.F_sq  += forces.F_sq_dir[mu];
      forces.F_avg += forces.F_avg_dir[mu];
    }

    // Divide F_sq and F_avg by the number of directions now (for average)
    Real fact_tmp = Real(1)/Real(Nd);
    forces.F_sq *= fact_tmp;
    forces.F_avg *= fact_tmp;

    // Find the maximum of the 4 directions.
    forces.F_max = forces.F_max_dir[0];
    for(int mu=1; mu < F.size(); ++mu)
    {
      if(  toBool( forces.F_max < forces.F_max_dir[mu] ) ) {
	forces.F_max = forces.F_max_dir[mu];
      }
    }

    END_CODE();

    return;
  }


  namespace ForceMonitorEnv { 
    static bool monitorForcesP = true;
  }

  void setForceMonitoring(bool monitorP) 
  {
    ForceMonitorEnv::monitorForcesP = monitorP;
  }


  void monitorForces(XMLWriter& xml_out, const string& path, const multi1d<LatticeColorMatrix>& F)
  {
    if( ForceMonitorEnv::monitorForcesP == true ) { 
      QDPIO::cout << "Monitoring force" << endl;
      ForceMonitors mon;
      forceMonitorCalc(F, mon);
      write(xml_out, path, mon);
    }
  }
}  //end namespace Chroma


