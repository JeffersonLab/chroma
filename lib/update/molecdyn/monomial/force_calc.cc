// $Id: force_calc.cc,v 1.1 2007-03-15 13:35:48 edwards Exp $
/*! @file
 * @brief Helper function for calculating forces
 */

#include "update/molecdyn/monomial/force_calc.h"

namespace Chroma 
{ 
 
  //! Writes a ForceCalc_t
  /*! @ingroup molecdyn */
  void write(XMLWriter& xml_out, const string& path, const ForceCalc_t& param)
  {
    push(xml_out, path);

    write(xml_out, "F_sq", param.F_sq);
    write(xml_out, "F_avg", param.F_avg);
    write(xml_out, "F_max", param.F_max);

    pop(xml_out);
  }


  //! Helper function for calculating forces
  /*! @ingroup molecdyn */
  ForceCalc_t forceCalc(const multi1d<LatticeColorMatrix>& F)
  {
    START_CODE();

    ForceCalc_t forces;

    forces.F_sq  = zero;
    forces.F_avg = zero;
    forces.F_max = zero;

    // Precompute these
    multi1d<LatticeReal> f2(F.size());
    multi1d<LatticeReal> f1(F.size());
    for(int mu=0; mu < F.size(); ++mu)
    {
      f2[mu] = localNorm2(F[mu]);
      f1[mu] = sqrt(f2[mu]);
    }

    // Standard kind of sums
    for(int mu=0; mu < F.size(); ++mu)
    {
      forces.F_sq  += sum(f2[mu]);
      forces.F_avg += sum(f1[mu]);
    }

    // Unroll for global max
    forces.F_max = globalMax(f1[0]);
    for(int mu=1; mu < F.size(); ++mu)
    {
      Real fmax = globalMax(f1[mu]);

      if (toBool(forces.F_max < fmax))
	forces.F_max = fmax;
    }

    END_CODE();

    return forces;
  }

}  //end namespace Chroma


