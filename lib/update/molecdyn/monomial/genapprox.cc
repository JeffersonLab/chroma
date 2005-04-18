// $Id: genapprox.cc,v 1.3 2005-04-18 16:23:23 edwards Exp $
/*! @file
 * @brief Wrapper for Remez code
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/remez.h"
#include "update/molecdyn/monomial/genapprox.h"

namespace Chroma 
{ 
 
  //! Wrapper for Remez code
  /*! Compute partial fraction expansions for force, action and heat-bath */
  void generateApprox(RemezCoeff_t& fpfe, RemezCoeff_t& spfe, RemezCoeff_t& sipfe,
		      const Real& lower, const Real& upper,
		      int num_frac, int den_frac, int force_degree, int action_degree,
		      int digit_precision)
  {
    unsigned long prec = abs(digit_precision);
    unsigned long power_num = abs(num_frac);
    unsigned long power_den = abs(den_frac);

    if (den_frac <=0)
      QDP_error_exit("%s: invalid params", __func__);

    if (num_frac > 0)
    {
      // Find approx to  x^(num_frac/den_frac)
      // Force PFE
      QDPIO::cout << "Compute partial fraction expansion for force calc." << endl;
      Remez  remez(lower, upper, prec);
      remez.generateApprox(force_degree, power_num, power_den);
      fpfe = remez.getPFE();

      // Action and heat-bath PFE 
      QDPIO::cout << "Compute partial fraction expansion for action and heatbath" << endl;
      remez.generateApprox(action_degree, power_num, 2*power_den);
      spfe = remez.getPFE();
      sipfe = remez.getIPFE();
    }
    else
    {
      // Find approx to  x^(-num_frac/den_frac)
      // Force PFE
      QDPIO::cout << "Compute partial fraction expansion for force calc." << endl;
      Remez  remez(lower, upper, prec);
      remez.generateApprox(force_degree, power_num, power_den);
      fpfe = remez.getIPFE();

      // Action and heat-bath PFE 
      QDPIO::cout << "Compute partial fraction expansion for action and heatbath" << endl;
      remez.generateApprox(action_degree, power_num, 2*power_den);
      spfe = remez.getIPFE();
      sipfe = remez.getPFE();
    }

  }

};  //end namespace Chroma


