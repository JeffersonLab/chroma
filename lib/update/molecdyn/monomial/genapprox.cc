// $Id: genapprox.cc,v 1.1 2005-02-03 03:16:41 edwards Exp $
/*! @file
 * @brief Wrapper for Remez code
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/remez.h"

namespace Chroma 
{ 
 
  //! Wrapper for Remez code
  /*! Compute partial fraction expansions for force, action and heat-bath */
  void generateApprox(RemezCoeff_t& fpfe, RemezCoeff_t& spfe, RemezCoeff_t& sipfe,
		      const Real& lower, const Real& upper,
		      int nth_root, int force_degree, int action_degree,
		      int digit_precision)
  {
    unsigned long prec = abs(digit_precision);
    unsigned long power_num = 1l;
    unsigned long power_den = abs(nth_root);

    if (nth_root > 0)
    {
      // Find approx to  x^(1/n)
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
      // Find approx to  x^(-1/n)
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


