/*
 * syssolver_mdagm_mrhs.h
 *
 *  Created on: Mar 11, 2019
 *      Author: bjoo
 */

#ifndef LIB_ACTIONS_FERM_INVERT_SYSSOLVER_MDAGM_MRHS_H_
#define LIB_ACTIONS_FERM_INVERT_SYSSOLVER_MDAGM_MRHS_H_


#include "linearop.h"
#include "handle.h"
#include "syssolver.h"

namespace Chroma
{
  //! SystemSolver disambiguator
  /*! This struct is solely to disambiguate the type of SystemSolvers */

  /* NB: Previously these were declared as 'virtual public SystemSolver<T>'
     BUT That seemed to break XLC in a Bad Way */
  template<typename T>
  class MdagMMRHSSystemSolver : public MultiRHSSystemSolver<T>
  {
  };

}



#endif /* LIB_ACTIONS_FERM_INVERT_SYSSOLVER_MDAGM_MRHS_H_ */
