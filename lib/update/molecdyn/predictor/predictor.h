// -*- C++ -*-

/*! \file
 * \brief HMC linear system solver initial guess 
 *
 * Monomials for HMC
 */

/*! \defgroup predictor HMC linear system solver initial guess 
 * \ingroup molecdyn
 *
 * Monomials for HMC
 */

#ifndef __predictor_h__
#define __predictor_h__

#include "update/molecdyn/predictor/chrono_predictor.h"
#include "update/molecdyn/predictor/chrono_predictor_factory.h"
#include "update/molecdyn/predictor/predictor_aggregate.h"
#include "update/molecdyn/predictor/zero_guess_predictor.h"
#include "update/molecdyn/predictor/last_solution_predictor.h"

#include "update/molecdyn/predictor/circular_buffer.h"
#include "update/molecdyn/predictor/linear_extrap_predictor.h"

#include "update/molecdyn/predictor/lu_solve.h"
#include "update/molecdyn/predictor/mre_extrap_predictor.h"
#endif
