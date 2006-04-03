// -*- C++ -*-
// $Id: integrator.h,v 3.0 2006-04-03 04:59:07 edwards Exp $

/*! \file
 * \brief HMD integrators
 *
 * Integrators
 */

/*! \defgroup integrator HMD integrators
 * \ingroup molecdyn
 *
 * Integrators
 */

#ifndef __integrator_h__
#define __integrator_h__

#include "update/molecdyn/integrator/abs_integrator.h"
#include "update/molecdyn/integrator/md_integrator_factory.h"
#include "update/molecdyn/integrator/lcm_integrator_leaps.h"
#include "update/molecdyn/integrator/lcm_pqp_leapfrog.h"
#include "update/molecdyn/integrator/lcm_sexton_weingarten.h"
#include "update/molecdyn/integrator/lcm_minimum_norm2_integrator.h"
#include "update/molecdyn/integrator/lcm_sw_min_mixed.h"

// Multi timescale integrators from Carsten
#include "update/molecdyn/integrator/lcm_minimum_norm2_integrator_mts.h"
#include "update/molecdyn/integrator/lcm_minimum_norm2_qpq_integrator_mts.h"
#include "update/molecdyn/integrator/lcm_pqp_leapfrog_mts.h"

#include "update/molecdyn/integrator/integrator_aggregate.h"
#endif
