// -*- C++ -*-
// $Id: smear.h,v 2.1 2005-11-07 06:40:55 edwards Exp $

/*! \file
 * \brief Smearing routines
 *
 * Central include file for all smearing routines
 */

/*! \defgroup smear Smearing routines
 * \ingroup meas
 *
 * Support for smearing of gauge and fermion fields.
 */

#ifndef __smear_h__
#define __smear_h__

#include "wvfkind.h"

#include "gaus_smear.h"
#include "laplacian.h"
#include "sink_smear2.h"
#include "stout_smear.h"
#include "hyp_smear.h"
#include "hyp_smear3d.h"
#include "ape_smear.h"
#include "displacement.h"

#include "quark_smearing.h"
#include "quark_source_sink.h"
#include "quark_smearing_factory.h"

#endif


