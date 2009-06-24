// -*- C++ -*-
// $Id: sources.h,v 3.4 2009-06-24 19:59:33 jbulava Exp $

/*! \file
 * \brief Quark sources
 *
 * Central include file for all quark source construction
 */

/*! \defgroup sources Quark sources 
 * \ingroup meas
 *
 *  This include file supports various types of quark sources
 */

#ifndef __sources_h__
#define __sources_h__

#include "srcfil.h"
#include "z2_src.h"
#include "srcfil.h"
#include "zN_src.h"

#include "source_construction.h"
#include "source_const_factory.h"
#include "source_smearing_factory.h"

#include "source_const_aggregate.h"
#include "source_smearing_aggregate.h"

#include "pt_source_const.h"
#include "sh_source_const.h"
#include "rndz2wall_source_const.h"
#include "dilutezN_source_const.h"
#include "dilute_zN_eigvec_source_const.h"
#include "diluteGrid_source_const.h"

#include "pt_source_smearing.h"
#include "sh_source_smearing.h"

#include "sources_s.h"
#include "sources_w.h"

#endif
