// -*- C++ -*-
// $Id: io.h,v 3.3 2007-02-04 22:08:39 edwards Exp $

#ifndef __io_h__
#define __io_h__

/*! \file
 * \brief IO routines
 *
 * Readers and writers of gauge fields and propagators
 */

/*! \defgroup io IO routines
 * \ingroup lib
 *
 * Readers and writers of gauge fields and propagators
 */

#include "enum_io/enum_io.h"

#include "xml_group_reader.h"
#include "cfgtype_io.h"

#include "readszin.h"
#include "qprop_io.h"
#include "szin_io.h"

#include "writeszin.h"

#include "gauge_io.h"
#include "milc_io.h"
#include "kyugauge_io.h"
#include "readmilc.h"
#include "writemilc.h"

#include "param_io.h"
#include "aniso_io.h"
#include "io_w.h"
#include "io_s.h"

#include "monomial_io.h"

#include "inline_io.h"

#endif
