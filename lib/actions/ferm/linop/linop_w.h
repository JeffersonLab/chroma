// -*- C++ -*-
// $Id: linop_w.h,v 1.16 2004-12-10 02:54:24 edwards Exp $

/*! \file
 * \brief Linear operators
 *
 * Various fermion linear operators
 */

/*! \defgroup linop Linear operators
 * \ingroup actions
 *
 * Various fermion linear operators
 */

#ifndef __linop_w_h__
#define __linop_w_h__

// the file dslash_w.h switches in the right dslash op for
// WilsonDslash
#include "dslash_w.h"

#include "unprec_wilson_linop_w.h"
#include "prec_wilson_linop_w.h"
#include "overlapbu_linop_w.h"

#include "unprec_parwilson_linop_w.h"
#include "prec_parwilson_linop_w.h"

#include "unprec_dwf4d_linop_w.h"
#include "unprec_pdwf4d_linop_w.h"
#include "unprec_ppdwf4d_linop_w.h"

#include "unprec_dwf_linop_base_array_w.h"
#include "unprec_dwf_linop_array_w.h"
#include "unprec_nef_linop_array_w.h"
#include "unprec_ovdwf_linop_array_w.h"
#include "unprec_ovext_linop_array_w.h"

#include "prec_dwf_linop_base_array_w.h"
#include "prec_dwf_linop_array_w.h"
#include "prec_nef_linop_array_w.h"
#include "prec_nef_general_linop_array_w.h"
#include "prec_ovdwf_linop_array_w.h"

#include "dwffld_w.h"

#include "lDeltaLs_w.h"

#include "lovddag_w.h"
#include "lovddag_double_pass_w.h"
#include "lovlapms_w.h"
#include "lovlap_double_pass_w.h"
#include "lg5eps_w.h"

#include "unprec_ovlap_contfrac5d_linop_array_w.h"
#include "unprec_ovlap_contfrac5d_nonhermop_array_w.h"
#include "prec_ovlap_contfrac5d_linop_array_w.h"
#include "lgherm_w.h"

#include "unprec_dwftransf_linop_w.h"


#endif


