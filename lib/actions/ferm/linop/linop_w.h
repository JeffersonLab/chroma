// -*- C++ -*-
// $Id: linop_w.h,v 3.4 2008-10-08 19:40:17 bjoo Exp $

/*! \file
 * \brief Linear operators
 *
 * Various fermion linear operators
 */

#ifndef __linop_w_h__
#define __linop_w_h__

// the file dslash_w.h switches in the right dslash op for
// WilsonDslash
#include "dslash_w.h"
#include "dslash_array_w.h"

#include "unprec_wilson_linop_w.h"
#include "unprec_parwilson_linop_w.h"
#include "unprec_dwf4d_linop_w.h"
#include "unprec_pdwf4d_linop_w.h"
#include "unprec_ppdwf4d_linop_w.h"
#include "unprec_clover_linop_w.h"

#include "eoprec_wilson_linop_w.h"
#include "eoprec_clover_linop_w.h"
#include "eoprec_parwilson_linop_w.h"

#include "unprec_s_cprec_t_wilson_linop_w.h"

#include "iluprec_s_cprec_t_wilson_linop_w.h"
#include "iluprec_s_cprec_t_clover_linop_w.h"

#include "ilu2prec_s_cprec_t_wilson_linop_w.h"
#include "ilu2prec_s_cprec_t_clover_linop_w.h"

#include "unprec_dwflike_linop_base_array_w.h"
#include "unprec_dwf_linop_array_w.h"
#include "unprec_nef_linop_array_w.h"
#include "unprec_ovdwf_linop_array_w.h"
#include "unprec_ovext_linop_array_w.h"
#include "eoprec_ovext_linop_array_w.h"

#include "eoprec_dwflike_linop_base_array_w.h"
#include "eoprec_dwf_linop_array_w.h"
#include "eoprec_nef_linop_array_w.h"
#include "eoprec_nef_general_linop_array_w.h"
#include "eoprec_ovdwf_linop_array_w.h"

#include "dwffld_w.h"

#include "lDeltaLs_w.h"

#include "lovddag_w.h"
#include "lovddag_double_pass_w.h"
#include "lovlapms_w.h"
#include "lovlap_double_pass_w.h"
#include "lg5eps_w.h"

#include "unprec_ovlap_contfrac5d_linop_array_w.h"
#include "eoprec_ovlap_contfrac5d_linop_base_array_w.h"
#include "unprec_ovlap_contfrac5d_nonhermop_array_w.h"
#include "unprec_ht_contfrac5d_linop_array_w.h"
#include "eoprec_ht_contfrac5d_linop_array_w.h"

#include "ovlap_contfrac5d_w.h"

#include "lgherm_w.h"
#include "lg5Rherm_w.h"

#include "unprec_dwftransf_linop_w.h"

#endif


