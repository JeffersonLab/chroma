// -*- C++ -*-
// $Id: fermacts_w.h,v 1.16 2004-11-02 10:33:49 bjoo Exp $

/*! \file
 * \brief Fermion actions
 *
 * Various fermion actions
 */

/*! \defgroup fermact Fermion actions
 * \ingroup actions
 *
 * Various fermion actions
 */

#ifndef __fermacts_w_h__
#define __fermacts_w_h__

#include "fermfactory_w.h"

#include "unprec_wilson_fermact_w.h"
#include "prec_wilson_fermact_w.h"

#include "unprec_parwilson_fermact_w.h"
#include "prec_parwilson_fermact_w.h"

#include "unprec_dwf_fermact_base_array_w.h"
#include "unprec_dwf_fermact_array_w.h"
#include "unprec_nef_fermact_array_w.h"
#include "unprec_zolo_nef_fermact_array_w.h"
#include "prec_zolo_nef_fermact_array_w.h"
#include "unprec_ovdwf_fermact_array_w.h"
#include "unprec_ovext_fermact_array_w.h"

#include "prec_dwf_fermact_base_array_w.h"
#include "prec_dwf_fermact_array_w.h"
#include "prec_nef_fermact_array_w.h"
#include "prec_ovdwf_fermact_array_w.h"

#if defined(BUILD_SSE_DWF_CG)
#include "prec_dwf_fermact_array_sse_w.h"
#endif

#include "overlap_fermact_base_w.h"
#include "overlap_state.h"
#include "ovlap_partfrac4d_fermact_w.h"
#include "unprec_ovlap_contfrac5d_fermact_array_w.h"
#include "prec_ovlap_contfrac5d_fermact_array_w.h"

#include "unprec_dwftransf_fermact_w.h"


#endif


