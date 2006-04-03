// -*- C++ -*-
// $Id: fermbcs.h,v 3.0 2006-04-03 04:58:48 edwards Exp $

/*! \file
 * \brief Fermion boundary conditions
 *
 * Support for fermion boundary conditions
 */

/*! \defgroup fermbcs Fermion boundary conditions
 * \ingroup fermact
 *
 * Support for fermion boundary conditions
 */

#ifndef __fermbcs_h__
#define __fermbcs_h__

#include "periodic_fermbc.h"
#include "simple_fermbc.h"

#include "fermbcs_reader_w.h"
#include "fermbc_factory_w.h"
#include "simple_fermbc_w.h"
#include "periodic_fermbc_w.h"
#include "twisted_fermbc_w.h"
#include "schroedinger_fermbc_w.h"
#include "schr_sf_fermbc_w.h"
#include "schr_triv_fermbc_w.h"
#include "schr_nonpert_fermbc_w.h"
#include "schr_coupling_fermbc_w.h"
#include "schr_chromomag_fermbc_w.h"
#include "schr_dirich_fermbc_w.h"

#include "fermbcs_reader_s.h"
#include "fermbc_factory_s.h"
#include "simple_fermbc_s.h"

#endif


