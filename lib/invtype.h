// -*- C++ -*-
// $Id: invtype.h,v 1.5 2004-05-25 21:47:39 bjoo Exp $

/*! @file
 * @brief Types of inverters
 */

#ifndef __invtype_h__
#define __invtype_h__

enum InvType {
  CG_INVERTER = 21, 
  MR_INVERTER = 22,
  BICG_INVERTER = 23,
  CR_INVERTER = 24,
  SUMR_INVERTER = 25,
  REL_CG_INVERTER = 26,
  REL_SUMR_INVERTER = 27,
  REL_GMRESR_SUMR_INVERTER = 28,
  REL_GMRESR_CG_INVERTER = 29,
  BICGSTAB_INVERTER = 30
};

#endif
