// -*- C++ -*-
// $Id: invtype.h,v 1.2 2004-05-12 15:45:10 bjoo Exp $

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
  SUMR_INVERTER = 25
};

#endif
