// $Id: readmilc.h,v 1.2 2003-10-10 03:46:46 edwards Exp $

/*! \file
 *  \brief Read a MILC gauge configuration written in the 1997 format
 */

#ifndef __readmilc_h__
#define __readmilc_h__

#include "io/milc_io.h"

//! Read a MILC gauge configuration written in the 1997 format
/*!
 * \ingroup io
 *
 * \param xml        xml reader holding config info ( Modify )
 * \param u          gauge configuration ( Modify )
 * \param cfg_file   path ( Read )
 */    

void readMILC(XMLReader& xml, multi1d<LatticeColorMatrix>& u, const string& cfg_file);

//! Read a MILC gauge configuration written in the 1997 format
/*!
 * \ingroup io
 *
 * \param header     structure holding config info ( Modify )
 * \param u          gauge configuration ( Modify )
 * \param cfg_file   path ( Read )
 */    

void readMILC(MILCGauge_t& header, multi1d<LatticeColorMatrix>& u, const string& cfg_file);

#endif
