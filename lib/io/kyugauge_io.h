// $Id: kyugauge_io.h,v 1.1 2004-04-15 03:23:49 edwards Exp $
/*! \file
 *  \brief Read/write a KYU gauge configuration
 */

#ifndef __kyugauge_io_h__
#define __kyugauge_io_h__

//! Read a Kentucky gauge configuration
/*!
 * \ingroup io
 *
 * \param u          gauge configuration ( Modify )
 * \param cfg_file   path ( Read )
 */    
void readKYU(multi1d<LatticeColorMatrix>& u, const string& cfg_file);

#endif
