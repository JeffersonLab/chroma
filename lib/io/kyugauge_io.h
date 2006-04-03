// $Id: kyugauge_io.h,v 3.0 2006-04-03 04:58:55 edwards Exp $
/*! \file
 *  \brief Read/write a KYU gauge configuration
 */

#ifndef __kyugauge_io_h__
#define __kyugauge_io_h__

namespace Chroma {

//! Read a Kentucky gauge configuration
/*!
 * \ingroup io
 *
 * \param u          gauge configuration ( Modify )
 * \param cfg_file   path ( Read )
 */    
void readKYU(multi1d<LatticeColorMatrix>& u, const string& cfg_file);


//! Write a Kentucky gauge configuration
/*!
 * \ingroup io
 *
 * \param u          gauge configuration ( Read )
 * \param cfg_file   path ( Read )
 */    

void writeKYU(const multi1d<LatticeColorMatrix>& u, const string& cfg_file);

}  // end namespace Chroma

#endif
