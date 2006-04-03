// $Id: writemilc.h,v 3.0 2006-04-03 04:58:56 edwards Exp $

/*! \file
 *  \brief  Write a MILC gauge configuration in the 1997 format
 */

#ifndef __writemilc_h__
#define __writemilc_h__

namespace Chroma {


//! Write a MILC gauge configuration in the 1997 format
/*!
 * \ingroup io
 *
 * \param xml        xml writer holding config info ( Read )
 * \param u          gauge configuration ( Read )
 * \param cfg_file   path ( Read )
 */    

void writeMILC(XMLBufferWriter& xml, const multi1d<LatticeColorMatrix>& u, 
	       const string& cfg_file);


//! Write a MILC gauge configuration in the 1997 format
/*!
 * \ingroup io
 *
 * \param header     structure holding config info ( Modify )
 * \param u          gauge configuration ( Read )
 * \param cfg_file   path ( Read )
 */    

void writeMILC(const MILCGauge_t& header, const multi1d<LatticeColorMatrix>& u, 
	       const string& cfg_file);

}  // end namespace Chroma

#endif
