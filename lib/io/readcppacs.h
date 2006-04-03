// $Id: readcppacs.h,v 3.0 2006-04-03 04:58:55 edwards Exp $

/*! \file
 *  \brief Read a CPPACS gauge configuration 
 */

#ifndef __readcppacs_h__
#define __readcppacs_h__

#include "io/cppacs_io.h"

namespace Chroma {

//! Read a CPPACS gauge configuration 
/*!
 * \ingroup io
 *
 * \param xml        xml reader holding config info ( Modify )
 * \param u          gauge configuration ( Modify )
 * \param cfg_file   path ( Read )
 */    

void readCPPACS(XMLReader& xml, multi1d<LatticeColorMatrix>& u, const string& cfg_file);

//! Read a CPPACS gauge configuration written in the 1997 format
/*!
 * \ingroup io
 *
 * \param header     structure holding config info ( Modify )
 * \param u          gauge configuration ( Modify )
 * \param cfg_file   path ( Read )
 */    

void readCPPACS(CPPACSGauge_t& header, multi1d<LatticeColorMatrix>& u, const string& cfg_file);

}  // end namespace Chroma

#endif
