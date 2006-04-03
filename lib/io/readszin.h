// $Id: readszin.h,v 3.0 2006-04-03 04:58:55 edwards Exp $

/*! \file
 *  \brief Read in a configuration written by SZIN up to configuration version 7.
 */

#ifndef __readszin_h__
#define __readszin_h__

#include "io/szin_io.h"

namespace Chroma {

//! Read a SZIN configuration file
/*!
 *   Gauge field layout is (fortran ordering)
 *     u(real/imag,color_row,color_col,site,cb,Nd)
 *         = u(2,Nc,Nc,VOL_CB,2,4)
 *
 *
 * \param xml        xml reader holding config info ( Modify )
 * \param u          gauge configuration ( Modify )
 * \param cfg_file   path ( Read )
 */    

void readSzin(XMLReader& xml, multi1d<LatticeColorMatrix>& u, const string& cfg_file);

//! Read a SZIN configuration file
/*!
 * \ingroup io
 *
 *   Gauge field layout is (fortran ordering)
 *     u(real/imag,color_row,color_col,site,cb,Nd)
 *         = u(2,Nc,Nc,VOL_CB,2,4)
 *
 *
 * \param header     structure holding config info ( Modify )
 * \param u          gauge configuration ( Modify )
 * \param cfg_file   path ( Read )
 */    

void readSzin(SzinGauge_t& header, multi1d<LatticeColorMatrix>& u, const string& cfg_file);

}  // end namespace Chroma

#endif
