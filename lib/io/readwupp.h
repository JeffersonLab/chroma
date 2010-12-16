// $Id: readszin.h,v 3.0 2006/04/03 04:58:55 edwards Exp $
/*! \file
 *  \brief Read in a configuration written by Wupp up to configuration version 7.
 */

#ifndef __readwupp_h__
#define __readwupp_h__

namespace Chroma {

//! Read an expanded BMW configuration file
/*!
 *
 *
 * \param xml        xml reader holding config info ( Modify )
 * \param u          gauge configuration ( Modify )
 * \param cfg_file   path ( Read )
 */    

void readWupp(XMLReader& xml, multi1d<LatticeColorMatrix>& u, const string& cfg_file);

//! Read a WUPP configuration file
/*!
 * \ingroup io
 *   !!!!!!!!INCORRECT, CHANGE!!!!!!!!!
 *   Gauge field layout is (fortran ordering)
 *     u(real/imag,color_row,color_col,site,cb,Nd)
 *         = u(2,Nc,Nc,VOL_CB,2,4)
 *
 *
 * \param header     structure holding config info ( Modify )
 * \param u          gauge configuration ( Modify )
 * \param cfg_file   path ( Read )
 */    


}  // end namespace Chroma

#endif
