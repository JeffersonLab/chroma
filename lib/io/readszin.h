// $Id: readszin.h,v 1.6 2003-08-27 22:08:41 edwards Exp $

#ifndef __readszin_h__
#define __readszin_h__

/*! \file
 *  \brief Read in a configuration written by SZIN up to configuration version 7.
 */

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

#endif
