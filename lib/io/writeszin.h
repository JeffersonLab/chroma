// $Id: writeszin.h,v 1.1 2003-09-25 22:20:28 edwards Exp $

#ifndef __writeszin_h__
#define __writeszin_h__

/*! \file
 *  \brief Write a SZIN configuration written at configuration version 7.
 */

//! Write a SZIN configuration file
/*!
 *   Gauge field layout is (fortran ordering)
 *     u(real/imag,color_row,color_col,site,cb,Nd)
 *         = u(2,Nc,Nc,VOL_CB,2,4)
 *
 *
 * \param xml        xml writer holding config info ( Read )
 * \param u          gauge configuration ( Read )
 * \param cfg_file   path ( Read )
 */    

void writeSzin(XMLBufferWriter& xml, const multi1d<LatticeColorMatrix>& u, const string& cfg_file);

#endif
