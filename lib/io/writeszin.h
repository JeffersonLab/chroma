// $Id: writeszin.h,v 1.2 2003-10-08 04:37:50 edwards Exp $

/*! \file
 *  \brief Write a SZIN configuration written at configuration version 7.
 */

#ifndef __writeszin_h__
#define __writeszin_h__

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

//! Write a SZIN configuration file
/*!
 * \ingroup io
 *
 *   Gauge field layout is (fortran ordering)
 *     u(real/imag,color_row,color_col,site,cb,Nd)
 *         = u(2,Nc,Nc,VOL_CB,2,4)
 *
 *
 * \param header     structure holding config info ( Modify )
 * \param u          gauge configuration ( Read )
 * \param cfg_file   path ( Read )
 */    

void writeSzin(const SzinGauge_t& header, const multi1d<LatticeColorMatrix>& u, const string& cfg_file);

#endif
