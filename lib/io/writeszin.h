// $Id: writeszin.h,v 3.0 2006-04-03 04:58:56 edwards Exp $

/*! \file
 *  \brief Write a SZIN configuration written at configuration version 7.
 */

#ifndef __writeszin_h__
#define __writeszin_h__

#include "io/szin_io.h"

namespace Chroma 
{

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



  //! Write a truncated SZIN configuration file
  /*!
   * \ingroup io
   *
   * \param header     structure holding config info ( Modify )
   * \param u          gauge configuration ( Read )
   * \param j_decay    direction which will be truncated ( Read )
   * \param t_start    starting slice in j_decay direction ( Read )
   * \param t_end      ending slice in j_decay direction ( Read )
   * \param cfg_file   path ( Read )
   */    

  void writeSzinTrunc(const SzinGauge_t& header, const multi1d<LatticeColorMatrix>& u, 
		      int j_decay, int t_start, int t_end, 
		      const string& cfg_file);


  //! Write a replicated (in time direction) SZIN configuration file
  /*!
   * \ingroup io
   *
   * \param header     structure holding config info ( Modify )
   * \param u          gauge configuration ( Read )
   * \param j_decay    direction for replication ( Read )
   * \param n_replica  number of replicas in j_decay direction ( Read )
   * \param cfg_file   path ( Read )
   */    

  void writeSzinReplica(SzinGauge_t& header, const multi1d<LatticeColorMatrix>& u, 
			int j_decay, int n_replica, 
			const string& cfg_file);

}  // end namespace Chroma

#endif

