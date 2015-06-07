// -*- C++ -*-
/*! \file
 *  \brief Hex smear a gauge field
 */

#ifndef __hex_smear_h__
#define __hex_smear_h__

namespace Chroma 
{ 
  //! Construct the "hex-smeared" links of Capitani et al. hep-lat/0607006
  /*!
   * \ingroup smear
   *
   * Arguments:
   *
   *  \param u		gauge field (Read)
   *  \param u_hex	"hex-smeared" gauge field (Write)
   *  \param nstep      number of hex smeared steps
   *
   */

  void Hex_Smear(const multi1d<LatticeColorMatrix>& u,
		 multi1d<LatticeColorMatrix>& u_hex, const int nstep) ;

}

#endif
