// -*- C++ -*-
/*! \file
 * \brief Hadron 2pt correlators
 */

#ifndef __key_hadron_2pt_corr_h__
#define __key_hadron_2pt_corr_h__

#include "chromabase.h"

namespace Chroma
{

  //----------------------------------------------------------------------------
  //! Key for Hadron 2pt corr
  struct KeyHadron2PtCorr_t
  {
    int          num_vecs;    /*!< Number of vectors used in this corr */

    std::string  src_name;    /*!< Some string label for the operator */
    std::string  src_smear;   /*!< Some string label for the smearing of this operator */
    multi1d<int> src_lorentz; /*!< Source Lorentz indices */
    int          src_spin;    /*!< Source Dirac spin indices */

    std::string  snk_name;    /*!< Some string label for the operator */
    std::string  snk_smear;   /*!< Some string label for the smearing of this operator */
    multi1d<int> snk_lorentz; /*!< Sink Lorentz indices */
    int          snk_spin;    /*!< Sink Dirac spin indices */

    multi1d<int> mom;         /*!< D-1 momentum of the sink operator */
    std::string  mass;        /*!< Some string label for the mass(es) in the corr */
    std::string  ensemble;    /*!< Label for the ensemble */
  };

  //----------------------------------------------------------------------------
  //! Used for error output
//  std::ostream& operator<<(std::ostream& os, const KeyHadron2PtCorr_t& d);
	
  //----------------------------------------------------------------------------
  //! KeyHadron2PtCorr reader
  void read(XMLReader& xml, const std::string& path, KeyHadron2PtCorr_t& param);

  //! KeyHadron2PtCorr writer
  void write(XMLWriter& xml, const std::string& path, const KeyHadron2PtCorr_t& param);

  //----------------------------------------------------------------------------
  //! KeyHadron2PtCorr reader
  void read(BinaryReader& bin, KeyHadron2PtCorr_t& param);

  //! KeyHadron2PtCorr writer
  void write(BinaryWriter& bin, const KeyHadron2PtCorr_t& param);

} // namespace ColorVec

#endif
