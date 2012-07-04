// -*- C++ -*-
/*! \file
 * \brief Hadron 3pt correlators
 */

#ifndef __key_hadron_3pt_corr_h__
#define __key_hadron_3pt_corr_h__

#include "chromabase.h"

namespace Chroma
{

  //----------------------------------------------------------------------------
  //! Hold momenta
  /*! \ingroup ferm */
  struct PiPf
  {
    multi1d<int>  p_i;      /*!< Source momentum */
    multi1d<int>  p_f;      /*!< Sink momentum */
  };

  //! Reader
  /*! \ingroup ferm */
  void read(XMLReader& xml, const std::string& path, PiPf& val);

  //! Writer
  /*! \ingroup ferm */
  void write(XMLWriter& xml, const std::string& path, const PiPf& val);


  //----------------------------------------------------------------------------
  //! Key for Hadron 3pt corr
  /*! \ingroup ferm */
  struct KeyHadron3PtCorr_t
  {
    int          num_vecs;    /*!< Number of vectors used in this corr */

    PiPf         pi_pf;       /*!< Source and sink momenta */
    int          gamma;       /*!< Gamma matrix index (0 .. 15). In DP basis */
    multi1d<int> links;       /*!< Gauge link insertions */

    int          dt;          /*!< Source-sink separation */
    int          quark;       /*!< Some number indicating which quark line */

    std::string  src_name;    /*!< Some string label for the operator */
    std::string  src_smear;   /*!< Some string label for the smearing of this operator */
    multi1d<int> src_lorentz; /*!< Source Lorentz indices */
    int          src_spin;    /*!< Source Dirac spin indices */

    std::string  snk_name;    /*!< Some string label for the operator */
    std::string  snk_smear;   /*!< Some string label for the smearing of this operator */
    multi1d<int> snk_lorentz; /*!< Sink Lorentz indices */
    int          snk_spin;    /*!< Sink Dirac spin indices */

    std::string  mass;        /*!< Some string label for the mass(es) in the corr */
    std::string  ensemble;    /*!< Label for the ensemble */
  };

  //----------------------------------------------------------------------------
  //! Used for error output
//  std::ostream& operator<<(std::ostream& os, const KeyHadron3PtCorr_t& d);
	
  //----------------------------------------------------------------------------
  //! KeyHadron3PtCorr reader
  /*! \ingroup ferm */
  void read(XMLReader& xml, const std::string& path, KeyHadron3PtCorr_t& param);

  //! KeyHadron3PtCorr writer
  /*! \ingroup ferm */
  void write(XMLWriter& xml, const std::string& path, const KeyHadron3PtCorr_t& param);

  //----------------------------------------------------------------------------
  //! KeyHadron3PtCorr reader
  /*! \ingroup ferm */
  void read(BinaryReader& bin, KeyHadron3PtCorr_t& param);

  //! KeyHadron3PtCorr writer
  /*! \ingroup ferm */
  void write(BinaryWriter& bin, const KeyHadron3PtCorr_t& param);

} // namespace ColorVec

#endif
