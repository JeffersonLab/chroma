// -*- C++ -*-
/*! \file
 * \brief Sparse matrix representation of spin matrices
 */

#ifndef __spin_rep_h__
#define __spin_rep_h__

#include "chromabase.h"

namespace Chroma
{
  //----------------------------------------------------------------------------------
  //! The sparse representation of a spin matrix
  struct MatrixSpinRep_t
  {
    int       left;              /*!< Spin row */
    int       right;             /*!< Spin column */
    ComplexD  op;                /*!< Value of data */
  };

  //----------------------------------------------------------------------------
  //! MatrixSpinRep reader
  void read(XMLReader& xml, const std::string& path, MatrixSpinRep_t& param);

  //! MatrixSpinRep writer
  void write(XMLWriter& xml, const std::string& path, const MatrixSpinRep_t& param);

  //----------------------------------------------------------------------------
  //! MatrixSpinRep reader
  void read(BinaryReader& bin, MatrixSpinRep_t& param);

  //! MatrixSpinRep writer
  void write(BinaryWriter& bin, const MatrixSpinRep_t& param);

  //----------------------------------------------------------------------------------
  //! Construct a spin matrix in the DR rep
  SpinMatrix constructSpinDR(int gamma);

  //----------------------------------------------------------------------------------
  //! Construct a spin matrix in the DP rep
  SpinMatrix constructSpinDP(int gamma);

  //----------------------------------------------------------------------------------
  //! Convert a DR gamma matrix indexed by gamma into a sparse spin representation
  std::vector<MatrixSpinRep_t> convertTwoQuarkSpinDR(int gamma);

  //----------------------------------------------------------------------------------
  //! Convert a DP gamma matrix indexed by gamma into a sparse spin representation
  std::vector<MatrixSpinRep_t> convertTwoQuarkSpinDP(int gamma);

  //----------------------------------------------------------------------------------
  //! Convert a generic spin matrix into a sparse spin representation
  std::vector<MatrixSpinRep_t> convertTwoQuarkSpin(const SpinMatrix& in);

  //----------------------------------------------------------------------------------
  //! Fold in gamma_4 for source ops
  MatrixSpinRep_t foldSourceDP(const MatrixSpinRep_t& disp, bool src);

  //! Fold in gamma_4 for source ops
  std::vector<MatrixSpinRep_t> foldSourceDP(const std::vector<MatrixSpinRep_t>& disp, bool src);


  //----------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------
  //! The sparse representation of a spin rank-3 tensor
  struct Rank3SpinRep_t
  {
    int       left;              /*!< Spin left */
    int       middle;            /*!< Spin middle */
    int       right;             /*!< Spin right */
    ComplexD  op;                /*!< Value of data */
  };

  //----------------------------------------------------------------------------
  //! Rank3SpinRep reader
  void read(XMLReader& xml, const std::string& path, Rank3SpinRep_t& param);

  //! Rank3SpinRep writer
  void write(XMLWriter& xml, const std::string& path, const Rank3SpinRep_t& param);

  //----------------------------------------------------------------------------
  //! Rank3SpinRep reader
  void read(BinaryReader& bin, Rank3SpinRep_t& param);

  //! Rank3SpinRep writer
  void write(BinaryWriter& bin, const Rank3SpinRep_t& param);

  //----------------------------------------------------------------------------------
  //! Convert a generic spin tensor into a sparse spin representation
  std::vector<Rank3SpinRep_t> convertThreeQuarkSpin(const Array3dO<Complex>& in);

  //----------------------------------------------------------------------------------
  //! Fold in gamma_4 for source ops
  Rank3SpinRep_t foldSourceDP(const Rank3SpinRep_t& disp, bool src);

  //! Fold in gamma_4 for source ops
  std::vector<Rank3SpinRep_t> foldSourceDP(const std::vector<Rank3SpinRep_t>& disp, bool src);

} // namespace Chroma

#endif
