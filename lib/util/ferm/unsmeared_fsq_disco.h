// -*- C++ -*-
/*! \file
 * \brief tr(F^2) disconnected blocks
 */

#ifndef __unsmeared_fsq_disco_h__
#define __unsmeared_fsq_disco_h__

#include "chromabase.h"

#include <string>
#include <vector>
#include <complex>

namespace Chroma
{
  //----------------------------------------------------------------------------
  //! tr(F^2) disco operator
  struct KeyFSqDiscoOperator_t
  {
    std::string        smear;          /*!< A smearing label */
    std::vector<int>   left_lorentz;   /*!< Lorentz index pairs of left F_munu */
    std::vector<int>   right_lorentz;  /*!< Lorentz index pairs of right F_munu */
    std::vector<int>   disp_list;      /*!< Displacement dirs of right F_munu */
    multi1d<int>       mom;            /*!< D-1 momentum of this operator */
  };

  bool operator<(const KeyFSqDiscoOperator_t& a, const KeyFSqDiscoOperator_t& b);
  std::ostream& operator<<(std::ostream& os, const KeyFSqDiscoOperator_t& d);
  
  //----------------------------------------------------------------------------
  //! tr(F^2) disco operator
  class ValFSqDiscoOperator_t
  {
  public:
    std::vector< std::complex<double> > op;  
    ~ValFSqDiscoOperator_t() {}
  };

  //----------------------------------------------------------------------------
  //! KeyOperator reader    
  void read(BinaryReader& bin, KeyFSqDiscoOperator_t& d);

  //! KeyOperator writer
  void write(BinaryWriter& bin, const KeyFSqDiscoOperator_t& d);

  //----------------------------------------------------------------------------
  //!  reader
  void read(XMLReader& xml, const std::string& path, KeyFSqDiscoOperator_t& param);

  //!  writer
  void write(XMLWriter& xml, const std::string& path, const KeyFSqDiscoOperator_t& param);

  //----------------------------------------------------------------------------
  //! KeyFSqDiscoOperator_t reader
  void read(BinaryReader& bin, KeyFSqDiscoOperator_t& param);

  //! KeyHadron1PtCorr writer
  void write(BinaryWriter& bin, const KeyFSqDiscoOperator_t& param);
 
  //----------------------------------------------------------------------------
  //! ValOperator reader    
  void read(BinaryReader& bin, ValFSqDiscoOperator_t& d);

  //! ValOperator writer
  void write(BinaryWriter& bin, const ValFSqDiscoOperator_t& d);

} // namespace Hadron

#endif
