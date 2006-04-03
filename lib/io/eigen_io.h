// -*- C++ -*-
// $Id: eigen_io.h,v 3.0 2006-04-03 04:58:55 edwards Exp $
/*! \file
 *  \brief Eigenvalue IO
 */

#ifndef EIGEN_IO_H
#define EIGEN_IO_H

#include "chromabase.h"
#include "io/cfgtype_io.h"
#include "io/enum_io/enum_eigenvectype_io.h"
#include "io/enum_io/enum_qdpvolfmt_io.h"

#include <iostream>

namespace Chroma 
{

  //! Struct for parameters needed for a Ritz KS type solve
  /*! \ingroup io */
  struct RitzParams_t
  {
    int  Neig;
    Real RsdR;
    Real RsdA;
    Real RsdZero;
    bool ProjApsiP;
    int  Ndummy;
    Real GammaFactor;
    int  MaxKS;
    int  MinKSIter;
    int  MaxKSIter;
    int  MaxCG;
    int  Nrenorm;
  };

  //! Struct for dumping the eigenvalues/vectors
  /*! \ingroup io */
  struct EigenIO_t 
  {
    EigenVecType eigen_filefmt;
    string eigen_file;
    QDP_volfmt_t eigen_volfmt;
  };


  //! Struct for the overall application. 
  /*! \ingroup io */
  struct ChromaWilsonRitz_t
  {
    int            version;
    std::string    fermact;
    std::string    state_info;
    multi1d<int>   nrow;
    QDP::Seed      seed;
    RitzParams_t    ritz_params;
    Cfg_t          cfg;
    EigenIO_t      eigen_io_params;
  };


  /*!
   * \ingroup io
   * @{
   */
  void read(XMLReader& xml, const string& path, RitzParams_t& header);
  void read(XMLReader& xml, const string& path, EigenIO_t& header);
  void read(XMLReader& xml, const string& path, ChromaWilsonRitz_t& header);

  void write(XMLWriter& xml, const string& path, const RitzParams_t& header);
  void write(XMLWriter& xml, const string& path, const EigenIO_t& header);
  void write(XMLWriter& xml, const string& path, const ChromaWilsonRitz_t& header);

  void writeEigen(const ChromaWilsonRitz_t& header, multi1d<Real>& lambda_lo,
		  multi1d<LatticeFermion>& eigv_lo, Real& lambda_hi,
		  QDP_serialparallel_t serpar);


  void readEigenPair(Real& lambda_lo, int& eig_index,
		     LatticeFermion& eigv, 
		     const string& filename,
		     QDP_serialparallel_t serpar,
		     XMLReader& file_xml);

  void readEigen(ChromaWilsonRitz_t& header, multi1d<Real>& lambda_lo,
		 multi1d<LatticeFermion>& eigv_lo, Real& lambda_hi,
		 const string& filename_stem, 
		 int Neig,
		 QDP_serialparallel_t serpar);

  void readEigenSzin(multi1d<Real>& lambda_lo,
		     multi1d<LatticeFermion>& eigv_lo, Real& lambda_hi,
		     const int Neig,
		     const string& filename_stem);
  /*! @} */  // end of group io
		
}  // end namespace Chroma

#endif
