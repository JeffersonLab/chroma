#ifndef EIGEN_IO_H
#define EIGEN_IO_H

#include "chroma.h"
#include "io/param_io.h"

#include <iostream>

using namespace std;
using namespace QDP;

// Struct for parameters needed for a Ritz KS type solve
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

// Struct for dumping the eigenvalues/vectors
struct EigenIO_t 
{
  string eigen_file;
  QDP_volfmt_t eigen_volfmt;
};


// Struct for the 
// Overall application. 
// Specialised to Wilson For now.
struct ChromaWilsonRitz_t
{
  int            version;
  Real           Mass;
  multi1d<int>   boundary;
  multi1d<int>   nrow;
  QDP::Seed      seed;
  RitzParams_t    ritz_params;
  Cfg_t          cfg;
  EigenIO_t      eigen_io_params;
};

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
#endif
