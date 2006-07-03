// $Id: param_io.cc,v 3.1 2006-07-03 15:26:09 edwards Exp $
/*! \file
 *  \brief Various parameter readers/writers for main programs
 */

#include "chromabase.h"
#include "io/param_io.h"

namespace Chroma 
{

  //! Convert a Kappa to a mass
  Real kappaToMass(const Real& Kappa)
  {
    return 1.0/(2*Kappa) - Nd;
  }


  //! Convert a Kappa to a mass
  multi1d<Real> kappaToMass(const multi1d<Real>& Kappa)
  {
    multi1d<Real> Mass(Kappa.size());

    for(int i=0; i < Kappa.size(); ++i)
      Mass[i] = 1.0/(2*Kappa[i]) - Nd;

    return Mass;
  }


  //! Convert a Kappa to a mass
  Real massToKappa(const Real& Mass)
  {
    return 0.5/(Nd + Mass);
  }


  //! Convert a mass to a Kappa
  multi1d<Real> massToKappa(const multi1d<Real>& Mass)
  {
    multi1d<Real> Kappa(Mass.size());

    for(int i=0; i < Kappa.size(); ++i)
      Kappa[i] = 0.5/(Nd + Mass[i]);

    return Kappa;
  }


  //! Read the input version
  void read(XMLReader& xml, const string& path, IO_version_t& param)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "version", param.version);
  }


#if 0

// THIS STUFF IS OBSOLETE AND WILL BE REMOVED!!!
//
//  //! Read inverter parameters
//  void read(XMLReader& xml, const string& path, InvertParam_t& param)
//  {
//    XMLReader paramtop(xml, path);
//
//    try {
//      read(paramtop, "invType", param.invType);
//      read(paramtop, "RsdCG", param.RsdCG);
//      read(paramtop, "MaxCG", param.MaxCG);
//      param.MROver = 1;
//    
//      if( paramtop.count("RsdCGPrec") == 1 ) {
//	read(paramtop, "RsdCGPrec", param.RsdCGPrec);
//      }
//      else {
//	param.RsdCGPrec = param.RsdCG;
//      }
//
//      if( paramtop.count("MaxCGPrec") == 1 ) {
//	read(paramtop, "MaxCGPrec", param.MaxCGPrec);
//      }
//      else {
//	param.MaxCGPrec = param.MaxCG;
//      }
//
//    }
//    catch( const string& e ) { 
//      QDPIO::cerr << "Caught exception : " << e << endl;
//      QDP_abort(1);
//    }
//
//  }
//
//  //! Read inverter parameters
//  void read(XMLReader& xml, const string& path, MultiInvertParam_t& param)
//  {
//    XMLReader paramtop(xml, path);
//
//    read(paramtop, "invType", param.invType);
//    read(paramtop, "RsdCG", param.RsdCG);
//    read(paramtop, "MaxCG", param.MaxCG);
//
//    param.MROver = 1;
//  }
//
//
//  //---------------------------- Writers -----------------------------
//  //! Write inverter parameters
//  void write(XMLWriter& xml, const string& path, const InvertParam_t& param)
//  {
//    push(xml, path);
//
//    write(xml, "invType", param.invType);
//    write(xml, "RsdCG", param.RsdCG);
//    write(xml, "MaxCG", param.MaxCG);
//    write(xml, "MROver", param.MROver);
//    write(xml, "RsdCGPrec", param.RsdCGPrec);
//    write(xml, "MaxCGPrec", param.MaxCGPrec);
//    pop(xml);
//  }
//
//  //! Write inverter parameters
//  void write(XMLWriter& xml, const string& path, const MultiInvertParam_t& param)
//  {
//    push(xml, path);
//
//    write(xml, "invType", param.invType);
//    write(xml, "RsdCG", param.RsdCG);
//    write(xml, "MaxCG", param.MaxCG);
//    write(xml, "MROver", param.MROver);
//
//    pop(xml);
//  }
#endif


}  // end namespace Chroma
