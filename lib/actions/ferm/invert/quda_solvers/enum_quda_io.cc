// -*- C++ -*-
// $Id: enum_quda_io.cc,v 1.1 2009-10-06 20:34:58 bjoo Exp $
/*! \file
 *  \brief QUDA enum readers 
 */

#include "actions/ferm/invert/quda_solvers/enum_quda_io.h"
#include <string>

namespace Chroma {

  namespace  QudaSolverTypeEnv {
    bool registerAll(void)
    {
      bool success;
      success = theQudaSolverTypeMap::Instance().registerPair(string("CG"),CG);
      success &= theQudaSolverTypeMap::Instance().registerPair(string("BICGSTAB"),BICGSTAB);
      success &= theQudaSolverTypeMap::Instance().registerPair(string("GCR"),GCR);
      success &= theQudaSolverTypeMap::Instance().registerPair(string("MR"),MR);
      return success;
    }
    const string typeIDString = "QudaSolverType";
    bool regisered = registerAll();
  };

  //! Read an QudaSolverType enum
  void read(XMLReader& xml_in, const string& path, QudaSolverType& t) 
  {
    theQudaSolverTypeMap::Instance().read(QudaSolverTypeEnv::typeIDString, xml_in, path, t);
  }

  //! Write an QudaSolverType enum
  void write(XMLWriter& xml_out, const string& path, const QudaSolverType& t)
  {
    theQudaSolverTypeMap::Instance().write(QudaSolverTypeEnv::typeIDString, xml_out, path, t);
  }

  namespace  QudaPrecisionTypeEnv {
    bool registerAll(void)
    {
      bool success;
      success = theQudaPrecisionTypeMap::Instance().registerPair(string("DEFAULT"),DEFAULT);
      success &= theQudaPrecisionTypeMap::Instance().registerPair(string("HALF"),HALF);
      success &= theQudaPrecisionTypeMap::Instance().registerPair(string("SINGLE"),SINGLE);
      success &= theQudaPrecisionTypeMap::Instance().registerPair(string("DOUBLE"),DOUBLE);
      return success;
    }
    const string typeIDString = "QudaPrecisionType";
    bool regisered = registerAll();
  };


  //! Read an QudaSolverType enum
  void read(XMLReader& xml_in, const string& path, QudaPrecisionType& t) 
  {
    theQudaPrecisionTypeMap::Instance().read(QudaPrecisionTypeEnv::typeIDString, xml_in, path, t);
  }

  //! Write an QudaSolverType enum
  void write(XMLWriter& xml_out, const string& path, const QudaPrecisionType& t)
  {
    theQudaPrecisionTypeMap::Instance().write(QudaPrecisionTypeEnv::typeIDString, xml_out, path, t);
  }

  namespace  QudaReconsTypeEnv {
    bool registerAll(void)
    {
      bool success;
      success = theQudaReconsTypeMap::Instance().registerPair(string("RECONS_NONE"),RECONS_NONE);
      success &= theQudaReconsTypeMap::Instance().registerPair(string("RECONS_8"),RECONS_8);
      success &= theQudaReconsTypeMap::Instance().registerPair(string("RECONS_12"),RECONS_12);
      return success;
    }
    const string typeIDString = "QudaReconsType";
    bool regisered = registerAll();
  };

  //! Read an QudaSolverType enum
  void read(XMLReader& xml_in, const string& path, QudaReconsType& t) 
  {
    theQudaReconsTypeMap::Instance().read(QudaReconsTypeEnv::typeIDString, xml_in, path, t);
  }

  //! Write an QudaSolverType enum
  void write(XMLWriter& xml_out, const string& path, const QudaReconsType& t)
  {
    theQudaReconsTypeMap::Instance().write(QudaReconsTypeEnv::typeIDString, xml_out, path, t);
  }

  namespace  QudaSchwarzMethodEnv {
    bool registerAll(void)
    {
      bool success;
      success = theQudaSchwarzMethodMap::Instance().registerPair(string("ADDITIVE_SCHWARZ"),ADDITIVE_SCHWARZ);
      success &= theQudaSchwarzMethodMap::Instance().registerPair(string("MULTIPLICATIVE_SCHWARZ"),MULTIPLICATIVE_SCHWARZ);
      return success;
    }
    const string typeIDString = "QudaSchwarzMethod";
    bool regisered = registerAll();
  };

  //! Read an QudaSolverType enum
  void read(XMLReader& xml_in, const string& path, QudaSchwarzMethod& t) 
  {
    theQudaSchwarzMethodMap::Instance().read(QudaSchwarzMethodEnv::typeIDString, xml_in, path, t);
  }

  //! Write an QudaSolverType enum
  void write(XMLWriter& xml_out, const string& path, const QudaSchwarzMethod& t)
  {
    theQudaSchwarzMethodMap::Instance().write(QudaSchwarzMethodEnv::typeIDString, xml_out, path, t);
  }

}
