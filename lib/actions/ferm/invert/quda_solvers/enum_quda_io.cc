// -*- C++ -*-
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
      success = theQudaSolverTypeMap::Instance().registerPair(std::string("CG"),CG);
      success &= theQudaSolverTypeMap::Instance().registerPair(std::string("BICGSTAB"),BICGSTAB);
      success &= theQudaSolverTypeMap::Instance().registerPair(std::string("GCR"),GCR);
      success &= theQudaSolverTypeMap::Instance().registerPair(std::string("MR"),MR);
      return success;
    }
    const std::string typeIDString = "QudaSolverType";
    bool regisered = registerAll();
  }

  //! Read an QudaSolverType enum
  void read(XMLReader& xml_in, const std::string& path, QudaSolverType& t) 
  {
    theQudaSolverTypeMap::Instance().read(QudaSolverTypeEnv::typeIDString, xml_in, path, t);
  }

  //! Write an QudaSolverType enum
  void write(XMLWriter& xml_out, const std::string& path, const QudaSolverType& t)
  {
    theQudaSolverTypeMap::Instance().write(QudaSolverTypeEnv::typeIDString, xml_out, path, t);
  }

  namespace  QudaPrecisionTypeEnv {
    bool registerAll(void)
    {
      bool success;
      success = theQudaPrecisionTypeMap::Instance().registerPair(std::string("DEFAULT"),DEFAULT);
      success &= theQudaPrecisionTypeMap::Instance().registerPair(std::string("QUARTER"),QUARTER);
      success &= theQudaPrecisionTypeMap::Instance().registerPair(std::string("HALF"),HALF);
      success &= theQudaPrecisionTypeMap::Instance().registerPair(std::string("SINGLE"),SINGLE);
      success &= theQudaPrecisionTypeMap::Instance().registerPair(std::string("DOUBLE"),DOUBLE);
      return success;
    }
    const std::string typeIDString = "QudaPrecisionType";
    bool regisered = registerAll();
  }


  //! Read an QudaSolverType enum
  void read(XMLReader& xml_in, const std::string& path, QudaPrecisionType& t) 
  {
    theQudaPrecisionTypeMap::Instance().read(QudaPrecisionTypeEnv::typeIDString, xml_in, path, t);
  }

  //! Write an QudaSolverType enum
  void write(XMLWriter& xml_out, const std::string& path, const QudaPrecisionType& t)
  {
    theQudaPrecisionTypeMap::Instance().write(QudaPrecisionTypeEnv::typeIDString, xml_out, path, t);
  }

  namespace  QudaReconsTypeEnv {
    bool registerAll(void)
    {
      bool success;
      success = theQudaReconsTypeMap::Instance().registerPair(std::string("RECONS_NONE"),RECONS_NONE);
      success &= theQudaReconsTypeMap::Instance().registerPair(std::string("RECONS_8"),RECONS_8);
      success &= theQudaReconsTypeMap::Instance().registerPair(std::string("RECONS_12"),RECONS_12);
      return success;
    }
    const std::string typeIDString = "QudaReconsType";
    bool regisered = registerAll();
  }

  //! Read an QudaSolverType enum
  void read(XMLReader& xml_in, const std::string& path, QudaReconsType& t) 
  {
    theQudaReconsTypeMap::Instance().read(QudaReconsTypeEnv::typeIDString, xml_in, path, t);
  }

  //! Write an QudaSolverType enum
  void write(XMLWriter& xml_out, const std::string& path, const QudaReconsType& t)
  {
    theQudaReconsTypeMap::Instance().write(QudaReconsTypeEnv::typeIDString, xml_out, path, t);
  }

  namespace  QudaSchwarzMethodEnv {
    bool registerAll(void)
    {
      bool success;
      success = theQudaSchwarzMethodMap::Instance().registerPair(std::string("ADDITIVE_SCHWARZ"),ADDITIVE_SCHWARZ);
      success &= theQudaSchwarzMethodMap::Instance().registerPair(std::string("MULTIPLICATIVE_SCHWARZ"),MULTIPLICATIVE_SCHWARZ);
      return success;
    }
    const std::string typeIDString = "QudaSchwarzMethod";
    bool regisered = registerAll();
  }

  //! Read an QudaSolverType enum
  void read(XMLReader& xml_in, const std::string& path, QudaSchwarzMethod& t) 
  {
    theQudaSchwarzMethodMap::Instance().read(QudaSchwarzMethodEnv::typeIDString, xml_in, path, t);
  }

  //! Write an QudaSolverType enum
  void write(XMLWriter& xml_out, const std::string& path, const QudaSchwarzMethod& t)
  {
    theQudaSchwarzMethodMap::Instance().write(QudaSchwarzMethodEnv::typeIDString, xml_out, path, t);
  }

}
