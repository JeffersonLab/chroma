// $Id: comp_approx.cc,v 3.1 2008-05-23 21:31:32 edwards Exp $
/*! @file
 * @brief Components of rational approximation
 */

#include "update/molecdyn/monomial/comp_approx.h"

namespace Chroma 
{ 

  //-----------------------------------------------------------------------------
  //! Remez param
  void read(XMLReader& xml, const string& path, 
	    TermApprox_t& param)
  {
    XMLReader paramtop(xml, path);

    param.ratApprox = readXMLGroup(paramtop, "RationalApprox", "ratApproxType");
    param.invParam  = readXMLGroup(paramtop, "InvertParam", "invType");
  }

  //! Remez param
  void write(XMLWriter& xml, const string& path, 
	     const TermApprox_t& param)
  {
    push(xml, path);

    xml << param.ratApprox.xml;
    xml << param.invParam.xml;

    pop(xml);
  }


  //-----------------------------------------------------------------------------
  //! Remez param
  void read(XMLReader& xml, const string& path, 
	    CompApprox_t& param)
  {
    XMLReader paramtop(xml, path);

    param.fermact = readXMLGroup(paramtop, "FermionAction", "FermAct");
    read(paramtop, "ActionApprox", param.action);
    read(paramtop, "ForceApprox", param.force);
  }


  //! Write Parameters
  void write(XMLWriter& xml, const std::string& path,
	     const CompApprox_t& params) 
  {
    push(xml, path);

    xml << params.fermact.xml;
    write(xml, "ActionApprox", params.action);
    write(xml, "ForceApprox", params.force);

    pop(xml);
  }
  

  //-----------------------------------------------------------------------------
  //! Remez param
  void read(XMLReader& xml, const string& path, 
	    CompAction_t& param)
  {
    XMLReader paramtop(xml, path);

    param.fermact = readXMLGroup(paramtop, "FermionAction", "FermAct");
  }


  //! Write Parameters
  void write(XMLWriter& xml, const std::string& path,
	     const CompAction_t& params) 
  {
    push(xml, path);

    xml << params.fermact.xml;

    pop(xml);
  }
  

  //-----------------------------------------------------------------------------
  //! Remez param
  void read(XMLReader& xml, const string& path, 
	    CompActionInv_t& param)
  {
    XMLReader paramtop(xml, path);

    param.fermact  = readXMLGroup(paramtop, "FermionAction", "FermAct");

    if( paramtop.count("InvertParam") == 0) { 
      // No invert param provided.
      param.invParam.xml = "";
      param.invParam.id = "NULL";
      param.invParam.path = "";
    }
    else { 
      param.invParam = readXMLGroup(paramtop, "InvertParam", "invType");
    }

  }


  //! Write Parameters
  void write(XMLWriter& xml, const std::string& path,
	     const CompActionInv_t& params) 
  {
    push(xml, path);

    xml << params.fermact.xml;
    xml << params.invParam.xml;

    pop(xml);
  }
  
} //end namespace Chroma


