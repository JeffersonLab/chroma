/*! \file
 *  \brief Params of CG inverter
 */
#include <string>
#include "actions/ferm/invert/syssolver_fgmres_dr_params.h"
using namespace std;

namespace Chroma
{

  // Read parameters
  void read(XMLReader& xml, const std::string& path, SysSolverFGMRESDRParams& p)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "RsdTarget", p.RsdTarget);
    read(paramtop, "NKrylov",   p.NKrylov);
    read(paramtop, "NDefl",     p.NDefl);
    read(paramtop, "MaxIter",   p.MaxIter);
    p.PrecondParams = readXMLGroup(paramtop, "PrecondParams", "invType");
    
  }

  // Writer parameters
  void write(XMLWriter& xml, const std::string& path, const SysSolverFGMRESDRParams& p)
  {
    push(xml, path);
    write(xml, "RsdTarget", p.RsdTarget);
    write(xml, "NKyrlov",   p.NKrylov);
    write(xml, "NDefl",     p.NDefl);
    write(xml, "MaxIter",   p.MaxIter);
    xml << p.PrecondParams.xml;
  }

  SysSolverFGMRESDRParams::SysSolverFGMRESDRParams()
  {
    RsdTarget = 0;
    NKrylov = 0;
    NDefl = 0;
    MaxIter = 0;


    // Create a dummy XML
    XMLBufferWriter xml_buf;
    push(xml_buf, "root");
    push(xml_buf, "PrecondParams");
    write(xml_buf, "invType", "NULL");
    pop(xml_buf);
    pop(xml_buf);


    XMLReader read_back(xml_buf);
    PrecondParams = readXMLGroup(read_back, "PrecondParams", "invType");
  }

  //! Read parameters
  SysSolverFGMRESDRParams::SysSolverFGMRESDRParams(XMLReader& xml, const std::string& path)
  {
    read(xml, path, *this);
  }

}
