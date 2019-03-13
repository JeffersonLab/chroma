/*
 * syssolver_mrhs_proxy_params.cc
 *
 *  Created on: Mar 13, 2019
 *      Author: bjoo
 */

#include "actions/ferm/invert/syssolver_mrhs_proxy_params.h"

using namespace QDP;

namespace Chroma {

void read(XMLReader& xml_in, const std::string& path, SysSolverMRHSProxyParams& p)
{
	XMLReader paramtop(xml_in,path);
	read(paramtop, "BlockSize", p.BlockSize);
	p.SubInverterXML = readXMLGroup(paramtop, "SubInvertParam", "invType");
}

void write(XMLWriter& xml_out, const std::string& path, const SysSolverMRHSProxyParams& p) {
	push(xml_out, path);
	write(xml_out, "BlockSize", p.BlockSize);
	xml_out << p.SubInverterXML.xml;
	pop(xml_out);
}

SysSolverMRHSProxyParams::SysSolverMRHSProxyParams(XMLReader& xml_in, const std::string& path)
{
	read(xml_in, path,*this);
}


} // Chroma namespace


