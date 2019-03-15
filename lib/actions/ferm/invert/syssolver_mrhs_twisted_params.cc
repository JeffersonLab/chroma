/*
 * syssolver_mrhs_twisted_params.cc
 *
 *  Created on: Mar 14, 2019
 *      Author: bjoo
 */


#include "actions/ferm/invert/syssolver_mrhs_twisted_params.h"

using namespace QDP;

namespace Chroma {

void read(XMLReader& xml_in, const std::string& path, SysSolverMRHSTwistedParams& p)
{
	XMLReader paramtop(xml_in,path);
	read(paramtop, "BlockSize", p.BlockSize);
	read(paramtop, "Twists", p.Twists);
	if( p.Twists.size() != p.BlockSize ) {
		QDPIO::cout << "ERROR: SysSolverMRHSTwistedParams Twists has different size from BlockSize" << std::endl;
		QDP_abort(1);
	}
	p.SubInverterXML = readXMLGroup(paramtop, "SubInvertParam", "invType");
}

void write(XMLWriter& xml_out, const std::string& path, const SysSolverMRHSTwistedParams& p) {
	push(xml_out, path);
	write(xml_out, "BlockSize", p.BlockSize);
	write(xml_out, "Twists", p.Twists );
	xml_out << p.SubInverterXML.xml;
	pop(xml_out);
}

SysSolverMRHSTwistedParams::SysSolverMRHSTwistedParams(XMLReader& xml_in, const std::string& path)
{
	read(xml_in, path,*this);
}

}// namespace


