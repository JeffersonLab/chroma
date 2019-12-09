/*
 * syssolver_mrhs_twisted_params.h
 *
 *  Created on: Mar 14, 2019
 *      Author: bjoo
 */

#ifndef LIB_ACTIONS_FERM_INVERT_SYSSOLVER_MRHS_TWISTED_PARAMS_H_
#define LIB_ACTIONS_FERM_INVERT_SYSSOLVER_MRHS_TWISTED_PARAMS_H_

#include "chromabase.h"
#include <string>
#include "io/xml_group_reader.h"

using namespace QDP;

namespace Chroma {

	struct SysSolverMRHSTwistedParams {
		int BlockSize;   // Number of right hand sides
		multi1d<Real> Twists;
		GroupXML_t  SubInverterXML; // XML to create the sub-integrator from
		SysSolverMRHSTwistedParams() = default;
		SysSolverMRHSTwistedParams(XMLReader& xml_in, const std::string& path);
	};

	void read(XMLReader& xml, const std::string& path, SysSolverMRHSTwistedParams& p);
	void write(XMLWriter& xml, const std::string& path, const SysSolverMRHSTwistedParams& p);

}





#endif /* LIB_ACTIONS_FERM_INVERT_SYSSOLVER_MRHS_TWISTED_PARAMS_H_ */
