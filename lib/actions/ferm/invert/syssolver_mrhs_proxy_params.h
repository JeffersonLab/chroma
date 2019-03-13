/*
 * syssolver_mrhs_proxy_params.h
 *
 *  Created on: Mar 13, 2019
 *      Author: bjoo
 */

#ifndef LIB_ACTIONS_FERM_INVERT_SYSSOLVER_MRHS_PROXY_PARAMS_H_
#define LIB_ACTIONS_FERM_INVERT_SYSSOLVER_MRHS_PROXY_PARAMS_H_

#include "chromabase.h"
#include <string>
#include "io/xml_group_reader.h"

using namespace QDP;

namespace Chroma {

	struct SysSolverMRHSProxyParams {
		int BlockSize;   // Number of right hand sides
		GroupXML_t  SubInverterXML; // XML to create the sub-integrator from
		SysSolverMRHSProxyParams(XMLReader& xml_in, const std::string& path);
	};

	void read(XMLReader& xml, const std::string& path, SysSolverMRHSProxyParams& p);
	void write(XMLWriter& xml, const std::string& path, const SysSolverMRHSProxyParams& p);

}



#endif /* LIB_ACTIONS_FERM_INVERT_SYSSOLVER_MRHS_PROXY_PARAMS_H_ */
