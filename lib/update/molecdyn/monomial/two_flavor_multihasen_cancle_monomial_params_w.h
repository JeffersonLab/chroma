// -*- C++ -*-
/*! @file
 * @brief Two-flavor monomial params
 * cancle term for multi-Hasenbusch prec 
 */

#ifndef __TWO_FLAVOR_MULTIHASEN_CANCLE_MONOMIAL_PARAMS_W_H_
#define __TWO_FLAVOR_MULTIHASEN_CANCLE_MONOMIAL_PARAMS_W_H_

#include "chromabase.h"
#include "update/molecdyn/monomial/comp_approx.h"
#include "io/param_io.h"
#include "io/xml_group_reader.h"

namespace Chroma
{
	struct TwoFlavorMultihasenCancleMonomialParams
	{
		TwoFlavorMultihasenCancleMonomialParams();

		TwoFlavorMultihasenCancleMonomialParams(XMLReader& in, const std::string& path);
		Real mu;	// Shifted mass 
		GroupXML_t inv_param;
		GroupXML_t fermact;
		GroupXML_t predictor;
	};

	void read(XMLReader& xml, const std::string& path, TwoFlavorMultihasenCancleMonomialParams& param);
	
	void write(XMLReader& xml, const std::string& path, const TwoFlavorMultihasenCancleMonomialParams& params);

}

#endif
