/*! @file
 * @brief Two-flavor RatioConvConv monomial params
 */

#include "update/molecdyn/monomial/two_flavor_ratio_conv_conv_multihasen_monomial_params_w.h"
#include "io/param_io.h"

namespace Chroma
{
	// Read the parameters
	TwoFlavorRatioConvConvMultihasenWilsonTypeFermMonomialParams::
		TwoFlavorRatioConvConvMultihasenWilsonTypeFermMonomialParams(XMLReader& xml_in, const std::string& path){
			// Get the top of the parameter XML tree
			XMLReader paramtop(xml_in, path);
			try{
				read(paramtop, "Action", fermactInv);
				read(paramtop, "ShiftedMass", mu);
				read(paramtop, "NumofHasenTerms", numHasenTerms);
				if(paramtop.count("./ChronologicalPredictor") == 0){
					predictor.xml = "";
				}else{
					predictor = readXMLGroup(paramtop, "ChronologicalPredictor", "Name");
				}
			}catch(const std::string& s){
				QDPIO::cerr<<"Caught Exception while reading parameters: "<<s<<std::endl;
				QDP_abort(1);
			}
		}

	//! Read Parameters
	void read(XMLReader& xml, const std::string& path,
			TwoFlavorRatioConvConvMultihasenWilsonTypeFermMonomialParams& params){
		TwoFlavorRatioConvConvMultihasenWilsonTypeFermMonomialParams tmp(xml, path);
		params = tmp;
	}

	//! Write Parameters
	void write(XMLWriter& xml, const std::string& path,
			const TwoFlavorRatioConvConvMultihasenWilsonTypeFermMonomialParams& params){
		write(xml, "Action", params.fermactInv);
		write(xml, "ShiftedMass", params.mu);
		write(xml, "NumofHasenTerms", params.numHasenTerms);
	}

}// end namespace Chroma
