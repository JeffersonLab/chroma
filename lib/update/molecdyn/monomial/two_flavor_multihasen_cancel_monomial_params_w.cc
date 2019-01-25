/*! @file
 * @brief Two flavor Monomials - gauge action or fermion binlinear contributions for HMC
 */
#include "update/molecdyn/monomial/two_flavor_multihasen_cancel_monomial_params_w.h"

namespace Chroma
{
	// Read the parameters
	TwoFlavorMultihasenCancelMonomialParams::
		TwoFlavorMultihasenCancelMonomialParams(XMLReader& xml_in, const std::string& path)
		{
			// Get the top of the parameter XML tree
			XMLReader paramtop(xml_in, path);
			try{
				// Read the inverter Parameters
				read(paramtop, "ShiftedMass", mu);
				inv_param = readXMLGroup(paramtop, "InvertParam", "invType");
				fermact = readXMLGroup(paramtop, "FermionAction", "FermAct");

				if(paramtop.count("./ChronologicalPredictor") == 0)
				{
					predictor.xml = "";
				}else{
					predictor = readXMLGroup(paramtop, "ChronologicalPredictor", "Name");
				}
			}catch(const std::string& s){
				QDPIO::cerr<<"Caught Exception while reading parameters: "<<s<<std::endl;
				QDP_abort(1);
			}

			QDPIO::cout<<"TwoFlavorMultihasenCancelMonomialParams: read \n"<<fermact.id<<std::endl;
		}
	
	//! Read Parameters
	void read(XMLReader& xml, const std::string& path,
			TwoFlavorMultihasenCancelMonomialParams& params)
	{
		TwoFlavorMultihasenCancelMonomialParams tmp(xml, path);
		params = tmp;
	}

	//! Write Parameters
	void write(XMLWriter& xml, const std::string& path,
			const TwoFlavorMultihasenCancelMonomialParams& params)
	{
		write(xml, "ShiftedMass", params.mu);
		xml<<params.fermact.xml;
		xml<<params.inv_param.xml;
	}
}
