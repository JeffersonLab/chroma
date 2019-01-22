/*! @file
 * @brief Symmetric even-odd preconditioned multi-Hasenbusch monomial
 */

#include "update/molecdyn/monomial/two_flavor_ratio_conv_conv_multihasen_monomial_w.h"

namespace Chroma 
{ 
	namespace SymEvenOddPrecConstDetTwoFlavorRatioConvConvMultihasenWilsonTypeFermMonomialEnv
	{
		namespace
		{
			//! Callback function for the factory
			Monomial< multi1d<LatticeColorMatrix>,
				multi1d<LatticeColorMatrix> >* createMonomial(XMLReader& xml, const std::string& path) 
				{
					using T = LatticeFermion;
					using P = multi1d<LatticeColorMatrix>;
					using Q = multi1d<LatticeColorMatrix>;

					return new PrecConstDetTwoFlavorRatioConvConvMultihasenWilsonTypeFermMonomial<T,P,Q,
						   SymEvenOddPrecLogDetWilsonTypeFermAct, SymEvenOddPrecLogDetLinearOperator>(
							TwoFlavorRatioConvConvMultihasenWilsonTypeFermMonomialParams(xml, path));
				}

			//! Local registration flag
			bool registered = false;
		}

		const std::string name("TWO_FLAVOR_SEOPREC_CONSTDET_RATIO_CONV_CONV_MULTIHASEN_FERM_MONOMIAL");

		//! Register all the factories
		bool registerAll() 
		{
			bool success = true; 
			if (! registered)
			{
				success &= WilsonTypeFermActs4DEnv::registerAll();
				success &= TheMonomialFactory::Instance().registerObject(name, createMonomial);
				registered = true;
			}
			return success;
		}
	} //end namespace SeoPrec TwoFlavorRatioConvConvWilsonFermMonomialEnv
} //end namespace Chroma
