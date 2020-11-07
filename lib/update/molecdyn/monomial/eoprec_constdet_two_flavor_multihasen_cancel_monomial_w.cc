/*! @file
 * @brief cancel monomial for the multi-Hasenbusch term
 * with shifted linop
 */

#include "update/molecdyn/monomial/two_flavor_multihasen_cancel_monomial_w.h"
namespace Chroma
{
	namespace EvenOddPrecConstDetTwoFlavorWilsonMultihasenCancelMonomialEnv
	{
		namespace
		{
			//! Callback function for the factory
			Monomial<multi1d<LatticeColorMatrix>,
				multi1d<LatticeColorMatrix> >* createMonomial(XMLReader& xml, const std::string& path)
				{
					using T = LatticeFermion;
					using P = multi1d<LatticeColorMatrix>;
					using Q = multi1d<LatticeColorMatrix>;
					
					return new PrecConstDetTwoFlavorWilsonMultihasenCancelMonomial<T,P,Q,
						   EvenOddPrecWilsonTypeFermAct, EvenOddPrecLinearOperator>(
							TwoFlavorMultihasenCancelMonomialParams(xml, path));
				}
			bool registered = false;
		}

		const std::string name("TWO_FLAVOR_EOPREC_CONSTDET_MULTIHASEN_CANCEL_FERM_MONOMIAL");
		bool registerAll()
		{
			bool success = true;
			if(!registered)
			{
				success &= WilsonTypeFermActs4DEnv::registerAll();
				success &= TheMonomialFactory::Instance().registerObject(name, createMonomial);
				registered = true;
			}
			return success;
		}
	}
}


