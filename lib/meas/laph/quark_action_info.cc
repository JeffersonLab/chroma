//The definitions for the GaugeConfiguration class

#include "meas/laph/quark_action_info.h"
#include "chroma.h"
#include <sstream>
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"

namespace Chroma
{
   namespace LaphEnv
   {

      QuarkActionInfo::QuarkActionInfo(
				XMLReader& act_rdr)
      {

				XMLReader temp(act_rdr, "//FermionAction");

				if (temp.count("Mass") != 0)
				{
					read(temp, "Mass", mass);
				}
				else
				{
					read(temp, "Kappa", mass);
				}
					
				ostringstream strm;
				temp.print(strm);

				ferm_act_xml = strm.str();

				QDPIO::cout << "Ferm Act XML = XX" << ferm_act_xml << "XX" << endl;

				std::string type;
				read(temp, "//FermAct", type);
				QDPIO::cout << "Type = " <<  type << endl;

				try
				{
					StopWatch swatch;
					swatch.reset();
					QDPIO::cout << "Try the various factories" << endl;

					// Typedefs to save typing
					typedef LatticeFermion               T;
					typedef multi1d<LatticeColorMatrix>  P;
					typedef multi1d<LatticeColorMatrix>  Q;

					//
					// Initialize fermion action
					//
					std::istringstream  xml_s(ferm_act_xml);
					XMLReader  fermacttop(xml_s);
					QDPIO::cout << "FermAct = " << type << endl;

					// Generic Wilson-Type stuff
					Handle< FermionAction<T,P,Q> >
						S_f(TheFermionActionFactory::Instance().createObject(type,
									fermacttop,
									"/FermionAction"));

				}
				catch( std::bad_cast ) 
				{
					QDPIO::cerr << __func__ << ": caught dynamic cast error" 
						<< endl;
					QDP_abort(1);
				}
				catch (const string& e) 
				{
					QDPIO::cerr << __func__ << ": map call failed: " << e 
						<< endl;
					QDP_abort(1);
				}

				QDPIO::cout << "Mass = " << mass << endl;
				QDPIO::cout << "Action Info Initialized" << endl;

			}

	 }
}

