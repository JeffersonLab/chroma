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

				try {
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

				QDPIO::cout << "Mass/Kappa = " << mass << endl;
				QDPIO::cout << "Action Info Initialized" << endl;

				}
				catch(std::string& err)
				{
					QDPIO::cerr << "ERROR: " << err 
						<< ", Fermion Action Info could not be initiailized" << endl;

					QDP_abort(1);
				}
			}

	 }
}

