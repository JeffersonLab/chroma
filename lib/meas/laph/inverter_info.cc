//The definitions for the GaugeConfiguration class

#include "meas/laph/inverter_info.h"
#include "chroma.h"
#include <sstream>

namespace Chroma
{
   namespace LaphEnv
   {

		 InverterInfo::InverterInfo(XMLReader& act_rdr)
		 {

			 try{

				 XMLReader temp(act_rdr, "//InvertParam");

				 read(temp, "RsdCG", tol);

				 ostringstream strm;
				 temp.print(strm);

				 inverter_xml = strm.str();

				 QDPIO::cout << "Inverter XML = XX" << inverter_xml << "XX" << endl;

				 QDPIO::cout << "Tol = " << tol << endl;
				 QDPIO::cout << "Inverter Info Initialized" << endl;

			 }
			 catch(std::string& err)
			 {
				 QDPIO::cerr << "ERROR: " << err << ", Inverter Info invalid" << endl;
				 QDP_abort(1);
			 }

		 }
	 }

}
