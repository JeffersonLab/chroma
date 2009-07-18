//The definitions for the GaugeConfiguration class

#include "meas/laph/gauge_configuration_info.h"
#include "chroma.h"
#include <sstream>

namespace Chroma
{
	namespace LaphEnv
	{

		GaugeConfigurationInfo::GaugeConfigurationInfo(
				XMLReader& xml_rdr)
		{

			if (xml_rdr.count("//gauge_id") == 1)
			{

				read(xml_rdr, "//gauge_id", gauge_id);

				//Test reference to the cfg in the named object map
				XMLBufferWriter gauge_xml_buff;
				try
				{   
					TheNamedObjMap::Instance().get(gauge_id).getRecordXML(gauge_xml_buff);
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


				XMLReader gauge_rdr(gauge_xml_buff);

				XMLReader temp_rdr(gauge_rdr, "//MCControl");
				ostringstream temp_strm;
				temp_rdr.print(temp_strm);
				gauge_xml = temp_strm.str();

				if (gauge_rdr.count("/Params/MCControl/StartUpdateNum") != 0)
				{
					read(gauge_rdr, "/Params/MCControl/StartUpdateNum" , traj_num);
				}
				else
				{

					QDPIO::cerr << "WARNING: No trajectory number found!" << endl;
					traj_num = 1000;
				}

			}
			else
			{

				try
				{

					XMLReader xml_temp(xml_rdr, "//MCControl");
					ostringstream strm;
					xml_temp.print(strm);

					gauge_xml = strm.str();

					traj_num = 0;
					gauge_id = "";

				}
				catch (std::string& err)
				{
					QDPIO::cout << "Invalid gauge xml" << endl;
					QDP_abort(1);
				}

			}

			QDPIO::cout << "TrajNum = " << traj_num << endl;
			QDPIO::cout << "Gauge Cfg Initialized" << endl;

		}

	}
}

