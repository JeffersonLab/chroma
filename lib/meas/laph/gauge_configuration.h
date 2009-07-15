//This class manages the gauge configuration and its info.

#ifndef laph_gauge_configuration_h
#define laph_gauge_configuration_h

#include "meas/inline/io/named_objmap.h"
#include "handle.h"
#include "chromabase.h"

namespace Chroma
{

	namespace LaphEnv
	{

		class GaugeConfiguration
		{


			public:

				//This constructor will get the cfg from the named_obj map and
				//set the static ref to the gauge_field
				GaugeConfiguration(const std::string& gauge_id);

				const multi1d<LatticeColorMatrix>& getData() const {return *(cfg);}

				const std::string& output() const {return gauge_xml;}

				const int getTrajNum() const {return traj_num;}
			
			private:
				//Static refercence to the gauge field
				static multi1d<LatticeColorMatrix>* cfg;

				static std::string gauge_xml;

				static int traj_num;

		};

	}

}



#endif
