
//    Class "GaugeConfiguration" gives access to the actual gauge configuration
//    but also does checks when combining subsequent quantities that depend
//    on the gauge field that they came from the same configuration.
//
//      Usage:
//
//         GaugeConfiguration u(gauge_id);   
//               --> first call, use Chroma gauge id; set up static parameters.
//               --> any subsequent call to the constructor causes an abort!!
//               --> basically this is a simple singleton
//
//         u.check(gauge_xml);  
//               --> in subsequent computations, pass in the gauge_xml
//                   from a quantity (such as a quark propagator) and
//                   this function issues an abort if gauge_xml does not 
//                   match that stored in "u"
//
//         u.GetData();       <-- reference to the actual configuration
//         u.output();        <-- outputs the gauge xml description
//         u.getTrajNum();    <-- returns the integer number of the RHMC trajectory
//                                  for this configuration
//                                       


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

            std::string output() const {return gauge_xml;}

            int getTrajNum() const {return traj_num;}
         
            void check(const std::string& gauge_xml_in) const
	    {
               if (gauge_xml_in != gauge_xml)
               {
                  QDPIO::cerr << __func__ << 
                      "gauge xml does not match original XML." << endl;
                  QDP_abort(1);
               }
	    }
	    
         private:
	 
            //Static reference to the gauge field
  
            static multi1d<LatticeColorMatrix>* cfg;

            static std::string gauge_xml;

            static int traj_num;

      };

   }

}



#endif
