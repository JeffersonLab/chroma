//The definitions for the GaugeConfiguration class

#include "meas/laph/gauge_configuration_handler.h"
#include "chroma.h"

namespace Chroma
{
   namespace LaphEnv
   {

      GaugeConfigurationHandler::GaugeConfigurationHandler(
					const GaugeConfigurationInfo& gauge_in) : gauge_info(gauge_in) 
      {

				std::string gauge_id = gauge_in.getGaugeId();

         if (cfg == NULL)
         {
            //Assign the cfg from the named object map
            XMLBufferWriter gauge_xml_buff;
            try
            {
               TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(gauge_id);
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
            cfg = &(TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(gauge_id));

         }
         else
         {
            QDPIO::cerr << "Sorry, only one instance of GaugeConfiguration allowed!"<<endl;
            QDP_abort(1);
	 }

      }

           // initialization of static parameters
	   
      multi1d<LatticeColorMatrix>* GaugeConfigurationHandler::cfg = NULL;

   }

}

