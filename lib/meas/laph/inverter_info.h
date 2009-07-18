
//    Class "InverterInfo" holds information about the  
//    inverter.
//    It also can check that the tolerance on this 
//    inverter is the same when compared
//    against another InverterInfo. This is a very simple 
//    wrapper over an inverter_xml, that will also return the tolerance.
//
//      Usage:
//
//         InverterInfo u(XMLReader inv_rdr);   
//               -->Checks that this reader contains some valid info, 
//                  and extracts the tolerance
//
//         u.check(InverterInfo);  
//               --> checks that another inverter has the same tolerance 
//
//         u.output();        <-- outputs the inverter xml 
//         u.getTol();    		<-- returns the (Real) tolerance 
//                                


#ifndef laph_inverter_info_h
#define laph_inverter_info_h

#include "chromabase.h"

namespace Chroma
{

   namespace LaphEnv
   {

      class InverterInfo
      {
         public:

            InverterInfo(XMLReader& act_rdr);
           
						InverterInfo( 
								const InverterInfo& rhs) : inverter_xml(rhs.output()),
								  tol(rhs.getTol()) {}
						
						InverterInfo& operator=(const InverterInfo& rhs)
						{
							inverter_xml = rhs.output();
							tol = rhs.getTol();
							
							return *this;
						}

            const std::string& output() const { return inverter_xml;}

            const Real getTol() const {return tol;}
      
						void check(const InverterInfo& rhs) const
						{
							if (toBool(rhs.getTol() != tol))
							{
								QDPIO::cerr << __func__ << 
									": Inverter tolerances do not match" << endl;
								QDP_abort(1);
							}
						}

         private:
	 
            std::string inverter_xml;

						Real tol;
      };

   }

}



#endif
