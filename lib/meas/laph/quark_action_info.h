
//    Class "QuarkActionInfo" holds information about the  
//    quark action.
//    It also can check that this action is the same when compared
//    against another QuarkActionInfo. This is a very simple 
//    wrapper over a ferm_act_xml, that will also return the mass.
//
//      Usage:
//
//         QuarkActionInfo u(XMLReader act_rdr);   
//               -->Checks that this reader contains some valid info, 
//                  and extracts the mass
//
//         u.check(QuarkActionInfo);  
//               --> checks that another action has the same xml 
//
//         u.output();        <-- outputs the action xml 
//         u.getMass();    <-- returns the (integer) mass 
//                                


#ifndef laph_quark_action_info_h
#define laph_quark_action_info_h

#include "chromabase.h"

namespace Chroma
{

   namespace LaphEnv
   {

      class QuarkActionInfo
      {
         public:

            QuarkActionInfo(XMLReader& act_rdr);
           
						QuarkActionInfo( 
								const QuarkActionInfo& rhs) : ferm_act_xml(rhs.output()),
								  mass(rhs.getMass()) {}
						
						QuarkActionInfo& operator=(const QuarkActionInfo& rhs)
						{
							ferm_act_xml = rhs.output();
							mass = rhs.getMass();
							
							return *this;
						}

            const std::string& output() const { return ferm_act_xml;}

            const Real getMass() const {return mass;}
      
						void check(const QuarkActionInfo& rhs) const
						{
							if (rhs.output() != ferm_act_xml)
							{
								QDPIO::cerr << __func__ << 
									": Fermion Action xml does not match" << endl;
								QDP_abort(1);
							}
						}

         private:
	 
            std::string ferm_act_xml;

						Real mass;
      };

   }

}



#endif
