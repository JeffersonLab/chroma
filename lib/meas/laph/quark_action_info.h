#ifndef LAPH_QUARK_ACTION_INFO_H
#define LAPH_QUARK_ACTION_INFO_H

#include "chromabase.h"
#include "xml_help.h"

namespace Chroma {
  namespace LaphEnv {

// ********************************************************************
// *                                                                  *
// *  Class "QuarkActionInfo" holds information about the             *
// *  quark action.  It also can check that this action is the same   *
// *  when compared as another QuarkActionInfo. This is a very        *
// *  simple wrapper over a ferm_act_xml, that will also return the   *
// *  mass.                                                           *
// *                                                                  *
// *  Usage:                                                          *
// *                                                                  *
// *    XMLReader xml_in;                                             *
// *    QuarkActionInfo S(xml_in);                                    *
// *           -->Checks that this reader contains some valid info,   *
// *              and extracts the mass                               *
// *                                                                  *
// *    QuarkActionInfo S2(...);                                      *
// *    S.checkEqual(S2);  // checks S=S2, throw exception if not     *
// *    S.matchXMLverbatim(S2);                                       *
// *    if (S==S2) ...                                                *
// *                                                                  *
// *    string sval = S.output();   // xml output                     *
// *    string sval = S.output(2);  // xml indented output            *
// *    double rval = S.getMass();  // returns the quark mass         *
// *                                                                  *
// ********************************************************************                             



class QuarkActionInfo
{
    std::string ferm_act_xml;
    double mass;
    std::string mass_name;
    std::string action_id;

  public:

    QuarkActionInfo(XMLReader& act_rdr);
           
    QuarkActionInfo(const QuarkActionInfo& rhs);

    QuarkActionInfo& operator=(const QuarkActionInfo& rhs);

    ~QuarkActionInfo(){}

    void checkEqual(const QuarkActionInfo& rhs) const;

    void matchXMLverbatim(const QuarkActionInfo& rhs) const;

    bool operator==(const QuarkActionInfo& rhs) const;



    std::string output(int indent = 0) const;

    double getMass() const {return mass;}

    std::string getMassName() const {return mass_name;}
    
    std::string getActionId() const {return action_id;}


};


// ****************************************************************
  }
}
#endif
