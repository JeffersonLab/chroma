#ifndef LAPH_QUARK_ACTION_INFO_H
#define LAPH_QUARK_ACTION_INFO_H

#include "chromabase.h"
#include "xml_help.h"
#include "gauge_configuration_info.h"

namespace Chroma {
  namespace LaphEnv {

// ********************************************************************
// *                                                                  *
// *  Class "QuarkInfo" holds information about the valence quark     *
// *  action.  It also can check that this action is the same as      *
// *  compared to another QuarkInfo. This will also return the quark  *
// *  mass (or kappa value).                                          *
// *                                                                  *
// *  There are three constructors for this class.                    *
// *                                                                  *
// *    XMLReader xml_in(...);                                        *
// *    QuarkInfo SQ(xml_in);                                         *
// *                                                                  *
// *      --> Looks for and extracts the <QuarkInfo> tag in the       *
// *          XML input, which must have the form                     *
// *                                                                  *
// *             <QuarkInfo>                                          *
// *                <FermionAction>                                   *
// *                   <Mass> 0.321 </Mass> (or <Kappa>..</Kappa>)    *
// *                   .....                                          *
// *                </FermionAction>                                  *
// *             </QuarkInfo>                                         *
// *                                                                  *
// *          This constructor is generally useful, but is the only   *
// *          way of specifying a valence quark action that is        *
// *          different from the action used to generate the          *
// *          configurations. Consult Chroma documentation for XML    *
// *          format of a quark action needed by the inverter.  The   *
// *          quark mass (or kappa value) is extracted from the XML.  *
// *                                                                  *
// *    XMLReader xml_in(...);                                        *
// *    GaugeConfigurationInfo U(...);                                *
// *    QuarkInfo SQ(xml_in,U);                                       *
// *                                                                  *
// *      --> The XML input expected here either has the form as for  *
// *          the constructor above, or a simplified version can be   *
// *          used:                                                   *
// *                                                                  *
// *             <QuarkInfo>                                          *
// *               <Dynamical>                                        *
// *                 <Mass> 0.332 </Mass>  (or <Kappa>...</Kappa>)    *
// *               </Dynamical>                                       *
// *             </QuarkInfo>                                         *
// *                                                                  *
// *          The XML header info for the configuration is then       *
// *          searched to extract the <FermionAction> tag with the    *
// *          appropriate mass.  This constructor is only useful if   *
// *          the desired valence quark action is the same as that    *
// *          of one of the quarks in the action used by HMC.         *
// *                                                                  *
// *    string header(...);                                           *
// *    QuarkInfo SQ(header);                                         *
// *                                                                  *
// *      --> This constructor takes a string from the header info    *
// *          of some file (such as a quark source/sink) and then     *
// *          extracts the <QuarkHeader> tag from it verbatim.        *
// *          The quark mass (or kappa value) is extracted.           *
// *                                                                  *
// *  Other usage:                                                    *
// *                                                                  *
// *    QuarkInfo S2(...);                                            *
// *    S.checkEqual(S2);  // checks S=S2, throw exception if not     *
// *    S.matchXMLverbatim(S2);                                       *
// *    if (S==S2) ...                                                *
// *                                                                  *
// *    string sval = S.output();   // xml output                     *
// *    string sval = S.output(2);  // xml indented output            *
// *    double rval = S.getMass();  // returns the quark mass         *
// *                                                                  *
// ********************************************************************                             


class QuarkInfo
{
    std::string quark_header;
    double mass;
    std::string mass_name;
    std::string action_id;

  public:

    QuarkInfo(XMLReader& xml_in);

    QuarkInfo(XMLReader& xml_in, const GaugeConfigurationInfo& U);

    QuarkInfo(const std::string& header);
           
    QuarkInfo(const QuarkInfo& rhs);

    QuarkInfo& operator=(const QuarkInfo& rhs);

    ~QuarkInfo(){}

    void checkEqual(const QuarkInfo& rhs) const;

    void matchXMLverbatim(const QuarkInfo& rhs) const;

    bool operator==(const QuarkInfo& rhs) const;



    std::string output(int indent = 0) const;

    std::string getQuarkHeader() const { return quark_header; }

    double getMass() const {return mass;}

    std::string getMassName() const {return mass_name;}
    
    std::string getActionId() const {return action_id;}


  private:

    void setMass(XMLReader& xmlr, string& massName, double& massValue);

};


// ****************************************************************
  }
}
#endif
