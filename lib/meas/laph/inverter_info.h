#ifndef LAPH_INVERTER_INFO_H
#define LAPH_INVERTER_INFO_H

#include "chromabase.h"
#include "xml_help.h"

namespace Chroma {
  namespace LaphEnv {

// *******************************************************************
// *                                                                 *
// *  Class "InverterInfo" holds information about the inverter.     *
// *  It also can check that the tolerance on this inverter is the   *
// *  same when compared against another InverterInfo.               *
// *  The expected XML input for this class is                       *
// *                                                                 *
// *   <InvertParam>                                                 *
// *       <RsdCG> 1e-6 </RsdCG>                                     *
// *   </InvertParam>                                                *
// *                                                                 *
// *  The tolerance of the inversion is extracted only and used      *
// *  for identification purposes.                                   *
// *                                                                 *
// *  Usage:                                                         *
// *                                                                 *
// *    XMLReader xmlr(...);                                         *
// *    InverterInfo inv(xmlr);                                      *
// *              --> checks that this reader contains some valid    *
// *                  info and extracts the tolerance                *
// *                                                                 *
// *    InverterInfo inv2(...);                                      *
// *    inv.checkEqual(inv2);                                        *
// *             -->  checks inv and inv2 have same XML content,     *
// *                  throwing string exception if not               *
// *    inv.checkEqualTolerance(inv2);                               *
// *             -->  checks inv and inv2 have same tolerance,       *
// *                  throwing string exception if not               *
// *    inv.matchXMLverbatim(inv2);                                  *
// *             -->  checks inv and inv2 have exactly same XML,     *
// *                  throwing string exception if not               *
// *                                                                 *
// *    string sval = inv.output();    <-- outputs the inverter xml  *
// *    sval = inv.output(2);          <-- indented xml output       *
// *    double val = inv.getTolerance(); <-- returns the tolerance   *
// *                                                                 *
// *******************************************************************



class InverterInfo
{

   std::string inverter_xml;
   std::string id;
   double tol;

  public:

   InverterInfo(const InverterInfo& rhs);
                                          
   InverterInfo& operator=(const InverterInfo& rhs);

   InverterInfo(XMLReader& act_rdr);
     
   ~InverterInfo(){}

   void checkEqual(const InverterInfo& rhs) const;

   void checkEqualTolerance(const InverterInfo& rhs) const;

   void matchXMLverbatim(const InverterInfo& rhs) const;

   bool operator==(const InverterInfo& rhs) const;



   std::string output(int indent = 0) const;

   double getTolerance() const {return tol;}

   std::string getId() const { return id;}


};


// *****************************************************************
  }
}
#endif
