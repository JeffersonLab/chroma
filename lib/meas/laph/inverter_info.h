#ifndef LAPH_INVERTER_INFO_H
#define LAPH_INVERTER_INFO_H

#include "chromabase.h"
#include "xml_help.h"

namespace Chroma {
  namespace LaphEnv {

// **********************************************************************
// *                                                                    *
// *  Class "InverterInfo" holds information about the inverter.        *
// *  The role of the inverters is to solve the linear system of        *
// *  equations   M*x=y  for given source y.  Conjugate gradient        *
// *  methods actually solve M^dagger M x = M^dagger y, but at the      *
// *  cost of a worse condition number.  The solution methods are       *
// *  iterative, so one must specify the maximum number of iterations   *
// *  before giving up and the relative residuum for adequate           *
// *  convergence to a solution.  Convergence is usually specified      *
// *  by   r = Mx-y  with  norm(r) <= eps^2 * norm(y).  This            *
// *  convergence criterion is applied to the **pre-conditioned**       *
// *  system (pre-conditioning is usually done, based on the            *
// *  fermion action specified).  However the top-level solver          *
// *  reports the achieved relative residuum for the original           *
// *  not-preconditioned system, which may be a little more (usually    *
// *  a factor of 2-3) than for the preconditioned system.              *
// *                                                                    *
// *  Details about the inverters can be obtained from Chroma           *
// *  documentation, but sample XML inputs are given below.             *
// *                                                                    *
// *    Conjugate-Gradient:                                             *
// *                                                                    *
// *   <InvertParam>                                                    *
// *     <invType>CG_INVERTER</invType>                                 *
// *     <RsdCG>1.0e-8</RsdCG>    <!-- Normalized Relative Residuum     *
// *     <MaxCG>10000</MaxCG>  <!-- Max no of iterations                *
// *   </InvertParam>                                                   *
// *                                                                    *
// *    Biconjugate-gradient stabilized: (does not use M^dagger M)      *
// *                                                                    *
// *   <InvertParam>                                                    *
// *     <invType>BICGSTAB_INVERTER</invType>                           *
// *     <RsdBiCGStab>1.0e-8</RsdBiCGStab>                              *
// *     <MaxBiCGStab>100000</MaxBiCGStab>                              *
// *   </InvertParam>                                                   *
// *                                                                    *
// *    Eigenvector deflated conjugate gradient                         *
// *                                                                    *
// *   <InvertParam>                                                    *
// *      <invType>EIG_CG_INVERTER</invType>                            *
// *      <RsdCG>1e-08</RsdCG>                                          *
// *      <RsdCGRestart>5e-07</RsdCGRestart>                            *
// *      <MaxCG>20000</MaxCG>                                          *
// *      <Nmax>60</Nmax>                                               *
// *      <Neig>8</Neig>                                                *
// *      <Neig_max>192</Neig_max>                                      *
// *      <updateRestartTol>0</updateRestartTol>                        *
// *      <PrintLevel>1</PrintLevel>                                    *
// *      <eigen_id>linop_evs.-0.0840</eigen_id>                        *
// *      <cleanUpEvecs>false</cleanUpEvecs>                            *
// *   </InvertParam>                                                   *
// *                                                                    *
// *                                                                    *
// *  You can also use IBICGSTAB_INVERTER  with the rest of the param   *
// *  tags the same as BiCGStab. This is an 'improved BiCGStab' with    *
// *  only 1 global reduction per iteration which may scale better      *
// *  than normal BiCGStab on architectures where global reductions     *
// *  are very expensive.                                               *
// *                                                                    *
// *  Usage:                                                            *
// *                                                                    *
// *    XMLReader xmlr(...);                                            *
// *    InverterInfo inv(xmlr);                                         *
// *              --> checks that this reader contains some valid       *
// *                  info and extracts the tolerance                   *
// *                                                                    *
// *    InverterInfo inv2(...);                                         *
// *    inv.checkEqual(inv2);                                           *
// *             -->  checks inv and inv2 have same XML content,        *
// *                  throwing string exception if not                  *
// *    inv.checkEqualTolerance(inv2);                                  *
// *             -->  checks inv and inv2 have same tolerance,          *
// *                  throwing string exception if not                  *
// *    inv.matchXMLverbatim(inv2);                                     *
// *             -->  checks inv and inv2 have exactly same XML,        *
// *                  throwing string exception if not                  *
// *                                                                    *
// *    string sval = inv.output();    <-- outputs the inverter xml     *
// *    sval = inv.output(2);          <-- indented xml output          *
// *    double val = inv.getTolerance(); <-- returns the tolerance      *
// *                                                                    *
// **********************************************************************



class InverterInfo
{

   std::string inverter_xml;
   std::string id;
   double tol;
   int max_iterations;

  public:

   InverterInfo(XMLReader& xml_in);

   InverterInfo(const InverterInfo& rhs);
                                          
   InverterInfo& operator=(const InverterInfo& rhs);
     
   ~InverterInfo(){}

   void checkEqual(const InverterInfo& rhs) const;

   void checkEqualTolerance(const InverterInfo& rhs) const;

   void matchXMLverbatim(const InverterInfo& rhs) const;

   bool operator==(const InverterInfo& rhs) const;



   std::string output(int indent = 0) const;

   void output(XMLWriter& xmlout) const;

   double getTolerance() const {return tol;}

   int getMaxIterations() const {return max_iterations;}

   std::string getId() const { return id;}


};


// *****************************************************************
  }
}
#endif
