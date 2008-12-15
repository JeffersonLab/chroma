// -*- C++ -*-
// $Id: syssolver_OPTeigcg_params.h,v 1.4 2008-12-15 05:02:06 kostas Exp $
/*! \file
 *  \brief Solve a CG1 system
 */

#ifndef __syssolver_OPTeigcg_params_h__
#define __syssolver_OPTeigcg_params_h__

#include "chromabase.h"


namespace Chroma
{

  //! Params for EigCG inverter
  /*! \ingroup invert */
  struct SysSolverOptEigCGParams
  {
    SysSolverOptEigCGParams();
    SysSolverOptEigCGParams(XMLReader& in, const std::string& path);
    
    Real RsdCG ;           /*!< CG residual */
    int  MaxCG ;           /*!< Maximum CG iterations */
    int  PrintLevel ;      /*!< Debugg level */
    int  Neig       ;      /*!< number of eigenvectors to compute  */
    int  Nmax       ;      /*!< number of basis vectors */
    int  esize      ;      /*!< 2 <= esize <= 2*Neig + 1 */ 
    int  Neig_max   ;      /*!< maximum number of eigenvectors to be refined*/

    Real restartTol ;      /*!< CG  restart tolerence: restart when 
			      |res|<restartTol*|b-A x(0)| */

    int updateRestartTol ;  /*!< Whether to update restartTol from eresids
			     * Expensive. Requires computation of residuals
			     * =0 Never update restartTol
			     * =1 Compute all eigenresiduals and update
			     *    when ncurEvals=ldh for the first time
			     * =2 Update based on up to 10 eres picked from
			     *    ncurEvals, on every rhs that adds evecs
			     * =3 Compute all eres and update on every rhs
			     *    that adds evecs--unnecessarily expensive
			     * If updateRestartTol>0 Cholesky is not used */
    Real NormAest ; /* an estimate of the norm2 of the Matrix */
    
    bool  cleanUpEvecs ; /*!< clean up evecs upon destruction of SystemSolver*/
    string eigen_id ; /*!< named buffer holding the eigenvectors */

    struct File_t
    {
      bool read ;
      bool write ;
      std::string   file_name;
      QDP_volfmt_t  file_volfmt;
    } file;

    
    void defaults(){
      RsdCG = 1.0e-8;
      MaxCG = 1000;
      PrintLevel=2;
      restartTol = zero;
      updateRestartTol = 1;
      Neig =0 ;
      Neig_max =0 ;
      esize = 4 ;
      NormAest = 25.0 ;
      
      cleanUpEvecs=false;
      eigen_id="NULL";

      //IO control                                                                                               
      file.file_name = eigen_id ;
      file.file_volfmt = QDPIO_SINGLEFILE ;

      file.read   = false;
      file.write  = false;

    }

  };


  // Reader/writers
  /*! \ingroup invert */
  void read(XMLReader& xml, const string& path, SysSolverOptEigCGParams& param);

  /*! \ingroup invert */
  void write(XMLWriter& xml, const string& path, const SysSolverOptEigCGParams& param);

} // End namespace

#endif 

