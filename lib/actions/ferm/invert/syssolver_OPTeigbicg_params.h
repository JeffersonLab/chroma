// -*- C++ -*-
// $Id: syssolver_OPTeigbicg_params.h,v 1.2 2009-10-22 20:57:26 kostas Exp $
/*! \file
 *  \brief Solve a biCG system
 */

#ifndef __syssolver_OPTeigbicg_params_h__
#define __syssolver_OPTeigbicg_params_h__

#include "chromabase.h"


namespace Chroma
{

  //! Params for EigBiCG inverter
  /*! \ingroup invert */
  struct SysSolverOptEigBiCGParams
  {
    SysSolverOptEigBiCGParams();
    SysSolverOptEigBiCGParams(XMLReader& in, const std::string& path);
    
    Real RsdCG ;           /*!< CG residual */
    int  MaxCG ;           /*!< Maximum CG iterations */
    int  PrintLevel ;      /*!< Debugg level */
    int  Neig       ;      /*!< number of eigenvectors to compute  */
    int  Nmax       ;      /*!< number of basis vectors */
    int  esize      ;      /*!< 2*lde+2*nev <= esize <= 2*(nev+1)*lde */
    int  Neig_max   ;      /*!< maximum number of eigenvectors to be refined*/

    Real restartTol ;      /*!< CG  restart tolerence: restart when 
			      |res|<restartTol*|b-A x(0)| */

    Real NormAest ; /* an estimate of the norm2 of the Matrix */

    string sort_option;   /* (IN) option for sorting eigenvalues computed by eigbicg
			     'M' for smallest magnitude and 'R' for smallest real part*/

    Real epsi;            /* (IN) threshold used to check if two eignvalues are conjugate pairs:
                                   if( imag(x)*imag(y) < 0  and
                                   abs(imag(x)+imag(y))/abs(x+y) < epsi and
                                   abs(real(x)-real(y)) / abs(x+y) < epsi )
                                   then x and y are considered to be conjugate, otherwise they are not.*/
    
    int ConvTestOpt;       /* (IN)option for how to determine convergence of the linear system
			      1  means norm(residual) < tol*norm(b)
			      2  means norm(residual) < MAX( tol*norm(b), MACHEPS*(AnormEst*norm(x)+norm(b)))
			      with MACHEPS=1e-7 for single precision version and 1e-16 for double precision version.*/
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

      Neig =0 ;
      Neig_max =0 ;
      esize = 4 ;
      NormAest = 5.0 ;

      sort_option = "M";

      epsi = 1.0e-8 ;

      ConvTestOpt = 1 ;

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
  void read(XMLReader& xml, const string& path, SysSolverOptEigBiCGParams& param);

  /*! \ingroup invert */
  void write(XMLWriter& xml, const string& path, const SysSolverOptEigBiCGParams& param);

} // End namespace

#endif 

