/*! \file
 *  \QUDA MULTIGRID MdagM Clover solver.
 */
// comment
#include "actions/ferm/invert/syssolver_mdagm_factory.h"
#include "actions/ferm/invert/syssolver_mdagm_aggregate.h"
#include "actions/ferm/invert/quda_solvers/syssolver_quda_multigrid_clover_params.h"
#include "actions/ferm/invert/quda_solvers/syssolver_mdagm_clover_quda_multigrid_w.h"
#include "io/aniso_io.h"


#include "handle.h"
#include "actions/ferm/fermstates/periodic_fermstate.h"
#include "actions/ferm/linop/lwldslash_w.h"
#include "meas/glue/mesplq.h"
// QUDA Headers
#include <quda.h>
// #include <util_quda.h>

#include <string>
#include <iomanip>
#include <ctime>
#include <cstring>
namespace Chroma
{
  namespace MdagMSysSolverQUDAMULTIGRIDCloverEnv
  {

    //! Anonymous namespace
    namespace
    {
      //! Name to be used
      const std::string name("QUDA_MULTIGRID_CLOVER_INVERTER");

      //! Local registration flag
      bool registered = false;
    }



    MdagMSystemSolver<LatticeFermion>* createFerm(XMLReader& xml_in,	
						  const std::string& path,
						  Handle< FermState< LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > state, 
						  
						  Handle< LinearOperator<LatticeFermion> > A)
    {
      return new MdagMSysSolverQUDAMULTIGRIDClover(A, state,SysSolverQUDAMULTIGRIDCloverParams(xml_in, path));
    }

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::TheMdagMFermSystemSolverFactory::Instance().registerObject(name, createFerm);
	registered = true;
      }
      return success;
    }
  }

  SystemSolverResults_t 
  MdagMSysSolverQUDAMULTIGRIDClover::qudaInvert(const CloverTermT<T, U>& clover,
				                const CloverTermT<T, U>& invclov,
				       		const T& chi_s,
				       		T& psi_s) const{

    SystemSolverResults_t ret;

    T mod_chi;

    // Copy source into mod_chi, and zero the off-parity
    mod_chi[rb[0]] = zero;
  
  
    // This solver always solves with the SYMMETRIC preconditioned
    // Operator. If we are working with Asymmetric preconditioning 
    // Then we must apply a clover inverse.
    if( invParam.asymmetricP) { 

      //
      // symmetric
      // Solve with M_symm = 1 - A^{-1}_oo D A^{-1}ee D 
      //
      // Chroma M =  A_oo ( M_symm )
      //
      //  So  M x = b => A_oo (M_symm) x = b 
      //              =>       M_symm x = A^{-1}_oo b = chi_mod
      invclov.apply(mod_chi, chi_s, PLUS, 1);
    }
    else {
      // If we work with symmetric preconditioning nothing else needs done
      mod_chi[rb[1]] = chi_s;
    }

#ifndef BUILD_QUDA_DEVIFACE_SPINOR
    void* spinorIn =(void *)&(mod_chi.elem(rb[1].start()).elem(0).elem(0).real());
    void* spinorOut =(void *)&(psi_s.elem(rb[1].start()).elem(0).elem(0).real());
#else
    // void* spinorIn = GetMemoryPtr( mod_chi.getId() );
    // void* spinorOut = GetMemoryPtr( psi_s.getId() );
    void* spinorIn;
    void* spinorOut;
    GetMemoryPtr2(spinorIn,spinorOut,mod_chi.getId(),psi_s.getId());
#endif

    // Do the solve here 
    StopWatch swatch1; 
    swatch1.reset();
    swatch1.start();
    invertQuda(spinorOut, spinorIn, (QudaInvertParam*)&quda_inv_param);
    swatch1.stop();



    QDPIO::cout << solver_string<<"time="<< quda_inv_param.secs <<" s" ;
    QDPIO::cout << "\tPerformance="<<  quda_inv_param.gflops/quda_inv_param.secs<<" GFLOPS" ; 
    QDPIO::cout << "\tTotal Time (incl. load gauge)=" << swatch1.getTimeInSeconds() <<" s"<<std::endl;

    ret.n_count =quda_inv_param.iter;

    return ret;

 }


  void  MdagMSysSolverQUDAMULTIGRIDClover::dumpYSolver(const LatticeFermion& chi, const LatticeFermion& Y) const
  {
    // Grab the time -  took this C++ way from stackoverflow
    auto time_value = std::time(nullptr);
    auto local_time = std::localtime(&time_value);
    std::ostringstream  time_strstream;
    time_strstream << "./failed_Y_solve_" << std::put_time(local_time, "%d-%m-%Y-%H-%M-%S");
    
   
    std::string file_prefix( time_strstream.str() );

    std::string gauge_filename = file_prefix + "_gauge_field.lime";
    std::string chi_filename = file_prefix + "_chi.lime";
    std::string Y_filename = file_prefix + "_Y.lime";

    int foo = 5; // Some rubbish for the XML Files
    // Dump gauge field
    {
      XMLBufferWriter filebuf;
      XMLBufferWriter recbuf;
      push( filebuf, "GaugeFile" );
      write( filebuf, "FilePrefix", file_prefix);
      pop( filebuf);

      push( recbuf, "GaugeRecord" );
      write( recbuf, "FilePrefix", file_prefix);
      pop( recbuf );

      QDPIO::cout << "Dumping gauge links to " << gauge_filename << std::endl;

      QDPFileWriter writer(filebuf,gauge_filename, QDPIO_SINGLEFILE, QDPIO_PARALLEL);
      write(writer, recbuf, gstate->getLinks());
      writer.close();
    }
  
    // Dump chi 
    {
      XMLBufferWriter filebuf;
      XMLBufferWriter recbuf;
      push( filebuf, "ChiFile" );
      write( filebuf, "FilePrefix", file_prefix);
      pop( filebuf);

      push( recbuf, "ChiRecord" );
      write( recbuf, "FilePrefix", file_prefix);
      pop( recbuf );

      QDPIO::cout << "Dumping chi (original source) vector to " << chi_filename << std::endl;

      QDPFileWriter writer(filebuf, chi_filename, QDPIO_SINGLEFILE, QDPIO_PARALLEL);
      write(writer, recbuf, chi);
      writer.close();


    }

    // Dump Y 
    {
      XMLBufferWriter filebuf;
      XMLBufferWriter recbuf;
      push( filebuf, "YFile" );
      write( filebuf, "FilePrefix", file_prefix);
      pop( filebuf);

      push( recbuf, "YRecord" );
      write( recbuf, "FilePrefix", file_prefix);
      pop( recbuf );

      QDPIO::cout << "Dumping Y (source) vector to " << Y_filename << std::endl;

      QDPFileWriter writer(filebuf, Y_filename, QDPIO_SINGLEFILE, QDPIO_PARALLEL);
      write(writer, recbuf, Y);
      writer.close();
    }

        // Dump MG state 
    {
      auto mg_params = *(invParam.MULTIGRIDParams);
      for(int l = 0; l < mg_params.mg_levels; ++l) {
	std::ostringstream subspace_prefix;
	subspace_prefix << file_prefix << "_subspace_l" << l;

	// Up to the length of the buffer (256) padded with zeros
	std::strncpy((subspace_pointers->mg_param).vec_outfile[l], (subspace_prefix.str()).c_str(), 256);
	// If source string is too long it will be truncated and not null terminated, so null terminate
	if(  subspace_prefix.str().size() > 255 ) { (subspace_pointers->mg_param).vec_outfile[l][255] = '\0'; }
	
      }
      // Make sure everyone has thei varibles set before calling dump Multigrid
      // Strictly speaking I am using a sum as a barrier;
      double i=10;
      QDPInternal::globalSum(i);

#ifdef QUDA_MG_DUMP_ENABLED
      dumpMultigridQuda(subspace_pointers->preconditioner, &(subspace_pointers->mg_param));
#endif

      for(int l = 0; l < mg_params.mg_levels; ++l) {
	(subspace_pointers->mg_param).vec_outfile[l][0] ='\0';
      }
      QDPInternal::globalSum(i);
    }
  }


  unsigned long MdagMSysSolverQUDAMULTIGRIDClover::seqno = 0;

  void  MdagMSysSolverQUDAMULTIGRIDClover::dumpXSolver(const LatticeFermion& chi,
		  const LatticeFermion& Y,
		  const LatticeFermion& X) const

  {
    // Grab the time -  took this C++ way from stackoverflow
    auto time_value = std::time(nullptr);
    auto local_time = std::localtime(&time_value);
    std::ostringstream  time_strstream;
    time_strstream << "./failed_X_solve_" << std::put_time(local_time, "%d-%m-%Y-%H-%M-%S");
    
   
    std::string file_prefix( time_strstream.str() );

    std::string gauge_filename = file_prefix + "_gauge_field.lime";
    std::string chi_filename = file_prefix + "_chi.lime";
    std::string Y_filename = file_prefix + "_Y.lime";
    std::string X_filename = file_prefix + "_X.lime";

    int foo = 5; // Some rubbish for the XML Files
    // Dump gauge field
    {
      XMLBufferWriter filebuf;
      XMLBufferWriter recbuf;
      push( filebuf, "GaugeFile" );
      write( filebuf, "FilePrefix", file_prefix);
      pop( filebuf);

      push( recbuf, "GaugeRecord" );
      write( recbuf, "FilePrefix", file_prefix);
      pop( recbuf );

      QDPIO::cout << "Dumping gauge links to " << gauge_filename << std::endl;

      QDPFileWriter writer(filebuf,gauge_filename, QDPIO_SINGLEFILE, QDPIO_PARALLEL);
      write(writer, recbuf, gstate->getLinks());
      writer.close();
    }
  
    // Dump chi 
    {
      XMLBufferWriter filebuf;
      XMLBufferWriter recbuf;
      push( filebuf, "ChiFile" );
      write( filebuf, "FilePrefix", file_prefix);
      pop( filebuf);

      push( recbuf, "ChiRecord" );
      write( recbuf, "FilePrefix", file_prefix);
      pop( recbuf );

      QDPIO::cout << "Dumping chi (original source) vector to " << chi_filename << std::endl;

      QDPFileWriter writer(filebuf, chi_filename, QDPIO_SINGLEFILE, QDPIO_PARALLEL);
      write(writer, recbuf, chi);
      writer.close();


    }

    // Dump Y 
    {
      XMLBufferWriter filebuf;
      XMLBufferWriter recbuf;
      push( filebuf, "YFile" );
      write( filebuf, "FilePrefix", file_prefix);
      pop( filebuf);

      push( recbuf, "YRecord" );
      write( recbuf, "FilePrefix", file_prefix);
      pop( recbuf );

      QDPIO::cout << "Dumping Y (source) vector to " << Y_filename << std::endl;

      QDPFileWriter writer(filebuf, Y_filename, QDPIO_SINGLEFILE, QDPIO_PARALLEL);
      write(writer, recbuf, Y);
      writer.close();
    }

    // Dump final X
    {
      XMLBufferWriter filebuf;
      XMLBufferWriter recbuf;
      push( filebuf, "XFile" );
      write( filebuf, "FilePrefix", file_prefix);
      pop( filebuf);

      push( recbuf, "XRecord" );
      write( recbuf, "FilePrefix", file_prefix);
      pop( recbuf );

      QDPIO::cout << "Dumping X (solution) vector to " << X_filename << std::endl;

      QDPFileWriter writer(filebuf, X_filename, QDPIO_SINGLEFILE, QDPIO_PARALLEL);
      write(writer, recbuf, X);
      writer.close();
    }


    // Dump MG state 
    {
      auto mg_params = *(invParam.MULTIGRIDParams);
      for(int l = 0; l < mg_params.mg_levels; ++l) {
	std::ostringstream subspace_prefix;
	subspace_prefix << file_prefix << "_subspace_l" << l;

	// Up to the length of the buffer (256) padded with zeros
	std::strncpy((subspace_pointers->mg_param).vec_outfile[l], (subspace_prefix.str()).c_str(), 256);
	// If source string is too long it will be truncated and not null terminated, so null terminate
	if(  subspace_prefix.str().size() > 255 ) { (subspace_pointers->mg_param).vec_outfile[l][255] = '\0'; }
	
      }
      // Make sure everyone has thei varibles set before calling dump Multigrid
      // I use a global sum as a barrier here
      double i=10;
      QDPInternal::globalSum(i);
    

#ifdef QUDA_MG_DUMP_ENABLED
     dumpMultigridQuda(subspace_pointers->preconditioner, &(subspace_pointers->mg_param));
#endif

     for(int l = 0; l < mg_params.mg_levels; ++l) {
       (subspace_pointers->mg_param).vec_outfile[l][0] ='\0';
     }
     // I use a global sum as a weak barrier here.
     QDPInternal::globalSum(i);  // Make sure everyone is done
    }
  }

  

} // namespace 

