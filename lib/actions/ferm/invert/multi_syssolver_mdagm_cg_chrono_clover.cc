// $Id: multi_syssolver_mdagm_cg.cc,v 3.2 2006-09-20 20:28:00 edwards Exp $
/*! \file
 *  \brief Solve a MdagM*psi=chi linear system by CG2
 */

#include "actions/ferm/invert/multi_syssolver_mdagm_factory.h"
#include "actions/ferm/invert/multi_syssolver_mdagm_aggregate.h"

#include "actions/ferm/invert/multi_syssolver_mdagm_cg_chrono_clover.h"

namespace Chroma
{

#if 0
  //! CG2 system solver namespace
  namespace MdagMMultiSysSolverCGChronoEnv
  {
    //! Callback function
    MdagMMultiSystemSolver<LatticeFermion>* createFerm(XMLReader& xml_in,
						       const std::string& path,
						       Handle< LinearOperator<LatticeFermion> > A)
    {
      return new MdagMMultiSysSolverCGChrono<LatticeFermion>(A, MultiSysSolverCGParams(xml_in, path));
    }

    //! Name to be used
    const std::string name("CG_CHRONO_INVERTER");

    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::TheMdagMFermMultiSystemSolverFactory::Instance().registerObject(name, createFerm);
	registered = true;
      }
      return success;
    }
  }
#endif

  MultiSysSolverCGChronoCloverParams::MultiSysSolverCGChronoCloverParams(XMLReader& xml, 
						       const std::string& path)
  {
    XMLReader paramtop(xml, path);
    read(paramtop, "CloverParams", clovParams);
    read(paramtop, "MaxIter", MaxIter);
    read(paramtop, "MaxChrono", MaxChrono);
    read(paramtop, "Delta", Delta);
    read(paramtop, "CutoffRsd", CutoffRsd);
    read(paramtop, "RsdTarget", RsdTarget);
  }

  void read(XMLReader& xml, const std::string& path, 
	    MultiSysSolverCGChronoCloverParams& p)
  {
    MultiSysSolverCGChronoCloverParams tmp(xml, path);
    p = tmp;
  }

  void write(XMLWriter& xml, const std::string& path, 
	     const MultiSysSolverCGChronoCloverParams& p) {
    push(xml, path);
    write(xml, "CloverParams", p.clovParams);
    write(xml, "MaxIter", p.MaxIter);
    write(xml, "MaxChrono", p.MaxChrono);
    write(xml, "Delta", p.Delta);
    write(xml, "CutoffRsd", p.CutoffRsd);
    write(xml, "RsdTarget", p.RsdTarget);

    
    pop(xml);

  }

}
